#!/usr/bin/env python3
"""
Persistent queue worker that processes GPR analysis jobs one at a time.

Polls the database for jobs with status='queued', picks the oldest one,
runs the analysis, detects ambiguous cases, then moves to the next.
Only ONE job runs at a time to prevent memory exhaustion on the server.

Start with:
    python manage.py shell < GPR/run_gpr_worker.py
    # or directly:
    python GPR/run_gpr_worker.py

The worker writes its PID to media/gpr_queue_worker.pid
so that services.py can detect if it's already running.
"""
import os
import sys
import time
import signal
import json

# Bootstrap Django
project_root = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'config.settings')

import django
django.setup()

from django.conf import settings
from django.db import transaction

POLL_INTERVAL = 5  # seconds between queue checks
PID_FILE = os.path.join(settings.MEDIA_ROOT, 'gpr_queue_worker.pid')
SHUTDOWN = False


def _handle_signal(signum, frame):
    global SHUTDOWN
    print(f"GPR worker received signal {signum}, shutting down after current job...")
    SHUTDOWN = True


def write_pid():
    os.makedirs(os.path.dirname(PID_FILE), exist_ok=True)
    with open(PID_FILE, 'w') as f:
        f.write(str(os.getpid()))


def remove_pid():
    try:
        os.remove(PID_FILE)
    except OSError:
        pass


def process_next_job():
    """Pick the oldest queued job and run it. Returns True if a job was processed."""
    from GPR.models import GPRAnalysisJob
    from GPR.services import (
        analyze_gpr_with_complexes,
        extract_genes_from_model,
        extract_reactions_from_model,
        collect_ambiguous_cases,
        populate_blast_results_for_complexes,
        cleanup_blast_db,
        update_model_with_gpr_rules,
        build_gpr_rules_from_unique_cases,
        build_complex_summary,
        find_all_complexes_for_uniprot,
        parse_stoichiometry_string
    )
    import pickle

    # select_for_update prevents two workers from grabbing the same job
    with transaction.atomic():
        job = (GPRAnalysisJob.objects
               .select_for_update(skip_locked=True)
               .filter(status='queued')
               .order_by('created_at')
               .first())
        if job is None:
            return False
        
        job.status = 'processing'
        job.save()

    print(f"[gpr-worker] Processing job {job.id} (organism: {job.organism_name})")
    
    try:
        # Load complex data from pickle
        complex_data_path = os.path.join(settings.BASE_DIR, 'data', 'gpr_reference', 'complex_stoichiometry.pkl')
        
        if not os.path.exists(complex_data_path):
            raise FileNotFoundError(f"Complex stoichiometry pickle not found: {complex_data_path}")
        
        with open(complex_data_path, 'rb') as f:
            complex_to_stoichiometry = pickle.load(f)
        
        # Extract genes from metabolic model
        existing_model_genes = extract_genes_from_model(job.metabolic_model.path)
        
        # Step 1: Run main GPR analysis
        print(f"[gpr-worker] Job {job.id}: Starting analysis...")
        gpr_results = analyze_gpr_with_complexes(
            query_genome_fasta=job.genome_fasta.path,
            complex_to_stoichiometry=complex_to_stoichiometry,
            existing_model_genes=existing_model_genes,
            evalue_threshold=1e-30,  # User requirement: only hits with e-value < 1e-30
            job=job  # Pass job for progress updates
        )
        
        # Keep reference to genome database until we're done
        genome_db_path = gpr_results.get('genome_db_path')
        genome_db_dir = gpr_results.get('genome_db_dir')
        
        # Save BLAST results to database
        print(f"[gpr-worker] Job {job.id}: Saving BLAST results to database...")
        job.message = "Saving BLAST results..."
        job.save(update_fields=['message'])
        
        from GPR.models import GPRBlastResult
        all_blast_results = gpr_results.get('blast_results', [])
        print(f"[gpr-worker] Job {job.id}: Saving {len(all_blast_results)} BLAST results...")
        
        for result in all_blast_results:
            try:
                GPRBlastResult.objects.create(
                    gpr_job=job,
                    gene_id=result.get('gene_id', 'unknown'),
                    complex_id=result.get('uniprot_id', 'unknown'),
                    evalue=float(result.get('evalue', 0)),
                    bitscore=float(result.get('bitscore', 0)),
                    identity=float(result.get('identity', 0)),
                    query_coverage=float(result.get('query_coverage', 0))
                )
            except Exception as e:
                print(f"[gpr-worker] Warning saving BLAST result: {e}")
        
        # Save initial results
        best_hits = gpr_results.get('best_hits', {})
        job.gene_count = len(existing_model_genes)
        job.matched_genes = len(best_hits)
        job.message = "Analysis complete, checking for ambiguous cases..."
        job.save()
        
        # Extract reactions for context linking to genes
        print(f"[gpr-worker] Job {job.id}: Extracting reaction data...")
        job.message = "Extracting reaction data..."
        job.save(update_fields=['message'])
        try:
            model_reactions = extract_reactions_from_model(job.metabolic_model.path)
            print(f"[gpr-worker] Job {job.id}: Extracted {len(model_reactions)} reactions from model")
            if model_reactions:
                # Print sample of reactions
                sample_reactions = list(model_reactions.items())[:3]
                print(f"[gpr-worker] Job {job.id}: Sample reactions: {sample_reactions}")
        except Exception as e:
            print(f"[gpr-worker] Warning: Could not extract reactions: {e}")
            import traceback
            traceback.print_exc()
            model_reactions = {}
        
        # Step 2: Detect ambiguous cases
        print(f"[gpr-worker] Job {job.id}: Collecting ambiguous cases...")
        print(f"[gpr-worker] Job {job.id}: best_hits has {len(gpr_results['best_hits'])} entries")
        print(f"[gpr-worker] Job {job.id}: model_reactions has {len(model_reactions)} entries")
        print(f"[gpr-worker] Job {job.id}: Passing model_reactions to collect_ambiguous_cases: {bool(model_reactions)}")
        
        job.message = "Collecting ambiguous cases..."
        job.save(update_fields=['message'])
        ambiguous_data = collect_ambiguous_cases(
            gpr_results['best_hits'],
            complex_to_stoichiometry,
            model_reactions
        )
        print(f"[gpr-worker] Job {job.id}: Ambiguous cases found: {len(ambiguous_data['ambiguous_cases'])}")
        print(f"[gpr-worker] Job {job.id}: Unique cases found: {len(ambiguous_data['unique_cases'])}")

        # =====================================================================
        # DEBUG: Per-gene summary — BLAST hit, complex, and complex members
        # =====================================================================
        best_hits = gpr_results.get('best_hits', {})
        print(f"\n{'='*100}")
        print(f"[DEBUG GPR] Job {job.id} — PER-GENE BLAST & COMPLEX SUMMARY")
        print(f"[DEBUG GPR] Model genes: {len(existing_model_genes)} | "
              f"Genes with BLAST hit: {len(best_hits)} | "
              f"Genes without hit: {len(existing_model_genes) - len(best_hits)}")
        print(f"{'='*100}")

        # Iterate over ALL model genes (sorted for readability)
        for gene_id in sorted(existing_model_genes):
            if gene_id not in best_hits:
                print(f"[DEBUG GPR]  {gene_id:30s}  |  NO BLAST HIT")
                continue

            hit = best_hits[gene_id]
            hit_uniprot = hit.get('uniprot_id', '?')
            hit_evalue = hit.get('evalue', '?')
            hit_ident = hit.get('identity', '?')
            hit_qcov = hit.get('query_coverage', '?')

            # Find complexes containing this UniProt
            complexes_for_gene = find_all_complexes_for_uniprot(
                hit_uniprot, complex_to_stoichiometry
            )

            if not complexes_for_gene:
                print(f"[DEBUG GPR]  {gene_id:30s}  |  hit: {hit_uniprot}  "
                      f"e={hit_evalue:.1e}  id={hit_ident}%  qcov={hit_qcov}%  "
                      f"|  NO COMPLEX FOUND")
                continue

            for cpx_id in complexes_for_gene:
                stoich_str = complex_to_stoichiometry.get(cpx_id, '')
                stoich_dict = parse_stoichiometry_string(stoich_str)
                # Other members = all UniProt IDs in the complex except the hit itself
                other_members = [
                    f"{uid}(x{cnt})" if cnt > 1 else uid
                    for uid, cnt in stoich_dict.items()
                    if uid.upper() != hit_uniprot.upper()
                ]
                print(f"[DEBUG GPR]  {gene_id:30s}  |  hit: {hit_uniprot}  "
                      f"e={hit_evalue:.1e}  id={hit_ident}%  qcov={hit_qcov}%  "
                      f"|  complex: {cpx_id}  "
                      f"|  other members: {', '.join(other_members) if other_members else '(none — single subunit)'}")


        print(f"{'='*100}\n")
        
        # Step 1: Handle unique cases FIRST (auto-apply regardless of ambiguous cases)
        applied_unique_count = 0
        applied_unique_reactions = 0
        
        if ambiguous_data['unique_cases']:
            print(f"[gpr-worker] Job {job.id}: Found {len(ambiguous_data['unique_cases'])} unique cases (single complex)")
            job.message = f"Creating GPR rules for {len(ambiguous_data['unique_cases'])} genes with single complex match..."
            job.save(update_fields=['message'])
            
            # Build GPR rules using pre-built rules from DIAMOND complex member search
            print(f"[gpr-worker] Job {job.id}: Building GPR rules from unique cases...")
            prebuilt_gpr_rules = gpr_results.get('gpr_rules', {})
            print(f"[gpr-worker] Job {job.id}: Pre-built GPR rules available for {len(prebuilt_gpr_rules)} (gene, complex) pairs")
            gpr_rules = build_gpr_rules_from_unique_cases(
                ambiguous_data['unique_cases'],
                complex_to_stoichiometry,
                best_hits=best_hits,
                prebuilt_gpr_rules=prebuilt_gpr_rules
            )
            print(f"[gpr-worker] Job {job.id}: Generated {len(gpr_rules)} GPR rules from {len(ambiguous_data['unique_cases'])} unique cases")
            
            # Build and store complex summary for display
            complex_descriptions = {}
            desc_path = os.path.join(settings.BASE_DIR, 'data', 'gpr_reference', 'complex_metadata.pkl')
            try:
                with open(desc_path, 'rb') as f:
                    complex_descriptions = pickle.load(f)
            except Exception:
                pass
            summary_rows = build_complex_summary(
                ambiguous_data['unique_cases'],
                complex_to_stoichiometry,
                gpr_rules,
                best_hits=best_hits,
                complex_descriptions=complex_descriptions,
                prebuilt_gpr_rules=prebuilt_gpr_rules
            )
            job.complex_summary_json = json.dumps(summary_rows)
            print(f"[gpr-worker] Job {job.id}: Built {len(summary_rows)} summary rows")
            
            if gpr_rules:
                print(f"[gpr-worker] Job {job.id}: Created {len(gpr_rules)} GPR rules")
                applied_unique_reactions = len(gpr_rules)
                
                # Try to create updated model with same format as input
                try:
                    # Get the input model's file extension
                    input_model_path = job.metabolic_model.path
                    print(f"[gpr-worker] Job {job.id}: Input model path: {input_model_path}")
                    _, input_ext = os.path.splitext(input_model_path)
                    if not input_ext:
                        input_ext = '.json'
                    print(f"[gpr-worker] Job {job.id}: Input model extension: {input_ext}")
                    
                    output_model_path = os.path.join(
                        settings.MEDIA_ROOT,
                        'reconstructed_models',
                        f'model_gpr_updated_{job.id}{input_ext}'
                    )
                    os.makedirs(os.path.dirname(output_model_path), exist_ok=True)
                    print(f"[gpr-worker] Job {job.id}: Output model path: {output_model_path}")
                    
                    result = update_model_with_gpr_rules(
                        input_model_path,
                        gpr_rules,
                        output_model_path
                    )
                    
                    print(f"[gpr-worker] Job {job.id}: Update result: {result}")
                    
                    if result['success']:
                        # Store path and rules for download
                        job.gpr_rules_json = json.dumps(gpr_rules)
                        job.updated_model_path = output_model_path
                        applied_unique_count = result['updates']
                        print(f"[gpr-worker] Job {job.id}: Updated model saved to {output_model_path} ({applied_unique_count} reactions)")
                        print(f"[gpr-worker] Job {job.id}: Saved updated_model_path to DB: {job.updated_model_path}")
                    else:
                        job.gpr_rules_json = json.dumps(gpr_rules)
                        print(f"[gpr-worker] Job {job.id}: Model update failed: {result.get('error')}")
                        
                except Exception as e:
                    import traceback
                    print(f"[gpr-worker] Job {job.id}: Error updating model: {e}")
                    traceback.print_exc()
                    job.gpr_rules_json = json.dumps(gpr_rules)
            else:
                print(f"[gpr-worker] Job {job.id}: No GPR rules generated from {len(ambiguous_data['unique_cases'])} unique cases")
        
        # Step 2: Handle ambiguous cases (requires user review)
        if ambiguous_data['ambiguous_cases']:
            print(f"[gpr-worker] Job {job.id}: Found {len(ambiguous_data['ambiguous_cases'])} ambiguous cases")
            job.message = f"Found {len(ambiguous_data['ambiguous_cases'])} ambiguous cases, populating BLAST results..."
            job.save(update_fields=['message'])
            
            # Populate BLAST results for ambiguous cases
            print(f"[gpr-worker] Job {job.id}: Populating BLAST results for complexes...")
            try:
                complexes_fasta_path = os.path.join(settings.BASE_DIR, 'data', 'complexes_blast_db', 'complexes.fasta')
                populate_blast_results_for_complexes(
                    ambiguous_data,
                    genome_blast_db=genome_db_path,
                    complexes_fasta_path=complexes_fasta_path,
                    existing_genes=existing_model_genes
                )
            except Exception as e:
                print(f"[gpr-worker] Job {job.id}: Warning - could not populate all BLAST results: {e}")
            
            # Store ambiguous cases and mark for user review
            job.ambiguous_cases_json = json.dumps(ambiguous_data['ambiguous_cases'])
            
            # Determine final status - ALWAYS set to 'completed' to not block other jobs
            # Use has_unresolved_ambiguities flag to track whether user review is needed
            if applied_unique_count > 0:
                # Some rules applied, some need review
                try:
                    job.has_unresolved_ambiguities = True
                except AttributeError:
                    pass  # Column not in DB yet (migration pending)
                job.status = 'completed'  # Job completes even with ambiguous cases
                job.message = (
                    f"Applied {applied_unique_count} GPR rules from {applied_unique_reactions} reactions (auto). "
                    f"Review needed: {len(ambiguous_data['ambiguous_cases'])} ambiguous cases"
                )
            else:
                # No automatic rules applied, only ambiguous cases
                try:
                    job.has_unresolved_ambiguities = True
                except AttributeError:
                    pass  # Column not in DB yet (migration pending)
                job.status = 'completed'  # Job completes - user can resolve ambiguities asynchronously
                job.message = f"Ready for user review: {len(ambiguous_data['ambiguous_cases'])} cases need complex selection"
        
        elif applied_unique_count > 0:
            # Only unique cases, already applied
            try:
                job.has_unresolved_ambiguities = False
            except AttributeError:
                pass  # Column not in DB yet (migration pending)
            job.status = 'completed'
            job.message = f"Complete: Applied {applied_unique_count} GPR rules from {applied_unique_reactions} reactions"
        
        else:
            print(f"[gpr-worker] Job {job.id}: No cases found (ambiguous: {len(ambiguous_data['ambiguous_cases'])}, unique: {len(ambiguous_data['unique_cases'])})")
            try:
                job.has_unresolved_ambiguities = False
            except AttributeError:
                pass  # Column not in DB yet (migration pending)
            job.status = 'completed'
            job.message = "Analysis complete: no matches found"
        
        job.save()
        print(f"[gpr-worker] Job {job.id} finished (status: {job.status})")
        return True
        
    except Exception as exc:
        print(f"[gpr-worker] Job {job.id} failed: {exc}", flush=True)
        import traceback
        traceback.print_exc()
        job.status = 'error'
        job.error_message = str(exc)
        job.save()
        return True
    
    finally:
        # Cleanup genome database
        if 'genome_db_dir' in locals() and genome_db_dir:
            try:
                cleanup_blast_db(genome_db_dir)
            except:
                pass


def main():
    signal.signal(signal.SIGTERM, _handle_signal)
    signal.signal(signal.SIGINT, _handle_signal)

    write_pid()
    print(f"[gpr-worker] Started (pid={os.getpid()}), polling every {POLL_INTERVAL}s")

    try:
        while not SHUTDOWN:
            try:
                had_work = process_next_job()
            except Exception as exc:
                print(f"[gpr-worker] Unexpected error: {exc}", flush=True)
                import traceback
                traceback.print_exc()
                had_work = False
            
            if not had_work:
                time.sleep(POLL_INTERVAL)
    finally:
        remove_pid()
        print("[gpr-worker] Stopped")


if __name__ == '__main__':
    main()
