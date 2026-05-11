#!/usr/bin/env python3
"""
Persistent queue worker that processes GPR analysis jobs.

Polls the database for queued GPR jobs and runs up to MAX_CONCURRENT
jobs simultaneously.

Start with:
    python run_unified_worker.py

The worker writes its PID to  media/queue_worker.pid
so that services.py can detect if it's already running.
"""
import os
import sys
import time
import signal
import json
import threading

# Bootstrap Django
project_root = os.path.dirname(os.path.abspath(__file__))
if project_root not in sys.path:
    sys.path.insert(0, project_root)
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'config.settings')

import django
django.setup()

from django.conf import settings
from django.db import transaction, close_old_connections

POLL_INTERVAL = 5  # seconds between queue checks
MAX_CONCURRENT = 4  # maximum simultaneous jobs
PID_FILE = os.path.join(settings.MEDIA_ROOT, 'queue_worker.pid')
SHUTDOWN = False


def _handle_signal(signum, frame):
    global SHUTDOWN
    print(f"[worker] Received signal {signum}, shutting down after current job...")
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


# ---------------------------------------------------------------------------
# GPR job processing
# ---------------------------------------------------------------------------

def process_gpr_job_by_id(job_id):
    """Process a GPR job that has already been claimed (status set to 'processing').
    Called from a worker thread."""
    from GPR.models import GPRAnalysisJob
    from GPR.services import (
        analyze_gpr_with_complexes,
        extract_genes_from_model,
        extract_reactions_from_model,
        collect_ambiguous_cases,
        populate_blast_results_for_complexes,
        deduplicate_ambiguous_cases,
        cleanup_blast_db,
        update_model_with_gpr_rules,
        build_gpr_rules_from_unique_cases,
        build_complex_summary,
        find_all_complexes_for_uniprot,
        parse_stoichiometry_string
    )
    import pickle

    job = GPRAnalysisJob.objects.get(id=job_id)

    print(f"[worker] Processing GPR job {job.id} (organism: {job.organism_name})")

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
        print(f"[worker] GPR job {job.id}: Starting analysis...")
        gpr_results = analyze_gpr_with_complexes(
            query_genome_fasta=job.genome_fasta.path,
            complex_to_stoichiometry=complex_to_stoichiometry,
            existing_model_genes=existing_model_genes,
            evalue_threshold=1e-30,
            job=job
        )

        genome_db_path = gpr_results.get('genome_db_path')
        genome_db_dir = gpr_results.get('genome_db_dir')

        # Save BLAST results to database
        print(f"[worker] GPR job {job.id}: Saving BLAST results to database...")
        job.message = "Saving BLAST results..."
        job.save(update_fields=['message'])

        from GPR.models import GPRBlastResult
        all_blast_results = gpr_results.get('blast_results', [])
        print(f"[worker] GPR job {job.id}: Saving {len(all_blast_results)} BLAST results...")

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
                print(f"[worker] Warning saving BLAST result: {e}")

        # Save initial results
        best_hits = gpr_results.get('best_hits', {})
        job.gene_count = len(existing_model_genes)
        job.matched_genes = len(best_hits)
        job.message = "Analysis complete, checking for ambiguous cases..."
        job.save()

        # Extract reactions for context linking to genes
        print(f"[worker] GPR job {job.id}: Extracting reaction data...")
        job.message = "Extracting reaction data..."
        job.save(update_fields=['message'])
        try:
            model_reactions = extract_reactions_from_model(job.metabolic_model.path)
            print(f"[worker] GPR job {job.id}: Extracted {len(model_reactions)} reactions from model")
        except Exception as e:
            print(f"[worker] GPR job {job.id}: Warning: Could not extract reactions: {e}")
            import traceback
            traceback.print_exc()
            model_reactions = {}

        # Step 2: Detect ambiguous cases
        print(f"[worker] GPR job {job.id}: Collecting ambiguous cases...")
        job.message = "Collecting ambiguous cases..."
        job.save(update_fields=['message'])
        ambiguous_data = collect_ambiguous_cases(
            gpr_results['best_hits'],
            complex_to_stoichiometry,
            model_reactions
        )
        print(f"[worker] GPR job {job.id}: Ambiguous: {len(ambiguous_data['ambiguous_cases'])}, "
              f"Unique: {len(ambiguous_data['unique_cases'])}")

        # DEBUG: Per-gene summary
        best_hits = gpr_results.get('best_hits', {})
        print(f"\n{'='*100}")
        print(f"[DEBUG GPR] Job {job.id} — PER-GENE BLAST & COMPLEX SUMMARY")
        print(f"[DEBUG GPR] Model genes: {len(existing_model_genes)} | "
              f"Genes with BLAST hit: {len(best_hits)} | "
              f"Genes without hit: {len(existing_model_genes) - len(best_hits)}")
        print(f"{'='*100}")

        for gene_id in sorted(existing_model_genes):
            if gene_id not in best_hits:
                print(f"[DEBUG GPR]  {gene_id:30s}  |  NO BLAST HIT")
                continue

            hit = best_hits[gene_id]
            hit_uniprot = hit.get('uniprot_id', '?')
            hit_evalue = hit.get('evalue', '?')
            hit_ident = hit.get('identity', '?')
            hit_qcov = hit.get('query_coverage', '?')

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
                other_members = [
                    f"{uid}(x{cnt})" if cnt > 1 else uid
                    for uid, cnt in stoich_dict.items()
                    if uid.upper() != hit_uniprot.upper()
                ]
                print(f"[DEBUG GPR]  {gene_id:30s}  |  hit: {hit_uniprot}  "
                      f"e={hit_evalue:.1e}  id={hit_ident}%  qcov={hit_qcov}%  "
                      f"|  complex: {cpx_id}  "
                      f"|  other members: {', '.join(other_members) if other_members else '(none)'}")

        print(f"{'='*100}\n")

        # Step 3: Handle unique cases FIRST
        applied_unique_count = 0
        applied_unique_reactions = 0

        if ambiguous_data['unique_cases']:
            print(f"[worker] GPR job {job.id}: {len(ambiguous_data['unique_cases'])} unique cases")
            job.message = f"Creating GPR rules for {len(ambiguous_data['unique_cases'])} genes..."
            job.save(update_fields=['message'])

            prebuilt_gpr_rules = gpr_results.get('gpr_rules', {})
            print(f"[worker] GPR job {job.id}: Pre-built rules for {len(prebuilt_gpr_rules)} (gene, complex) pairs")
            gpr_rules = build_gpr_rules_from_unique_cases(
                ambiguous_data['unique_cases'],
                complex_to_stoichiometry,
                best_hits=best_hits,
                prebuilt_gpr_rules=prebuilt_gpr_rules
            )
            print(f"[worker] GPR job {job.id}: Generated {len(gpr_rules)} GPR rules")

            # Build complex summary
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
            print(f"[worker] GPR job {job.id}: Built {len(summary_rows)} summary rows")

            if gpr_rules:
                applied_unique_reactions = len(gpr_rules)
                try:
                    input_model_path = job.metabolic_model.path
                    _, input_ext = os.path.splitext(input_model_path)
                    if not input_ext:
                        input_ext = '.json'

                    output_model_path = os.path.join(
                        settings.MEDIA_ROOT,
                        'reconstructed_models',
                        f'model_gpr_updated_{job.id}{input_ext}'
                    )
                    os.makedirs(os.path.dirname(output_model_path), exist_ok=True)

                    result = update_model_with_gpr_rules(
                        input_model_path,
                        gpr_rules,
                        output_model_path
                    )

                    if result['success']:
                        final_rules = result.get('merged_rules', gpr_rules)
                        gpr_rules.update(final_rules)
                        job.gpr_rules_json = json.dumps(gpr_rules)
                        job.updated_model_path = output_model_path
                        applied_unique_count = result['updates']
                        print(f"[worker] GPR job {job.id}: Model updated ({applied_unique_count} reactions)")

                        # Update summary rows with merged rules
                        for row in summary_rows:
                            rxn_id = row.get('reaction_id', '')
                            if rxn_id in final_rules:
                                row['gpr_rule'] = final_rules[rxn_id]
                        job.complex_summary_json = json.dumps(summary_rows)
                    else:
                        job.gpr_rules_json = json.dumps(gpr_rules)
                        print(f"[worker] GPR job {job.id}: Model update failed: {result.get('error')}")

                except Exception as e:
                    import traceback
                    print(f"[worker] GPR job {job.id}: Error updating model: {e}")
                    traceback.print_exc()
                    job.gpr_rules_json = json.dumps(gpr_rules)

        # Step 4: Handle ambiguous cases
        if ambiguous_data['ambiguous_cases']:
            print(f"[worker] GPR job {job.id}: {len(ambiguous_data['ambiguous_cases'])} ambiguous cases")
            job.message = f"Populating BLAST results for {len(ambiguous_data['ambiguous_cases'])} ambiguous cases..."
            job.save(update_fields=['message'])

            try:
                complexes_fasta_path = os.path.join(settings.BASE_DIR, 'data', 'complexes_blast_db', 'complexes.fasta')
                populate_blast_results_for_complexes(
                    ambiguous_data,
                    genome_blast_db=genome_db_path,
                    complexes_fasta_path=complexes_fasta_path
                )
            except Exception as e:
                print(f"[worker] GPR job {job.id}: Warning - BLAST results: {e}")

            # Deduplicate: if all complexes for a gene produce the same GPR rule,
            # auto-resolve with the first complex and move to unique_cases.
            before_ambig = len(ambiguous_data['ambiguous_cases'])
            deduplicate_ambiguous_cases(ambiguous_data)
            after_ambig = len(ambiguous_data['ambiguous_cases'])
            promoted = before_ambig - after_ambig
            if promoted > 0:
                print(f"[worker] GPR job {job.id}: Deduplicated {promoted} ambiguous → unique")
                # Build rules and apply model for newly-promoted unique cases
                newly_unique = ambiguous_data['unique_cases'][-promoted:]
                prebuilt_gpr_rules = gpr_results.get('gpr_rules', {})
                extra_rules = build_gpr_rules_from_unique_cases(
                    newly_unique,
                    complex_to_stoichiometry,
                    best_hits=best_hits,
                    prebuilt_gpr_rules=prebuilt_gpr_rules
                )
                if extra_rules:
                    gpr_rules_all = {}
                    if job.gpr_rules_json:
                        try:
                            gpr_rules_all = json.loads(job.gpr_rules_json)
                        except Exception:
                            pass
                    gpr_rules_all.update(extra_rules)
                    # Re-apply all rules to model
                    try:
                        input_model_path = job.metabolic_model.path
                        _, input_ext = os.path.splitext(input_model_path)
                        if not input_ext:
                            input_ext = '.json'
                        output_model_path = os.path.join(
                            settings.MEDIA_ROOT, 'reconstructed_models',
                            f'model_gpr_updated_{job.id}{input_ext}'
                        )
                        os.makedirs(os.path.dirname(output_model_path), exist_ok=True)
                        result = update_model_with_gpr_rules(input_model_path, gpr_rules_all, output_model_path)
                        if result['success']:
                            final_rules = result.get('merged_rules', gpr_rules_all)
                            gpr_rules_all.update(final_rules)
                            job.gpr_rules_json = json.dumps(gpr_rules_all)
                            job.updated_model_path = output_model_path
                            applied_unique_count = result['updates']
                            applied_unique_reactions = len(gpr_rules_all)
                    except Exception as e:
                        print(f"[worker] GPR job {job.id}: Dedup model update error: {e}")
                    # Append to complex summary
                    complex_descriptions = {}
                    desc_path = os.path.join(settings.BASE_DIR, 'data', 'gpr_reference', 'complex_metadata.pkl')
                    try:
                        with open(desc_path, 'rb') as f:
                            complex_descriptions = pickle.load(f)
                    except Exception:
                        pass
                    extra_summary = build_complex_summary(
                        newly_unique, complex_to_stoichiometry, gpr_rules_all,
                        best_hits=best_hits, complex_descriptions=complex_descriptions,
                        prebuilt_gpr_rules=prebuilt_gpr_rules
                    )
                    existing_summary = []
                    if job.complex_summary_json:
                        try:
                            existing_summary = json.loads(job.complex_summary_json)
                        except Exception:
                            pass
                    # Update existing summary rows with merged rules too
                    for row in existing_summary:
                        rxn_id = row.get('reaction_id', '')
                        if rxn_id in gpr_rules_all:
                            row['gpr_rule'] = gpr_rules_all[rxn_id]
                    existing_summary.extend(extra_summary)
                    existing_summary.sort(key=lambda r: r.get('reaction_id', ''))
                    job.complex_summary_json = json.dumps(existing_summary)

            job.ambiguous_cases_json = json.dumps(ambiguous_data['ambiguous_cases'])

            if ambiguous_data['ambiguous_cases']:
                # Still have ambiguous cases after dedup
                try:
                    job.has_unresolved_ambiguities = True
                except AttributeError:
                    pass
                job.status = 'completed'
                if applied_unique_count > 0:
                    job.message = (
                        f"Applied {applied_unique_count} GPR rules from {applied_unique_reactions} reactions (auto). "
                        f"Review needed: {len(ambiguous_data['ambiguous_cases'])} ambiguous cases"
                    )
                else:
                    job.message = f"Ready for user review: {len(ambiguous_data['ambiguous_cases'])} cases need complex selection"
            else:
                # All ambiguous cases were resolved by deduplication
                try:
                    job.has_unresolved_ambiguities = False
                except AttributeError:
                    pass
                job.status = 'completed'
                job.message = f"Complete: Applied {applied_unique_count} GPR rules from {applied_unique_reactions} reactions (all ambiguities auto-resolved)"

        elif applied_unique_count > 0:
            try:
                job.has_unresolved_ambiguities = False
            except AttributeError:
                pass
            job.status = 'completed'
            job.message = f"Complete: Applied {applied_unique_count} GPR rules from {applied_unique_reactions} reactions"

        else:
            try:
                job.has_unresolved_ambiguities = False
            except AttributeError:
                pass
            job.status = 'completed'
            job.message = "Analysis complete: no matches found"

        job.save()
        print(f"[worker] GPR job {job.id} finished (status: {job.status})")
        return True

    except Exception as exc:
        print(f"[worker] GPR job {job.id} failed: {exc}", flush=True)
        import traceback
        traceback.print_exc()
        job.status = 'error'
        job.error_message = str(exc)
        job.save()
        return True

    finally:
        if 'genome_db_dir' in locals() and genome_db_dir:
            try:
                cleanup_blast_db(genome_db_dir)
            except Exception:
                pass


def _cleanup_stale_jobs():
    """Mark GPR jobs stuck in 'processing' for >2 hours as error."""
    from GPR.models import GPRAnalysisJob
    from django.utils import timezone
    from datetime import timedelta

    cutoff = timezone.now() - timedelta(hours=2)
    stale = GPRAnalysisJob.objects.filter(status='processing', updated_at__lt=cutoff)
    count = stale.count()
    if count:
        stale.update(status='error', error_message='Job timed out (worker crashed or stalled)')
        print(f"[worker] Cleaned up {count} stale GPRAnalysisJob(s)")


def main():
    signal.signal(signal.SIGTERM, _handle_signal)
    signal.signal(signal.SIGINT, _handle_signal)

    write_pid()
    print(f"[worker] Unified queue worker started (pid={os.getpid()}), "
          f"max {MAX_CONCURRENT} concurrent jobs, polling every {POLL_INTERVAL}s")

    cleanup_counter = 0
    job_semaphore = threading.Semaphore(MAX_CONCURRENT)
    active_threads = []

    def _run_in_thread(func):
        """Wrapper: run job function, then release semaphore slot."""
        try:
            close_old_connections()
            func()
        except Exception as exc:
            print(f"[worker] Job thread error: {exc}", flush=True)
            import traceback
            traceback.print_exc()
        finally:
            close_old_connections()
            job_semaphore.release()

    try:
        while not SHUTDOWN:
            try:
                # Periodically clean up stale jobs (every ~5 minutes)
                cleanup_counter += 1
                if cleanup_counter >= 60:  # 60 * 5s = 300s = 5 min
                    _cleanup_stale_jobs()
                    cleanup_counter = 0

                # Clean up finished threads
                active_threads[:] = [t for t in active_threads if t.is_alive()]

                # Try to acquire a concurrency slot (non-blocking)
                if not job_semaphore.acquire(blocking=False):
                    # All slots busy, wait and retry
                    time.sleep(POLL_INTERVAL)
                    continue

                # We have a slot — try to find and dispatch a job
                job_func = _claim_next_job()
                if job_func is None:
                    # No queued job, release slot and sleep
                    job_semaphore.release()
                    time.sleep(POLL_INTERVAL)
                    continue

                # Run the job in a background thread (semaphore released in wrapper)
                t = threading.Thread(target=_run_in_thread, args=(job_func,), daemon=True)
                t.start()
                active_threads.append(t)

            except Exception as exc:
                print(f"[worker] Unexpected error: {exc}", flush=True)
                import traceback
                traceback.print_exc()
                time.sleep(POLL_INTERVAL)
    finally:
        # Wait for running jobs to finish before exiting
        alive = [t for t in active_threads if t.is_alive()]
        if alive:
            print(f"[worker] Waiting for {len(alive)} active job(s) to finish...")
            for t in alive:
                t.join(timeout=300)
        remove_pid()
        print("[worker] Stopped")


def _claim_next_job():
    """Find and claim the oldest queued GPR job.
    Returns a callable that runs the job, or None if no work found.
    The job status is already set to 'processing' before returning."""
    return _claim_gpr_job()


def _claim_gpr_job():
    """Claim a GPR job and return a callable to process it, or None."""
    from GPR.models import GPRAnalysisJob

    with transaction.atomic():
        job = (GPRAnalysisJob.objects
               .select_for_update(skip_locked=True)
               .filter(status='queued')
               .order_by('created_at')
               .first())
        if job is None:
            return None
        job.status = 'processing'
        job.save()
        job_id = job.id
        job_name = job.organism_name

    print(f"[worker] Claimed GPR job {job_id} ({job_name})")

    def _run():
        process_gpr_job_by_id(job_id)

    return _run


if __name__ == '__main__':
    main()
