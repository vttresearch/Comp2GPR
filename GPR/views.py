from django.shortcuts import render, redirect, get_object_or_404
from django.http import JsonResponse
from django.views.decorators.http import require_http_methods
from django.views.decorators.csrf import ensure_csrf_cookie
from django.contrib.auth.decorators import login_required
from django.utils import timezone
from .models import GPRAnalysisJob, GPRBlastResult, ComplexMetadata, GPRAmbiguousCaseResolution
from .forms import GPRAnalysisForm, ComplexSelectionForm
from .services import (
    extract_genes_from_model, 
    extract_protein_sequences,
    run_blast_search,
    generate_gpr_summary,
    parse_stoichiometry_string,
    normalize_uniprot_id,
    update_model_with_gpr_rules,
    build_complex_summary
)
import uuid
import json
import os
from django.http import Http404


def _get_own_gpr_job_or_404(request, job_id):
    """Return the GPRAnalysisJob only if it belongs to this user/session."""
    job = get_object_or_404(GPRAnalysisJob, id=job_id)
    if request.user.is_authenticated and job.user_id and job.user_id == request.user.id:
        return job
    sk = request.session.session_key
    if sk and job.session_key == sk:
        return job
    raise Http404


def gpr_analysis_index(request):
    """Main GPR analysis page - form available to all users"""
    jobs = []
    form = None
    
    # Show form to all users
    if request.method == 'POST':
        form = GPRAnalysisForm(request.POST, request.FILES)
        if form.is_valid():
            job = form.save(commit=False)
            job.job_id = str(uuid.uuid4())[:8]
            job.organism_name = form.cleaned_data['organism_name']
            job.status = 'queued'
            # Only set user if authenticated
            if request.user.is_authenticated:
                job.user = request.user
            if not request.session.session_key:
                request.session.create()
            job.session_key = request.session.session_key
            job.save()
            
            return redirect('gpr:job_detail', job_id=job.id)
    else:
        form = GPRAnalysisForm()
    
    # Show jobs belonging to current user or session
    if request.user.is_authenticated:
        jobs = GPRAnalysisJob.objects.filter(user=request.user).order_by('-created_at')[:10]
    elif request.session.session_key:
        jobs = GPRAnalysisJob.objects.filter(session_key=request.session.session_key).order_by('-created_at')[:10]
    
    context = {
        'form': form,
        'jobs': jobs,
        'page_title': 'GPR Analysis'
    }
    return render(request, 'GPR/index.html', context)


@ensure_csrf_cookie
def gpr_job_detail(request, job_id):
    """Display GPR analysis job details"""
    job = _get_own_gpr_job_or_404(request, job_id)
    all_blast_results = job.blast_results.all().order_by('gene_id', 'evalue')
    
    # Keep only the first (best) hit per gene
    seen_genes = set()
    blast_results_filtered = []
    for result in all_blast_results:
        if result.gene_id not in seen_genes:
            seen_genes.add(result.gene_id)
            blast_results_filtered.append(result)
    
    # Compute match rate
    if job.gene_count and job.gene_count > 0:
        match_rate = round(job.matched_genes / job.gene_count * 100, 1)
    else:
        match_rate = None

    # Parse complex summary if available
    complex_summary = []
    if job.complex_summary_json:
        try:
            complex_summary = json.loads(job.complex_summary_json)
        except Exception:
            pass

    # Enrich summary rows with reaction formulas (real compound names)
    if complex_summary:
        from .services import extract_reactions_from_model
        try:
            model_reactions = extract_reactions_from_model(job.metabolic_model.path)
            for row in complex_summary:
                rxn_data = model_reactions.get(row.get('reaction_id', ''), {})
                row['formula'] = rxn_data.get('formula', '')
                row['original_gpr'] = rxn_data.get('gpr', '')
        except Exception:
            pass

    # Queue position and estimated wait
    pos = job.queue_position
    queue_len = GPRAnalysisJob.queue_length()
    estimated_wait_min = pos * 25 if pos else 0  # ~25 min per job avg
    estimated_wait_max = pos * 30 if pos else 0

    context = {
        'job': job,
        'blast_results': blast_results_filtered,
        'complex_summary': complex_summary,
        'match_rate': match_rate,
        'evalue_threshold': 1e-30,
        'queue_position': pos,
        'queue_length': queue_len,
        'estimated_wait_min': estimated_wait_min,
        'estimated_wait_max': estimated_wait_max,
        'page_title': f'GPR Analysis - {job.organism_name}'
    }
    return render(request, 'GPR/job_detail.html', context)


def gpr_job_status(request, job_id):
    """API endpoint for job status updates"""
    job = _get_own_gpr_job_or_404(request, job_id)
    
    pos = job.queue_position
    return JsonResponse({
        'id': job.id,
        'job_id': job.job_id,
        'status': job.status,
        'organism_name': job.organism_name,
        'gene_count': job.gene_count or 0,
        'matched_genes': job.matched_genes or 0,
        'message': job.message or '',
        'error_message': job.error_message,
        'created_at': job.created_at.isoformat(),
        'updated_at': job.updated_at.isoformat(),
        'queue_position': pos,
        'queue_length': GPRAnalysisJob.queue_length(),
        'estimated_wait_min': pos * 25 if pos else 0,
        'estimated_wait_max': pos * 30 if pos else 0,
    })


def gpr_results_download(request, job_id):
    """Download GPR analysis results or updated model"""
    job = _get_own_gpr_job_or_404(request, job_id)
    
    # Prefer updated model if available
    if job.updated_model_path and os.path.exists(job.updated_model_path):
        from django.http import FileResponse
        # Get the file extension from the actual file
        _, file_ext = os.path.splitext(job.updated_model_path)
        filename = f"model_with_gpr_rules_{job.job_id}{file_ext}"
        response = FileResponse(open(job.updated_model_path, 'rb'), as_attachment=True, filename=filename)
        return response
    
    # Fallback to results file
    if job.results_file:
        from django.http import FileResponse
        return FileResponse(job.results_file.open('rb'), as_attachment=True)
    
    return JsonResponse({'error': 'Results not available'}, status=404)


def gpr_resolve_ambiguities(request, job_id):
    """Display ambiguous complex assignments for user resolution"""
    job = _get_own_gpr_job_or_404(request, job_id)
    
    # Check if there are unresolved ambiguities (job can be 'completed' but still have ambiguities)
    # Handle case where column might not exist yet (migration pending)
    try:
        has_ambiguities = job.has_unresolved_ambiguities
    except AttributeError:
        has_ambiguities = False  # Default to False if column doesn't exist
    
    if not has_ambiguities:
        return redirect('gpr:job_detail', job_id=job_id)
    
    # Parse stored ambiguous cases
    ambiguous_data = json.loads(job.ambiguous_cases_json) if job.ambiguous_cases_json else []
    
    if not ambiguous_data:
        job.status = 'completed'
        job.save()
        return redirect('gpr:job_detail', job_id=job_id)
    
    # Load complex descriptions from pickle file
    import pickle
    from django.conf import settings
    complex_descriptions = {}
    pkl_path = os.path.join(settings.BASE_DIR, 'data', 'gpr_reference', 'complex_metadata.pkl')
    try:
        with open(pkl_path, 'rb') as f:
            complex_descriptions = pickle.load(f)  # {complex_id: description}
    except Exception as e:
        print(f"[resolve_ambiguities] Warning: could not load complex_metadata.pkl: {e}")
    
    # Get complex metadata for all involved complexes
    all_complex_ids = set()
    case_details = []
    
    # First pass: build case_details as before (for backward compatibility with POST processing)
    for case in ambiguous_data:
        case_details.append({
            'genome_gene': case['genome_gene'],
            'uniprot_id': case['uniprot_id'],
            'blast_hit': case['blast_hit'],
            'reactions': case.get('reactions_with_gene', []),
            'complexes': case['possible_complexes']
        })
        for complex_info in case['possible_complexes']:
            all_complex_ids.add(complex_info['complex_id'])
    
    # REORGANIZE BY REACTION for template display
    # Group all genes by reaction, preserving indices for form fields
    reactions_with_genes = {}  # {reaction_id: {reaction_data, genes: [{gene_info_with_indices}, ...]}}

    # Extract reaction formulas (with real compound names) from the model
    from .services import extract_reactions_from_model
    try:
        model_reactions = extract_reactions_from_model(job.metabolic_model.path)
    except Exception as e:
        print(f"[resolve_ambiguities] Warning: could not extract reactions: {e}")
        model_reactions = {}

    for case_idx, case in enumerate(case_details):
        genome_gene = case['genome_gene']
        for reaction_idx, reaction in enumerate(case.get('reactions', [])):
            reaction_id = reaction['reaction_id']
            if reaction_id not in reactions_with_genes:
                # Look up formula from extracted model reactions
                rxn_data = model_reactions.get(reaction_id, {})
                reactions_with_genes[reaction_id] = {
                    'reaction_id': reaction_id,
                    'reaction_name': reaction.get('name', 'Unknown'),
                    'gpr': reaction.get('gpr', ''),
                    'formula': rxn_data.get('formula', ''),
                    'genes': []
                }
            
            # Add this gene's case to the reaction with indices for form field naming
            reactions_with_genes[reaction_id]['genes'].append({
                'genome_gene': genome_gene,
                'uniprot_id': case['uniprot_id'],
                'blast_hit': case['blast_hit'],
                'complexes': case['complexes'],  # All possible complexes for this gene
                'case_idx': case_idx,            # Index in case_details (for form field names)
                'reaction_idx': reaction_idx     # Index in case['reactions'] (for form field names)
            })
    
    # Convert to list and sort
    reactions_list = list(reactions_with_genes.values())
    
    # Merge DB metadata + pickle descriptions
    complex_metadata_dict = {}
    for complex_id in all_complex_ids:
        meta = {
            'name': complex_id,
            'aliases': [],
            'description': complex_descriptions.get(complex_id, '')
        }
        try:
            db_meta = ComplexMetadata.objects.get(complex_id=complex_id)
            meta['name'] = db_meta.name
            meta['aliases'] = db_meta.get_aliases_list()
            if db_meta.description and not meta['description']:
                meta['description'] = db_meta.description
        except ComplexMetadata.DoesNotExist:
            pass
        complex_metadata_dict[complex_id] = meta
    
    # Inject descriptions into each case's complex list so the template can access them
    for case in case_details:
        for cpx in case['complexes']:
            cpx_id = cpx['complex_id']
            if cpx_id in complex_metadata_dict:
                cpx['description'] = complex_metadata_dict[cpx_id].get('description', '')
                cpx['display_name'] = complex_metadata_dict[cpx_id].get('name', cpx_id)
                cpx['aliases'] = complex_metadata_dict[cpx_id].get('aliases', [])
    
    # Load stoichiometry data for each complex (if not already in complex_info)
    from .services import parse_stoichiometry_string
    complex_to_stoichiometry = {}
    stoich_path = os.path.join(settings.BASE_DIR, 'data', 'gpr_reference', 'complex_stoichiometry.pkl')
    try:
        with open(stoich_path, 'rb') as f:
            complex_to_stoichiometry = pickle.load(f)
        print(f"[resolve_ambiguities] Loaded stoichiometry for {len(complex_to_stoichiometry)} complexes")
    except Exception as e:
        print(f"[resolve_ambiguities] Warning: could not load complex_stoichiometry.pkl: {e}")
    
    # Ensure stoichiometry is in the right format for the template
    for case in case_details:
        for cpx in case['complexes']:
            cpx_id = cpx['complex_id']
            
            # Add stoichiometry if not already there
            if 'stoichiometry' not in cpx or not cpx['stoichiometry']:
                if cpx_id in complex_to_stoichiometry:
                    stoich_str = complex_to_stoichiometry[cpx_id]
                    cpx['stoichiometry'] = parse_stoichiometry_string(stoich_str)
                    print(f"[resolve_ambiguities] Complex {cpx_id}: loaded stoichiometry with {len(cpx['stoichiometry'])} proteins")
                else:
                    cpx['stoichiometry'] = {}
                    print(f"[resolve_ambiguities] Complex {cpx_id}: NOT FOUND in stoichiometry pickle!")
            
            # Ensure blast_results exists (should already be there from worker)
            if 'blast_results' not in cpx:
                cpx['blast_results'] = []
                print(f"[resolve_ambiguities] Complex {cpx_id}: no blast_results, creating empty list")
            else:
                print(f"[resolve_ambiguities] Complex {cpx_id}: has {len(cpx.get('blast_results', []))} blast results")
    
    if request.method == 'POST':
        # Process user selections (checkboxes → multiple values per field)
        selections = {}  # {(gene, reaction_id): [complex_id, ...]}
        for case_idx, case in enumerate(case_details):
            gene = case['genome_gene']
            
            # Get selected complex(es) for each reaction
            for reaction_idx, reaction in enumerate(case.get('reactions', [])):
                reaction_id = reaction['reaction_id']
                field_name = f"case_{case_idx}_reaction_{reaction_idx}"
                
                selected_complexes = [v for v in request.POST.getlist(field_name) if v]
                # Delete old resolutions for this gene+reaction, then create new ones
                GPRAmbiguousCaseResolution.objects.filter(
                    gpr_job=job, genome_gene=gene, reaction_id=reaction_id
                ).delete()
                if selected_complexes:
                    for cpx_id in selected_complexes:
                        GPRAmbiguousCaseResolution.objects.create(
                            gpr_job=job,
                            genome_gene=gene,
                            reaction_id=reaction_id,
                            selected_complex=cpx_id,
                        )
                    selections[(gene, reaction_id)] = selected_complexes
        
        # --- Build GPR rules from user selections ---
        # Load complex stoichiometry
        import pickle
        from django.conf import settings
        complex_to_stoichiometry = {}
        stoich_path = os.path.join(settings.BASE_DIR, 'data', 'gpr_reference', 'complex_stoichiometry.pkl')
        try:
            with open(stoich_path, 'rb') as f:
                complex_to_stoichiometry = pickle.load(f)
        except Exception as e:
            print(f"[resolve POST] Warning: could not load stoichiometry: {e}")

        # Group selections by reaction, building GPR rules
        # For each reaction: collect all genes assigned to the same complex → AND them
        reaction_complex_genes = {}  # {reaction_id: {complex_id: [gene_ids]}}
        for (gene, reaction_id), complex_ids in selections.items():
            for complex_id in complex_ids:
                reaction_complex_genes.setdefault(reaction_id, {}).setdefault(complex_id, []).append(gene)

        # Build reverse map: uniprot → genome gene from BLAST results
        uniprot_to_gene = {}
        for br in job.blast_results.all():
            # br.complex_id actually stores the UniProt ID of the hit
            uid = br.complex_id
            if uid:
                uid_norm = normalize_uniprot_id(uid).upper()
                # Keep best hit per uniprot (lowest evalue)
                if uid_norm not in uniprot_to_gene or br.evalue < uniprot_to_gene[uid_norm][1]:
                    uniprot_to_gene[uid_norm] = (br.gene_id, br.evalue)
        # Simplify to {uniprot: gene}
        uniprot_to_gene = {uid: val[0] for uid, val in uniprot_to_gene.items()}

        # Also include unique cases from the original analysis (already resolved automatically)
        unique_data = []
        if job.gpr_rules_json:
            try:
                existing_rules = json.loads(job.gpr_rules_json)
            except Exception:
                existing_rules = {}
        else:
            existing_rules = {}

        # Normalize ALL incoming rules: convert uppercase operators to lowercase
        def normalize_rule(rule):
            """Normalize rule string to use lowercase and/or operators"""
            if not rule:
                return rule
            import re as _re
            normalized = _re.sub(r'\bAND\b', 'and', rule, flags=_re.IGNORECASE)
            normalized = _re.sub(r'\bOR\b', 'or', normalized, flags=_re.IGNORECASE)
            normalized = _re.sub(r'\s+', ' ', normalized).strip()
            return normalized
        
        # Start with normalized existing rules from unique cases
        gpr_rules = {}
        for rxn_id, rule in existing_rules.items():
            gpr_rules[rxn_id] = normalize_rule(rule)

        for reaction_id, complexes in reaction_complex_genes.items():
            # If multiple complexes contribute to the same reaction => OR them
            or_parts = []
            for complex_id, genes in complexes.items():
                # Start with genes from user selections
                gene_set = set(genes)

                # Find the complex in the original ambiguous data to get
                # the per-member DIAMOND results
                for case in case_details:
                    for cpx in case.get('complexes', []):
                        if cpx.get('complex_id') == complex_id:
                            for br in cpx.get('blast_results', []):
                                g = br.get('genome_gene', '')
                                if g:
                                    # Clean any stoichiometry notation from genes
                                    g_clean = str(g).strip()
                                    if '(' in g_clean and ')' in g_clean:
                                        g_clean = g_clean[:g_clean.index('(')]
                                    g_clean = g_clean.strip()
                                    if g_clean:
                                        gene_set.add(g_clean)

                # Fallback: also cross-reference stoichiometry UniProts
                stoich_str = complex_to_stoichiometry.get(complex_id, '')
                stoich_dict = parse_stoichiometry_string(stoich_str) if stoich_str else {}
                for uid in stoich_dict:
                    # uid is the UniProt ID (parsed without stoichiometry notation)
                    uid_clean = str(uid).strip()
                    # Remove any trailing (n) notation if it somehow got there
                    if '(' in uid_clean:
                        uid_clean = uid_clean[:uid_clean.index('(')]
                    
                    mapped = uniprot_to_gene.get(uid_clean.upper(), '')
                    if mapped and mapped not in gene_set:
                        # Also clean the mapped gene in case it has (n) notation
                        mapped_clean = str(mapped).strip()
                        if '(' in mapped_clean:
                            mapped_clean = mapped_clean[:mapped_clean.index('(')]
                        if mapped_clean:
                            gene_set.add(mapped_clean)

                gene_parts = sorted(gene_set)
                # Filter out empty gene names and remove any stoichiometry notation
                cleaned_parts = []
                for g in gene_parts:
                    g_clean = str(g).strip()
                    # Remove stoichiometry notation like (4) if present
                    if '(' in g_clean and ')' in g_clean:
                        g_clean = g_clean[:g_clean.index('(')]
                    g_clean = g_clean.strip()
                    if g_clean:
                        cleaned_parts.append(g_clean)
                gene_parts = cleaned_parts
                
                if not gene_parts:
                    print(f"[resolve POST] Warning: No valid genes for complex {complex_id} in reaction {reaction_id}")
                    continue
                    
                and_rule = ' and '.join(str(g).strip() for g in gene_parts)
                if len(gene_parts) > 1:
                    and_rule = f"({and_rule})"
                or_parts.append(and_rule)
            rule = ' or '.join(or_parts) if or_parts else ''
            if rule:
                # Normalize the rule before adding it
                normalized_rule = normalize_rule(rule)
                print(f"[resolve POST] Setting {reaction_id} = {normalized_rule}")
                gpr_rules[reaction_id] = normalized_rule

        # --- Update the model file ---
        if gpr_rules:
            try:
                input_model_path = job.metabolic_model.path
                _, input_ext = os.path.splitext(input_model_path)
                if not input_ext:
                    input_ext = '.sbml'

                from django.conf import settings
                output_dir = os.path.join(settings.MEDIA_ROOT, 'reconstructed_models')
                os.makedirs(output_dir, exist_ok=True)
                output_path = os.path.join(output_dir, f'model_gpr_updated_{job.id}{input_ext}')

                result = update_model_with_gpr_rules(input_model_path, gpr_rules, output_path)

                if result['success']:
                    # Use the merged rules (which include preserved OR branches from the original GPR)
                    final_rules = result.get('merged_rules', gpr_rules)
                    # Update gpr_rules so build_complex_summary shows the actual final rules
                    gpr_rules.update(final_rules)
                    job.gpr_rules_json = json.dumps(gpr_rules)
                    job.updated_model_path = output_path
                    job.message = f"GPR rules applied to {result['updates']} reactions. Model ready for download."
                    print(f"[resolve POST] Updated model saved: {output_path}")
                else:
                    job.message = f"Selections saved. Model update failed: {result.get('error', 'unknown')}"
                    print(f"[resolve POST] Model update failed: {result.get('error')}")
            except Exception as e:
                import traceback
                traceback.print_exc()
                job.message = f"Selections saved but model update failed: {str(e)}"
                print(f"[resolve POST] Error: {e}")
        else:
            job.message = "No GPR rules generated from selections."

        # Build complex summary for display
        # Only include cases where a complex was actually selected
        resolved_cases = []
        for case_idx, case in enumerate(case_details):
            for reaction_idx, reaction in enumerate(case.get('reactions', [])):
                rxn_id = reaction['reaction_id']
                sel_cpxs = selections.get((case['genome_gene'], rxn_id), [])
                for sel_cpx in sel_cpxs:
                    resolved_cases.append({
                        'genome_gene': case['genome_gene'],
                        'uniprot_id': case.get('uniprot_id', ''),
                        'possible_complexes': case.get('complexes', []),
                        'reactions_with_gene': [reaction],
                        'selected_complex': sel_cpx,
                    })
        new_summary_rows = build_complex_summary(
            resolved_cases,
            complex_to_stoichiometry,
            gpr_rules,
            complex_descriptions=complex_descriptions,
            blast_results=job.blast_results.all()
        )

        # Merge with existing summary from unique (auto-resolved) cases
        existing_summary = []
        if job.complex_summary_json:
            try:
                existing_summary = json.loads(job.complex_summary_json)
            except Exception:
                pass

        # Combine: keep existing rows for reactions not in the new summary,
        # then append all new rows. Update gpr_rule on existing rows too.
        new_rxn_ids = {r['reaction_id'] for r in new_summary_rows}
        merged = [r for r in existing_summary if r['reaction_id'] not in new_rxn_ids]
        # Update gpr_rule for existing rows (the rules dict may have changed)
        for row in merged:
            rxn_id = row.get('reaction_id', '')
            if rxn_id in gpr_rules:
                row['gpr_rule'] = gpr_rules[rxn_id]
        merged.extend(new_summary_rows)
        merged.sort(key=lambda r: r.get('reaction_id', ''))
        job.complex_summary_json = json.dumps(merged)

        job.status = 'completed'
        try:
            job.has_unresolved_ambiguities = False  # Mark ambiguities as resolved
        except AttributeError:
            pass  # Column not in DB yet (migration pending)
        job.save()
        
        return redirect('gpr:job_detail', job_id=job_id)
    
    context = {
        'job': job,
        'case_details': case_details,
        'reactions_list': reactions_list,  # NEW: reactions grouped with all ambiguous genes
        'complex_metadata': complex_metadata_dict,
        'page_title': f'Resolve Ambiguities - {job.organism_name}'
    }
    return render(request, 'GPR/resolve_ambiguities.html', context)


def gpr_user_selections(request, job_id):
    """Display user's selected complexes for review"""
    job = _get_own_gpr_job_or_404(request, job_id)
    
    resolutions = GPRAmbiguousCaseResolution.objects.filter(gpr_job=job).values(
        'genome_gene', 'reaction_id', 'selected_complex'
    ).order_by('genome_gene', 'reaction_id')
    
    # Get complex metadata
    complex_metadata = {}
    for resolution in resolutions:
        complex_id = resolution['selected_complex']
        if complex_id not in complex_metadata:
            try:
                metadata = ComplexMetadata.objects.get(complex_id=complex_id)
                complex_metadata[complex_id] = {
                    'name': metadata.name,
                    'aliases': metadata.get_aliases_list()
                }
            except ComplexMetadata.DoesNotExist:
                complex_metadata[complex_id] = {
                    'name': complex_id,
                    'aliases': []
                }
    
    context = {
        'job': job,
        'resolutions': resolutions,
        'complex_metadata': complex_metadata,
        'page_title': f'Complex Assignments - {job.organism_name}'
    }
    return render(request, 'GPR/user_selections.html', context)


@require_http_methods(["POST"])
def edit_gpr_rule(request, job_id):
    """AJAX endpoint: validate and apply a manually edited GPR rule for one reaction."""
    job = _get_own_gpr_job_or_404(request, job_id)

    if job.status != 'completed':
        return JsonResponse({'ok': False, 'error': 'Job is not in completed state.'}, status=400)

    try:
        data = json.loads(request.body)
    except json.JSONDecodeError:
        return JsonResponse({'ok': False, 'error': 'Invalid JSON.'}, status=400)

    reaction_id = data.get('reaction_id', '').strip()
    new_rule = data.get('gpr_rule', '').strip()

    if not reaction_id:
        return JsonResponse({'ok': False, 'error': 'Missing reaction_id.'}, status=400)
    if not new_rule:
        return JsonResponse({'ok': False, 'error': 'GPR rule cannot be empty.'}, status=400)

    # --- Validate the rule is well-formed ---
    import re
    # Only allow gene identifiers, parentheses, and/or operators
    # Strip to tokens and verify each one
    tokens = new_rule.replace('(', ' ( ').replace(')', ' ) ').split()
    paren_depth = 0
    prev_token = None
    operators = {'and', 'or'}
    for tok in tokens:
        if tok == '(':
            paren_depth += 1
        elif tok == ')':
            paren_depth -= 1
            if paren_depth < 0:
                return JsonResponse({'ok': False, 'error': 'Unmatched closing parenthesis.'}, status=400)
        elif tok.lower() in operators:
            if prev_token is None or prev_token == '(' or (prev_token and prev_token.lower() in operators):
                return JsonResponse({'ok': False, 'error': f'Unexpected operator "{tok}".'}, status=400)
        else:
            # Gene identifier — must be alphanumeric/dot/underscore/dash
            if not re.match(r'^[\w.\-]+$', tok):
                return JsonResponse({'ok': False, 'error': f'Invalid character in gene name "{tok}".'}, status=400)
        prev_token = tok

    if paren_depth != 0:
        return JsonResponse({'ok': False, 'error': 'Unmatched opening parenthesis.'}, status=400)

    # Also try COBRApy's parser as a second validation layer
    try:
        import cobra
        cobra.core.gene.GPR.from_string(new_rule)
    except Exception as e:
        return JsonResponse({'ok': False, 'error': f'COBRApy rejected the rule: {str(e)}'}, status=400)

    # --- Apply the new rule ---
    # 1. Update the stored gpr_rules_json
    gpr_rules = {}
    if job.gpr_rules_json:
        try:
            gpr_rules = json.loads(job.gpr_rules_json)
        except Exception:
            pass
    gpr_rules[reaction_id] = new_rule

    # 2. Update the complex_summary_json so the UI stays in sync
    complex_summary = []
    if job.complex_summary_json:
        try:
            complex_summary = json.loads(job.complex_summary_json)
        except Exception:
            pass
    for row in complex_summary:
        if row.get('reaction_id') == reaction_id:
            row['gpr_rule'] = new_rule
    job.complex_summary_json = json.dumps(complex_summary)
    job.gpr_rules_json = json.dumps(gpr_rules)

    # 3. Re-write the model file with ALL current rules
    try:
        input_model_path = job.metabolic_model.path
        output_path = job.updated_model_path
        if not output_path:
            _, input_ext = os.path.splitext(input_model_path)
            if not input_ext:
                input_ext = '.sbml'
            from django.conf import settings as dj_settings
            output_dir = os.path.join(dj_settings.MEDIA_ROOT, 'reconstructed_models')
            os.makedirs(output_dir, exist_ok=True)
            output_path = os.path.join(output_dir, f'model_gpr_updated_{job.id}{input_ext}')

        # For a single manual edit we write the rule directly (no merge)
        # because the user is providing the final desired rule.
        result = update_model_with_gpr_rules(input_model_path, gpr_rules, output_path, merge=False)

        if result['success']:
            job.updated_model_path = output_path
            job.save()
            return JsonResponse({'ok': True, 'gpr_rule': new_rule, 'updates': result['updates']})
        else:
            return JsonResponse({'ok': False, 'error': f"Model update failed: {result.get('error', 'unknown')}"}, status=500)
    except Exception as e:
        import traceback
        traceback.print_exc()
        return JsonResponse({'ok': False, 'error': str(e)}, status=500)
