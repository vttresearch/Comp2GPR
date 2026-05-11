"""
GPR Analysis Services - BLAST operations and data processing
"""
import os
import re
import subprocess
import tempfile
import json
import csv
from pathlib import Path
from Bio import SeqIO
from django.conf import settings

def extract_genes_from_model(model_path):
    """
    Extract gene IDs from metabolic model (SBML or JSON format)
    Returns list of gene IDs and their sequences if available
    """
    genes = []
    
    try:
        # Try JSON format first
        if model_path.endswith('.json'):
            with open(model_path, 'r') as f:
                data = json.load(f)
                if 'genes' in data:
                    for gene in data['genes']:
                        genes.append(gene.get('id', gene))
                elif 'reactions' in data:
                    # Extract genes from reaction GPR
                    for reaction in data['reactions']:
                        gpr = reaction.get('gpr', '')
                        # Simple parsing - extract identifiers
                        gene_ids = parse_gpr_rule(gpr)
                        genes.extend(gene_ids)
        
        # Try SBML format
        elif model_path.endswith('.xml') or model_path.endswith('.sbml'):
            try:
                import cobra
                model = cobra.io.read_sbml_model(model_path)
                genes = [gene.id for gene in model.genes]
            except:
                # Fallback: extract from XML manually
                with open(model_path, 'r') as f:
                    content = f.read()
                    # Simple regex extraction
                    import re
                    gene_ids = re.findall(r'<fbc:geneProduct fbc:id="([^"]+)"', content)
                    genes.extend(gene_ids)
    except Exception as e:
        raise Exception(f"Error extracting genes from model: {str(e)}")
    
    return list(set(genes))  # Remove duplicates


def _extract_fasta_identifiers(record):
    """
    Extract all possible identifiers from a FASTA record header.
    Handles UniProt format: >tr|G0RHX3|G0RHX3_HYPJQ description GN=TRIREDRAFT_60676 ...
    Also handles plain headers: >gene_name

    Returns a set of candidate IDs (all uppercased for matching).
    """
    import re
    candidates = set()
    raw_id = record.id  # e.g. "tr|G0RHX3|G0RHX3_HYPJQ"

    # 1. Raw record.id
    candidates.add(raw_id)

    # 2. Pipe-delimited UniProt format: sp|ACCESSION|ENTRY_ORGANISM
    parts = raw_id.split('|')
    if len(parts) >= 2:
        candidates.add(parts[1])  # UniProt accession, e.g. "G0RHX3"
    if len(parts) >= 3:
        candidates.add(parts[2])  # Entry name, e.g. "G0RHX3_HYPJQ"
        # Also the part before the underscore
        entry_base = parts[2].split('_')[0]
        candidates.add(entry_base)  # e.g. "G0RHX3"

    # 3. GN= gene name from description
    full_header = record.description  # full header line
    gn_match = re.search(r'\bGN=(\S+)', full_header)
    if gn_match:
        candidates.add(gn_match.group(1))  # e.g. "TRIREDRAFT_60676"

    # 4. Version-suffix stripped variants (e.g. "xxxx.1" -> "xxxx")
    for c in list(candidates):
        stripped = re.sub(r'\.\d+$', '', c)
        if stripped != c:
            candidates.add(stripped)

    # Remove empty strings
    candidates.discard('')
    return {c.upper() for c in candidates}


def _extract_model_gene_identifiers(gene_id):
    """
    Extract candidate identifiers from a model gene ID.
    SBML models often prefix gene IDs (e.g. G_G0RHX3 for G0RHX3, 
    because SBML identifiers cannot start with a digit).

    Returns a set of candidate IDs (all uppercased for matching).
    """
    candidates = {gene_id}

    # Strip common SBML prefixes: G_, g_, gene_, GENE_
    import re
    stripped = re.sub(r'^[Gg]_', '', gene_id)
    if stripped != gene_id:
        candidates.add(stripped)

    stripped2 = re.sub(r'^(?:gene|GENE)_', '', gene_id, flags=re.IGNORECASE)
    if stripped2 != gene_id:
        candidates.add(stripped2)

    # Strip version suffixes (e.g. "xxxx.1" -> "xxxx")
    for c in list(candidates):
        stripped_ver = re.sub(r'\.\d+$', '', c)
        if stripped_ver != c:
            candidates.add(stripped_ver)

    candidates.discard('')
    return {c.upper() for c in candidates}


def extract_model_genes_from_fasta(genome_fasta_path, model_genes):
    """
    Extract ONLY the genes that are in the metabolic model from the genome FASTA file.
    Returns a temporary FASTA file containing only those sequences.

    Performs flexible ID matching:
    - UniProt pipe-delimited headers (tr|ACC|ENTRY_ORG)
    - SBML-prefixed model gene IDs (G_ACC, gene_ACC)
    - GN= gene names from FASTA descriptions
    
    Args:
        genome_fasta_path: Path to complete genome FASTA file
        model_genes: List/set of gene IDs to extract (from metabolic model)
    
    Returns:
        (temp_fasta_path, extracted_count, gene_id_map) -
            path to filtered FASTA, count of genes found,
            and dict mapping the FASTA label used → original model gene ID
    """
    # Build a lookup: {candidate_upper: original_model_gene_id}
    model_lookup = {}
    for gene_id in model_genes:
        for candidate in _extract_model_gene_identifiers(gene_id):
            model_lookup[candidate] = gene_id

    temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
    extracted_count = 0
    gene_id_map = {}  # {fasta_label: model_gene_id}
    
    try:
        for record in SeqIO.parse(genome_fasta_path, "fasta"):
            fasta_ids = _extract_fasta_identifiers(record)

            # Find the first matching model gene
            matched_model_gene = None
            for fid in fasta_ids:
                if fid in model_lookup:
                    matched_model_gene = model_lookup[fid]
                    break

            if matched_model_gene:
                # Write with the MODEL gene ID as the header so downstream
                # BLAST results link back to the model
                temp_fasta.write(f">{matched_model_gene}\n{record.seq}\n")
                gene_id_map[matched_model_gene] = matched_model_gene
                extracted_count += 1

        temp_fasta.close()
    except Exception as e:
        temp_fasta.close()
        os.unlink(temp_fasta.name)
        raise Exception(f"Error extracting model genes from genome FASTA: {str(e)}")
    
    if extracted_count == 0:
        os.unlink(temp_fasta.name)
        # Provide a diagnostic message
        sample_fasta_ids = []
        try:
            for i, record in enumerate(SeqIO.parse(genome_fasta_path, "fasta")):
                sample_fasta_ids.append(record.id)
                if i >= 4:
                    break
        except Exception:
            pass
        sample_model = list(model_genes)[:5]
        raise Exception(
            f"No metabolic model genes found in genome FASTA file. "
            f"Sample FASTA IDs: {sample_fasta_ids}. "
            f"Sample model genes: {sample_model}. "
            f"Check that gene naming conventions match."
        )
    
    print(f"[extract_model_genes] Matched {extracted_count}/{len(model_genes)} model genes in genome FASTA")
    return temp_fasta.name, extracted_count


def _build_named_reaction_string(reaction):
    """
    Build a reaction equation string using real metabolite names instead of IDs.
    E.g. "2.0 ATP + Water --> ADP + Phosphate" instead of "2.0 atp_c + h2o_c --> adp_c + pi_c"
    """
    reactants = []
    products = []
    for met, coeff in reaction.metabolites.items():
        # Use name if available, otherwise fall back to id
        label = met.name if met.name else met.id
        abs_coeff = abs(coeff)
        if abs_coeff != 1.0:
            entry = f"{abs_coeff:g} {label}"
        else:
            entry = label
        if coeff < 0:
            reactants.append(entry)
        else:
            products.append(entry)
    
    arrow = '<=>' if reaction.lower_bound < 0 else '-->'
    return f"{' + '.join(reactants)} {arrow} {' + '.join(products)}"


def extract_reactions_from_model(model_path):
    """
    Extract reactions with their GPR rules from metabolic model
    Returns dict: {reaction_id: {id, name, gpr, formula, genes: [list]}}
    """
    reactions = {}
    print(f"[extract_reactions] Reading model from: {model_path}")
    
    try:
        if model_path.endswith('.json'):
            print(f"[extract_reactions] Detected JSON format")
            with open(model_path, 'r') as f:
                data = json.load(f)
                print(f"[extract_reactions] JSON loaded, keys: {list(data.keys())}")
                
                if 'reactions' in data:
                    print(f"[extract_reactions] Found {len(data['reactions'])} reactions in JSON")
                    for reaction in data['reactions']:
                        rx_id = reaction.get('id', '')
                        gpr = reaction.get('gpr', '')
                        genes_in_reaction = parse_gpr_rule(gpr) if gpr else []
                        
                        reactions[rx_id] = {
                            'id': rx_id,
                            'name': reaction.get('name', ''),
                            'gpr': gpr,
                            'genes': genes_in_reaction,
                            'formula': reaction.get('formula', ''),
                            'metabolites': reaction.get('metabolites', {})
                        }
                else:
                    print(f"[extract_reactions] No 'reactions' key found in JSON")
        
        elif model_path.endswith('.xml') or model_path.endswith('.sbml'):
            print(f"[extract_reactions] Detected SBML format")
            try:
                import cobra
                model = cobra.io.read_sbml_model(model_path)
                print(f"[extract_reactions] SBML loaded, found {len(model.reactions)} reactions")
                
                for reaction in model.reactions:
                    gpr = str(reaction.gpr)
                    genes_in_reaction = parse_gpr_rule(gpr)
                    
                    # Build formula with real compound names
                    formula = _build_named_reaction_string(reaction)
                    
                    reactions[reaction.id] = {
                        'id': reaction.id,
                        'name': reaction.name,
                        'gpr': gpr,
                        'genes': genes_in_reaction,
                        'formula': formula,
                        'metabolites': {m.id: c for m, c in reaction.metabolites.items()}
                    }
            except Exception as e:
                print(f"[extract_reactions] COBRApy failed: {e}, trying XML fallback")
                # Fallback: basic XML extraction
                pass
        else:
            print(f"[extract_reactions] Unknown format: {model_path}")
            
    except Exception as e:
        print(f"[extract_reactions] Error extracting reactions: {str(e)}")
        import traceback
        traceback.print_exc()
    
    print(f"[extract_reactions] Extracted {len(reactions)} total reactions")
    return reactions


def parse_gpr_rule(gpr_string):
    """
    Parse GPR rule string to extract gene IDs
    Handles boolean operators (and, or, parentheses)
    """
    import re
    # Remove boolean operators and parentheses, extract identifiers
    identifiers = re.findall(r'\b[\w\.\-]+\b', gpr_string)
    # Filter out 'and' and 'or' keywords
    genes = [g for g in identifiers if g.lower() not in ['and', 'or']]
    return genes


def create_blast_db(fasta_file, db_name):
    """
    Create a BLAST database from a FASTA file
    Returns the path to the database
    """
    db_dir = tempfile.mkdtemp(prefix='blast_db_')
    db_path = os.path.join(db_dir, db_name)
    
    try:
        # Run makeblastdb
        cmd = [
            'makeblastdb',
            '-in', str(fasta_file),
            '-dbtype', 'prot',
            '-out', db_path,
            '-title', db_name
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return db_path, db_dir
    except subprocess.CalledProcessError as e:
        raise Exception(f"Failed to create BLAST database: {e.stderr}")


def run_blast_search(query_fasta, db_path, evalue_threshold=1e-5):
    """
    Run BLAST search to match genes against complexity database
    Returns list of results, sorted by evalue (best first)
    """
    results = []
    
    # Create output file
    output_file = tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False)
    output_path = output_file.name
    output_file.close()
    
    if db_path is None:
        db_path = str(settings.BLAST_COMPLEXES_DB)
    
    try:
        # Run blast
        cmd = [
            'blastp',
            '-query', str(query_fasta),
            '-db', db_path,
            '-evalue', str(evalue_threshold),
            '-outfmt', '6 qseqid sseqid evalue bitscore pident qcovs',
            '-out', output_path,
            '-num_threads', '4'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Parse results
        if os.path.getsize(output_path) > 0:
            with open(output_path, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t', 
                                      fieldnames=['query', 'subject', 'evalue', 'bitscore', 'pident', 'qcovs'])
                for row in reader:
                    # Extract just the UniProt accession ID from full descriptor
                    uniprot_id = extract_uniprot_accession(row['subject'])
                    results.append({
                        'gene_id': row['query'],
                        'uniprot_id': uniprot_id,
                        'evalue': float(row['evalue']),
                        'bitscore': float(row['bitscore']),
                        'identity': float(row['pident']),
                        'query_coverage': float(row['qcovs'])
                    })
        
        # Sort by evalue (ascending - smallest first = best matches)
        results.sort(key=lambda x: x['evalue'])
        return results
    
    except subprocess.CalledProcessError as e:
        raise Exception(f"BLAST search failed: {e.stderr}")
    finally:
        if os.path.exists(output_path):
            os.unlink(output_path)


def get_best_blast_hit(blast_results, evalue_threshold=1e-5):
    """
    Filter BLAST results to get the best hit for each query gene, ensuring
    each complex protein (UniProt ID) is assigned to at most ONE gene.

    When multiple genes share the same best UniProt hit, the gene with the
    lowest e-value keeps it; the others fall back to their next-best
    *unassigned* UniProt hit.  This guarantees a 1-to-1 mapping.

    Returns dict: {gene_id: best_hit_dict}
    """
    # Step 1: collect ALL valid hits per gene, sorted by e-value (ascending)
    hits_per_gene = {}  # gene_id -> [hit, hit, ...]
    for result in blast_results:
        if result['evalue'] <= evalue_threshold:
            gene_id = result['gene_id']
            hits_per_gene.setdefault(gene_id, []).append(result)

    # Each gene's list is already sorted (blast_results is sorted by evalue),
    # but let's be safe
    for gene_id in hits_per_gene:
        hits_per_gene[gene_id].sort(key=lambda r: r['evalue'])

    # Step 2: greedy 1-to-1 assignment
    # Priority: genes whose top hit has the lowest e-value go first
    genes_sorted = sorted(
        hits_per_gene.keys(),
        key=lambda g: hits_per_gene[g][0]['evalue']
    )

    assigned_uniprots = set()   # UniProt IDs already claimed
    best_hits = {}              # gene_id -> hit

    for gene_id in genes_sorted:
        for hit in hits_per_gene[gene_id]:
            uid = hit['uniprot_id']
            if uid not in assigned_uniprots:
                best_hits[gene_id] = hit
                assigned_uniprots.add(uid)
                break
        # If every hit for this gene is already claimed, skip it

    # Debug: report reassignments
    reassigned = 0
    for gene_id in best_hits:
        original_top = hits_per_gene[gene_id][0]['uniprot_id']
        if best_hits[gene_id]['uniprot_id'] != original_top:
            reassigned += 1
            print(f"[get_best_blast_hit] {gene_id}: reassigned from "
                  f"{original_top} (taken) -> {best_hits[gene_id]['uniprot_id']} "
                  f"(e={best_hits[gene_id]['evalue']:.1e})")

    skipped = len(hits_per_gene) - len(best_hits)
    print(f"[get_best_blast_hit] {len(best_hits)} genes assigned, "
          f"{reassigned} reassigned to next-best hit, "
          f"{skipped} genes skipped (all hits taken)")

    return best_hits


# ============================================================================
# DIAMOND functions (faster alternative to BLAST) - 20-100x speedup
# ============================================================================

def create_diamond_db(fasta_file, db_name):
    """
    Create a DIAMOND database from a FASTA file (much faster than BLAST)
    Returns the path to the database (.dmnd file)
    """
    db_dir = tempfile.mkdtemp(prefix='diamond_db_')
    db_path = os.path.join(db_dir, db_name)
    
    try:
        # Run diamond makedb
        cmd = [
            'diamond', 'makedb',
            '--in', str(fasta_file),
            '--db', db_path,
            '--threads', '4'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        # Diamond adds .dmnd extension automatically
        dmnd_path = db_path + '.dmnd'
        return dmnd_path, db_dir
    except subprocess.CalledProcessError as e:
        raise Exception(f"Failed to create DIAMOND database: {e.stderr}")
    except FileNotFoundError:
        raise Exception("DIAMOND not installed. Install with: sudo apt-get install diamond-aligner")


def run_diamond_search(query_fasta, db_path, evalue_threshold=1e-5):
    """
    Run DIAMOND search (20-100x faster than BLAST, similar sensitivity)
    Returns list of results in same format as run_blast_search
    """
    results = []
    
    # Create output file
    output_file = tempfile.NamedTemporaryFile(mode='w', suffix='.m8', delete=False)
    output_path = output_file.name
    output_file.close()
    
    try:
        # Run diamond blastp
        cmd = [
            'diamond', 'blastp',
            '--db', str(db_path),
            '--query', str(query_fasta),
            '--out', output_path,
            '--evalue', str(evalue_threshold),
            '--outfmt', '6', 'qseqid', 'sseqid', 'evalue', 'bitscore', 'pident', 'qcovhsp',
            '--threads', '4',
            '--sensitive'  # Balance speed and sensitivity
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Parse results (same format as BLAST)
        if os.path.getsize(output_path) > 0:
            with open(output_path, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t', 
                                      fieldnames=['query', 'subject', 'evalue', 'bitscore', 'pident', 'qcovs'])
                for row in reader:
                    # Extract just the UniProt accession ID from full descriptor
                    uniprot_id = extract_uniprot_accession(row['subject'])
                    results.append({
                        'gene_id': row['query'],
                        'uniprot_id': uniprot_id,
                        'evalue': float(row['evalue']),
                        'bitscore': float(row['bitscore']),
                        'identity': float(row['pident']),
                        'query_coverage': float(row['qcovs'])
                    })
        
        # Sort by evalue (ascending - smallest first = best matches)
        results.sort(key=lambda x: x['evalue'])
        return results
    
    except subprocess.CalledProcessError as e:
        raise Exception(f"DIAMOND search failed: {e.stderr}")
    except FileNotFoundError:
        raise Exception("DIAMOND not installed. Install with: sudo apt-get install diamond-aligner")
    finally:
        if os.path.exists(output_path):
            os.unlink(output_path)


def extract_uniprot_accession(uniprot_descriptor):
    """
    Extract UniProt accession ID from various formats:
    - 'P00904' -> 'P00904'
    - 'sp|P00904|TRPGD_ECOLI' -> 'P00904'
    - 'tr|Q8W123|ABC1_HUMAN' -> 'Q8W123'
    """
    if not uniprot_descriptor:
        return ''
    
    # If it contains pipes, it's a full descriptor - extract the accession (2nd element)
    if '|' in uniprot_descriptor:
        parts = uniprot_descriptor.split('|')
        if len(parts) >= 2:
            return parts[1]  # Return the accession ID
    
    # Otherwise return as-is (it's already just the accession)
    return uniprot_descriptor


def normalize_uniprot_id(uid):
    """Strip chain/isoform suffixes like -PRO_0000006048 from UniProt accessions."""
    if not uid:
        return uid
    if '-' in uid:
        return uid.split('-')[0]
    return uid


def parse_stoichiometry_string(stoich_string):
    """
    Parse stoichiometry string like 'P15790(2)|P25368(0)|P38930(1)|P39520(0)|P43639(1)|P53254(0)'
    Returns dict: {uniprot_id: stoichiometry}
    """
    stoichiometry_dict = {}
    
    if not stoich_string:
        return stoichiometry_dict
    
    entries = stoich_string.split('|')
    for entry in entries:
        if '(' in entry and ')' in entry:
            raw_id = entry[:entry.index('(')]
            uniprot_id = normalize_uniprot_id(raw_id)
            stoich_value = entry[entry.index('(')+1:entry.index(')')]
            try:
                stoichiometry_dict[uniprot_id] = int(stoich_value)
            except ValueError:
                pass
    
    return stoichiometry_dict


def find_all_complexes_for_uniprot(uniprot_id, complex_to_stoichiometry):
    """
    Find ALL complexes that contain a specific UniProt ID
    
    Args:
        uniprot_id: UniProt protein ID (e.g., 'P15790')
        complex_to_stoichiometry: Dict mapping {complex_id: 'P15790(2)|P25368(0)|...'}
    
    Returns:
        List of complex IDs that contain this UniProt protein
    """
    matching_complexes = []
    
    for complex_id, stoich_string in complex_to_stoichiometry.items():
        # Parse stoichiometry string to extract all UniProt IDs in this complex
        if stoich_string:
            entries = stoich_string.split('|')
            for entry in entries:
                if '(' in entry:
                    raw_id = entry[:entry.index('(')]
                    protein_id = normalize_uniprot_id(raw_id)
                    # Case-insensitive comparison
                    if protein_id.upper() == uniprot_id.upper():
                        matching_complexes.append(complex_id)
                        break  # Found in this complex, move to next complex
    
    return matching_complexes


def extract_sequences_from_complexes_db(uniprot_ids, complexes_fasta_path):
    """
    Extract protein sequences from the complexes FASTA database.
    FASTA headers use the format: sp|Q9W1F4|THOC5_DROME
    Returns: {uniprot_id: sequence}
    """
    sequences = {}
    # Normalize input IDs for matching
    target_ids = {normalize_uniprot_id(uid).upper(): uid for uid in uniprot_ids}
    
    try:
        for record in SeqIO.parse(complexes_fasta_path, "fasta"):
            # Extract accession from sp|P04037|COX5A_YEAST  → P04037
            parts = record.id.split('|')
            if len(parts) >= 2:
                accession = parts[1].upper()
            else:
                accession = record.id.upper()
            accession = normalize_uniprot_id(accession).upper()

            if accession in target_ids:
                original_uid = target_ids[accession]
                sequences[original_uid] = str(record.seq)
    
    except Exception as e:
        raise Exception(f"Error extracting sequences from complexes DB: {str(e)}")
    
    return sequences


def blast_sequences_against_genome(sequences_dict, genome_blast_db, evalue_threshold=1e-5, existing_genes=None):
    """
    Search a dictionary of sequences against the genome database using DIAMOND (faster than BLAST)
    Returns: {sequence_id: [list of hit dicts with genome gene, evalue, bitscore, identity, query_coverage]}
    
    If existing_genes is provided, genome gene IDs will be normalized to match
    model gene IDs (e.g. 'xxxx.1' -> 'xxxx' if 'xxxx' is in existing_genes).
    """
    # Build a lookup for normalizing genome gene IDs to model gene IDs
    _gene_normalize_map = {}
    if existing_genes:
        for mg in existing_genes:
            _gene_normalize_map[mg.upper()] = mg
            # Also index without version suffix
            stripped = re.sub(r'\.\d+$', '', mg)
            if stripped != mg:
                _gene_normalize_map[stripped.upper()] = mg

    def _normalize_genome_gene(raw_id):
        """Map a raw genome hit ID back to the model gene ID if possible,
        or strip version suffix (e.g. .1) to get a clean gene name."""
        upper = raw_id.upper()
        if upper in _gene_normalize_map:
            return _gene_normalize_map[upper]
        stripped = re.sub(r'\.\d+$', '', raw_id)
        if stripped.upper() in _gene_normalize_map:
            return _gene_normalize_map[stripped.upper()]
        # Even if the gene isn't in the model, strip version suffix
        # so that e.g. ETS00265.1 becomes ETS00265
        if stripped != raw_id:
            return stripped
        return raw_id

    # Create temporary FASTA file with sequences
    temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
    try:
        for seq_id, sequence in sequences_dict.items():
            temp_fasta.write(f">{seq_id}\n{sequence}\n")
        temp_fasta.close()
        
        # Run DIAMOND search against genome database (20-100x faster)
        blast_results = run_diamond_search(temp_fasta.name, genome_blast_db, evalue_threshold)
        
        # Group results by query (complex protein)
        matches_by_protein = {}
        for result in blast_results:
            complex_protein = result['gene_id']
            genome_gene = _normalize_genome_gene(result['uniprot_id'])
            
            if complex_protein not in matches_by_protein:
                matches_by_protein[complex_protein] = []
            
            # Store full hit information including evalue
            hit_info = {
                'genome_gene': genome_gene,
                'evalue': result['evalue'],
                'bitscore': result['bitscore'],
                'identity': result['identity'],
                'query_coverage': result['query_coverage']
            }
            matches_by_protein[complex_protein].append(hit_info)
        
        return matches_by_protein
    
    finally:
        if os.path.exists(temp_fasta.name):
            os.unlink(temp_fasta.name)


def deduplicate_ambiguous_cases(ambiguous_data):
    """
    After BLAST results have been populated, check each ambiguous case:
    if all complexes for a gene would produce the same GPR rule (same set
    of matched genome genes), auto-resolve by keeping the first complex
    and moving the case to unique_cases.

    Modifies ambiguous_data in place and returns it.
    """
    still_ambiguous = []
    promoted_to_unique = 0

    for case in ambiguous_data['ambiguous_cases']:
        complexes = case.get('possible_complexes', [])
        if len(complexes) <= 1:
            # Already unique
            ambiguous_data['unique_cases'].append(case)
            promoted_to_unique += 1
            continue

        # Build the "gene signature" each complex produces
        signatures = []
        for cpx in complexes:
            matched_genes = frozenset(
                br['genome_gene']
                for br in cpx.get('blast_results', [])
                if br.get('genome_gene')
            )
            signatures.append(matched_genes)

        # Keep only complexes with distinct gene signatures
        seen = {}  # signature → first complex index
        distinct_indices = []
        for idx, sig in enumerate(signatures):
            if sig not in seen:
                seen[sig] = idx
                distinct_indices.append(idx)

        if len(distinct_indices) <= 1:
            # All complexes produce the same GPR rule → auto-resolve with first
            case['possible_complexes'] = [complexes[0]]
            ambiguous_data['unique_cases'].append(case)
            promoted_to_unique += 1
            print(f"[dedup] Gene {case['genome_gene']}: all {len(complexes)} complexes "
                  f"produce same GPR → auto-resolved with {complexes[0]['complex_id']}")
        else:
            # Keep only the distinct complexes for user review
            case['possible_complexes'] = [complexes[i] for i in distinct_indices]
            still_ambiguous.append(case)
            print(f"[dedup] Gene {case['genome_gene']}: {len(complexes)} complexes → "
                  f"{len(distinct_indices)} distinct GPR rules remain ambiguous")

    ambiguous_data['ambiguous_cases'] = still_ambiguous
    print(f"[dedup] SUMMARY: promoted {promoted_to_unique} cases to unique, "
          f"{len(still_ambiguous)} remain ambiguous")
    return ambiguous_data


def collect_ambiguous_cases(best_hits, complex_to_stoichiometry, model_reactions=None):
    """
    Collect cases where a gene matches multiple complexes
    Returns structured data for user review and selection
    
    Args:
        best_hits: Dict from get_best_blast_hit() {uniprot_id: best_hit_dict}
        complex_to_stoichiometry: Dict {complex_id: stoichiometry_string}
        model_reactions: Optional dict of reaction_id → reaction info with GPR rules
    
    Returns:
        {
            'ambiguous_cases': [
                {
                    'genome_gene': 'gene1',
                    'uniprot_id': 'P12345',
                    'blast_hit': {...},
                    'reactions_with_gene': [
                        {
                            'reaction_id': 'R_REACTION1',
                            'gpr': 'gene1 OR gene2',
                            'name': 'Reaction name',
                            'formula': 'A + B -> C'
                        }
                    ],
                    'possible_complexes': [
                        {
                            'complex_id': 'CPX-544',
                            'name': 'Complex name',
                            'aliases': ['alias1', 'alias2'],
                            'stoichiometry': {uniprot_id: count},
                            'blast_results': [
                                {uniprot_id, genome_gene, stoichiometry, evalue, ...}
                            ]
                        }
                    ]
                }
            ],
            'unique_cases': [similar but single complex]
        }
    """
    ambiguous_cases = []
    unique_cases = []
    
    print(f"[collect_ambiguous_cases] Starting: {len(best_hits)} best hits")
    hits_processed = 0
    hits_skipped = 0
    cases_created = 0
    
    for gene_id, best_hit in best_hits.items():
        # best_hits dict is keyed by gene_id
        # best_hit dict contains {gene_id, uniprot_id, evalue, bitscore, identity, query_coverage}
        genome_gene = gene_id  # Use the key directly - it's the genome gene ID
        hit_uniprot_id = best_hit.get('uniprot_id')  # The UniProt complex protein it matched
        
        print(f"[collect_ambiguous_cases] Processing gene: {genome_gene}, matched to UniProt: {hit_uniprot_id}")
        
        if not hit_uniprot_id:
            print(f"[collect_ambiguous_cases]   -> SKIPPED: no uniprot_id")
            hits_skipped += 1
            continue
        
        # Find all complexes containing this UniProt
        complexes_containing = find_all_complexes_for_uniprot(
            hit_uniprot_id, 
            complex_to_stoichiometry
        )
        
        print(f"[collect_ambiguous_cases]   -> Lookup: UniProt ID '{hit_uniprot_id}' in {len(complex_to_stoichiometry)} complexes")
        print(f"[collect_ambiguous_cases]   -> Found {len(complexes_containing)} complexes for this UniProt")
        
        if not complexes_containing:
            print(f"[collect_ambiguous_cases]   -> SKIPPED: no complexes found")
            hits_skipped += 1
            continue
        
        hits_processed += 1
        
        # Get reactions where this gene is annotated (if model provided)
        reactions_with_gene = []
        if model_reactions:
            print(f"[collect_ambiguous_cases]   -> Searching {len(model_reactions)} reactions for gene matches...")
            match_count = 0
            for reaction_id, reaction_info in model_reactions.items():
                gpr = reaction_info.get('gpr', '')
                genes_in_rxn = reaction_info.get('genes', [])
                
                # Check if gene appears in GPR string
                found_in_gpr = False
                if genome_gene in gpr or hit_uniprot_id in gpr:
                    found_in_gpr = True
                else:
                    # Also check if gene is in parsed genes list
                    if genome_gene in genes_in_rxn:
                        found_in_gpr = True
                
                if found_in_gpr:
                    reactions_with_gene.append({
                        'reaction_id': reaction_id,
                        'gpr': gpr,
                        'name': reaction_info.get('name', reaction_id),
                        'formula': reaction_info.get('formula', '')
                    })
                    match_count += 1
            
            print(f"[collect_ambiguous_cases]   -> Found {match_count} reactions involving this gene")
        else:
            print(f"[collect_ambiguous_cases]   -> Skipping reaction lookup: model_reactions not provided")
        
        # Build complex information
        complex_info_list = []
        for complex_id in complexes_containing:
            stoichiometry_str = complex_to_stoichiometry[complex_id]
            stoichiometry = parse_stoichiometry_string(stoichiometry_str)
            
            complex_info_list.append({
                'complex_id': complex_id,
                'name': complex_id,  # Will be replaced by actual name from database
                'aliases': [],       # Will be filled from ComplexMetadata
                'stoichiometry': stoichiometry,
                'blast_results': [],  # Will be filled with BLAST info
                'reactions_with_gene': reactions_with_gene
            })
        
        # Sort cases by number of complexes
        case = {
            'genome_gene': genome_gene,
            'uniprot_id': hit_uniprot_id,
            'blast_hit': best_hit,
            'possible_complexes': complex_info_list,
            'reactions_with_gene': reactions_with_gene
        }
        
        if len(complex_info_list) > 1:
            ambiguous_cases.append(case)
            print(f"[collect_ambiguous_cases]   -> Added to AMBIGUOUS cases ({len(complex_info_list)} complexes)")
        else:
            unique_cases.append(case)
            print(f"[collect_ambiguous_cases]   -> Added to UNIQUE cases (1 complex)")
        
        cases_created += 1
    
    print(f"[collect_ambiguous_cases] SUMMARY: {hits_processed} hits processed, {hits_skipped} skipped")
    print(f"[collect_ambiguous_cases] SUMMARY: {len(ambiguous_cases)} ambiguous cases, {len(unique_cases)} unique cases")
    
    return {
        'ambiguous_cases': ambiguous_cases,
        'unique_cases': unique_cases
    }


def populate_blast_results_for_complexes(ambiguous_data, genome_blast_db, complexes_fasta_path, evalue_threshold=1e-5, existing_genes=None):
    """
    For each ambiguous case, BLAST complex members and populate BLAST results
    
    Args:
        ambiguous_data: Output from collect_ambiguous_cases()
        genome_blast_db: Path to genome BLAST database
        complexes_fasta_path: Path to complexes FASTA
        evalue_threshold: E-value threshold
        existing_genes: Set of model gene IDs for normalizing genome gene names
    
    Returns: Updated ambiguous_data with populated BLAST results
    """
    for case in ambiguous_data['ambiguous_cases']:
        # For each possible complex
        for complex_info in case['possible_complexes']:
            # Extract sequences for complex members
            sequences = extract_sequences_from_complexes_db(
                list(complex_info['stoichiometry'].keys()),
                complexes_fasta_path
            )
            
            if not sequences:
                continue
            
            # BLAST complex members against genome
            genome_matches = blast_sequences_against_genome(
                sequences,
                genome_blast_db,
                evalue_threshold,
                existing_genes=existing_genes
            )
            
            # Populate BLAST results - now includes alternative hits
            blast_results = []
            for uniprot_id, stoich_value in complex_info['stoichiometry'].items():
                matches = genome_matches.get(uniprot_id, [])
                if matches:
                    # Best match (first/lowest evalue)
                    best_match = matches[0]
                    
                    # Collect alternative hits (if any)
                    alternative_hits = []
                    if len(matches) > 1:
                        for alt_match in matches[1:]:
                            alt_str = (
                                f"{alt_match['genome_gene']} "
                                f"(e-val: {alt_match['evalue']:.2e}, "
                                f"identity: {alt_match['identity']:.1f}%, "
                                f"coverage: {alt_match['query_coverage']:.1f}%)"
                            )
                            alternative_hits.append(alt_str)
                    
                    blast_results.append({
                        'uniprot_id': uniprot_id,
                        'genome_gene': best_match['genome_gene'],
                        'stoichiometry': stoich_value,
                        'evalue': best_match['evalue'],
                        'bitscore': best_match['bitscore'],
                        'identity': best_match['identity'],
                        'query_coverage': best_match['query_coverage'],
                        'is_best': True,  # Mark as the selected hit
                        'alternative_hits': alternative_hits  # Show alternatives
                    })
            
            complex_info['blast_results'] = blast_results
    
    return ambiguous_data


def build_gpr_rule_from_complex(best_hit, complex_to_stoichiometry, 
                                genomic_blast_db, existing_genes, complexes_fasta_path, 
                                complex_id, evalue_threshold=1e-5):
    """
    Given a BLAST hit to a uniprot protein and a specific complex ID, 
    build a GPR rule with all complex members that match in the genome.
    
    All debug output is written to media/debug_build_gpr_rule.log
    
    Returns: {
        'gpr_rule': 'gene1 AND gene2(2) AND gene3(1)',
        'new_genes': ['gene2', 'gene3'],
        'complex_id': 'CPX-544',
        'query_gene': 'original_gene',
        'stoichiometry': {uniprot_id: count}
    }
    """
    from datetime import datetime
    
    # Debug log file path
    debug_log_path = os.path.join(settings.MEDIA_ROOT, 'debug_build_gpr_rule.log')
    
    def write_debug(msg):
        """Write debug message to file"""
        try:
            with open(debug_log_path, 'a') as f:
                f.write(f"[{datetime.now().isoformat()}] {msg}\n")
        except Exception as e:
            print(f"[ERROR] Could not write to debug log: {e}")
    
    query_gene = best_hit['gene_id']
    write_debug(f"Starting build_gpr_rule for complex_id={complex_id}, query_gene={query_gene}")
    
    # Get stoichiometry for this specific complex
    stoichiometry_string = complex_to_stoichiometry.get(complex_id)
    if not stoichiometry_string:
        msg = f"[build_gpr_rule] No stoichiometry found for {complex_id}"
        write_debug(msg)
        return None
    
    # Parse stoichiometry
    stoichiometry = parse_stoichiometry_string(stoichiometry_string)
    if not stoichiometry:
        msg = f"[build_gpr_rule] Could not parse stoichiometry for {complex_id}: {stoichiometry_string}"
        write_debug(msg)
        return None
    
    # Filter out non-protein entries (CHEBI compounds, etc.)
    protein_stoichiometry = {uid: count for uid, count in stoichiometry.items()
                            if not uid.startswith('CHEBI:') and not uid.startswith('chebi:')}
    
    if not protein_stoichiometry:
        msg = f"[build_gpr_rule] No protein members in {complex_id} (only CHEBI entries)"
        write_debug(msg)
        return None
    
    msg = (f"[build_gpr_rule] {complex_id}: {len(protein_stoichiometry)} protein members "
           f"(filtered {len(stoichiometry) - len(protein_stoichiometry)} non-protein entries)")
    write_debug(msg)
    
    # Extract sequences for all proteins in the complex
    sequences = extract_sequences_from_complexes_db(
        list(protein_stoichiometry.keys()), 
        complexes_fasta_path
    )
    
    if not sequences:
        msg = (f"[build_gpr_rule] No sequences found in FASTA for {complex_id} proteins: "
               f"{list(protein_stoichiometry.keys())}")
        write_debug(msg)
        return None
    
    missing_seqs = set(protein_stoichiometry.keys()) - set(sequences.keys())
    if missing_seqs:
        msg = f"[build_gpr_rule] {complex_id}: {len(missing_seqs)} proteins NOT found in FASTA: {missing_seqs}"
        write_debug(msg)
    msg = f"[build_gpr_rule] {complex_id}: extracted {len(sequences)}/{len(protein_stoichiometry)} sequences"
    write_debug(msg)
    
    # DIAMOND each complex protein against genome
    genome_matches = blast_sequences_against_genome(
        sequences, 
        genomic_blast_db, 
        evalue_threshold,
        existing_genes=existing_genes
    )
    
    msg = (f"[build_gpr_rule] {complex_id}: DIAMOND found genome matches for "
           f"{len(genome_matches)}/{len(sequences)} proteins: {genome_matches}")
    write_debug(msg)
    
    # Build GPR rule with stoichiometry
    gpr_components = []
    new_genes = []
    
    for uniprot_id, stoich_value in protein_stoichiometry.items():
        # Get best genome match for this complex protein
        if uniprot_id in genome_matches and genome_matches[uniprot_id]:
            genome_gene = genome_matches[uniprot_id][0]['genome_gene']  # First/best match
            
            # Do NOT add stoichiometry notation to gene names
            # COBRApy doesn't support "gene(4)" syntax in GPR rules
            gpr_components.append(genome_gene)
            
            # Track new genes not already in metabolic model
            if genome_gene not in existing_genes:
                new_genes.append(genome_gene)
        else:
            msg = f"[build_gpr_rule] {complex_id}: no genome match for {uniprot_id}"
            write_debug(msg)
    
    if not gpr_components:
        msg = f"[build_gpr_rule] {complex_id}: no genome matches found for any protein, rule is empty"
        write_debug(msg)
        return None
    
    # Deduplicate gene names (same genome gene may match multiple complex proteins)
    seen = set()
    unique_components = []
    for g in gpr_components:
        if g not in seen:
            seen.add(g)
            unique_components.append(g)
    gpr_components = unique_components
    
    # Discard complex if it has >3 proteins and less than half matched
    total_proteins = len(protein_stoichiometry)
    matched_proteins = len(gpr_components)
    if total_proteins > 3 and matched_proteins < total_proteins / 2:
        msg = (f"[build_gpr_rule] {complex_id}: DISCARD – only {matched_proteins}/{total_proteins} "
               f"proteins matched (need at least half for complexes with >3 proteins)")
        write_debug(msg)
        print(msg)
        return None

    # Discard if the query gene (the original model gene) is not in the final rule.
    # This happens when a different genome gene scores better for the same complex
    # member, replacing the original gene.  The complex should not be assigned to
    # reactions of a gene that is no longer part of the rule.
    if query_gene and query_gene not in gpr_components:
        msg = (f"[build_gpr_rule] {complex_id}: DISCARD – query gene {query_gene} "
               f"is not in the final rule components {gpr_components}")
        write_debug(msg)
        print(msg)
        return None
    
    # Build a mapping of uniprot_id → genome_gene for downstream use
    member_gene_map = {}
    for uniprot_id in protein_stoichiometry:
        if uniprot_id in genome_matches and genome_matches[uniprot_id]:
            member_gene_map[uniprot_id] = genome_matches[uniprot_id][0]['genome_gene']

    msg = f"[build_gpr_rule] {complex_id}: GPR rule = {' and '.join(gpr_components)}"
    write_debug(msg)
    write_debug(f"SUCCESS: build_gpr_rule completed for {complex_id}\n")
    
    return {
        'gpr_rule': ' and '.join(gpr_components),
        'new_genes': new_genes,
        'complex_id': complex_id,
        'query_gene': query_gene,
        'stoichiometry': stoichiometry,  # Keep original stoichiometry (incl CHEBI) for display
        'member_gene_map': member_gene_map,  # {uniprot_id: genome_gene}
    }


def extract_protein_sequences(genome_fasta, gene_ids):
    """
    Extract protein sequences for genes from genome FASTA
    If model doesn't have protein sequences, use nucleotide translation
    """
    query_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
    
    try:
        # Parse genome FASTA and extract matching sequences
        sequence_count = 0
        for record in SeqIO.parse(genome_fasta, "fasta"):
            # Try to match gene ID
            gene_id = record.id.split()[0]  # Use first part before space
            
            if gene_id in gene_ids or any(gene_id in gid for gid in gene_ids):
                # Translate if nucleotide
                seq = record.seq
                if len(seq) % 3 == 0 and all(str(c) in 'ATGCN' for c in str(seq).upper()):
                    # Translate nucleotide to protein
                    from Bio.Seq import Seq
                    seq = Seq(str(seq)).translate(to_stop=True)
                
                query_file.write(f">{gene_id}\n{seq}\n")
                sequence_count += 1
        
        query_file.close()
        return query_file.name, sequence_count
    
    except Exception as e:
        query_file.close()
        raise Exception(f"Error extracting sequences: {str(e)}")


def cleanup_blast_db(db_dir):
    """Clean up temporary BLAST database files"""
    try:
        import shutil
        if os.path.exists(db_dir):
            shutil.rmtree(db_dir)
    except:
        pass


def generate_gpr_summary(blast_results, gene_ids, gpr_rules=None):
    """
    Generate summary of GPR analysis
    Includes GPR rule generation for genes with complex assignments
    
    Args:
        blast_results: List of BLAST results
        gene_ids: List of original gene IDs from model
        gpr_rules: Optional dict mapping gene_id to GPR rule
    
    Returns dict with analysis summary
    """
    # Get genes that have BLAST hits
    matched_genes = set(r['gene_id'] for r in blast_results)
    unmatched_genes = set(gene_ids) - matched_genes
    
    summary = {
        'total_genes': len(gene_ids),
        'matched_genes': len(matched_genes),
        'unmatched_genes': len(unmatched_genes),
        'match_percentage': (len(matched_genes) / len(gene_ids) * 100) if gene_ids else 0,
        'unmatched_list': list(unmatched_genes),
        'results_by_gene': {},
        'gpr_improvements': {}  # Gene -> improved GPR rule
    }
    
    # Group results by gene
    for result in blast_results:
        gene_id = result['gene_id']
        if gene_id not in summary['results_by_gene']:
            summary['results_by_gene'][gene_id] = []
        summary['results_by_gene'][gene_id].append(result)
    
    # Add GPR rule improvements if provided
    if gpr_rules:
        summary['gpr_improvements'] = gpr_rules
    
    return summary


def analyze_gpr_with_complexes(query_genome_fasta, complex_to_stoichiometry,
                                existing_model_genes, complexes_fasta_path=None,
                                evalue_threshold=1e-5, job=None):
    """
    Comprehensive GPR analysis workflow:
    1. Create BLAST database from query genome
    2. BLAST complexes database against genome
    3. For each best hit, find ALL complexes containing that UniProt protein
    4. For each complex, BLAST all members and build GPR rules with stoichiometry
    
    Args:
        query_genome_fasta: Path to genome FASTA file
        complex_to_stoichiometry: Dict mapping {complex_id: stoichiometry_string}
                                 Format: 'P15790(2)|P25368(0)|P38930(1)|...'
        existing_model_genes: Set/list of gene IDs already in metabolic model
        complexes_fasta_path: Path to complexes FASTA (default: data/complexes_blast_db/complexes.fasta)
        evalue_threshold: E-value threshold for filtering BLAST hits
        job: GPRAnalysisJob instance for progress updates (optional)
    
    Returns:
        {
            'gpr_rules': {(gene_id, complex_id): gpr_rule_dict},
            'new_genes_found': [list of genes],
            'summary': summary_dict
        }
    """
    def update_job_message(msg):
        """Update job message if job instance provided"""
        if job:
            job.message = msg
            job.save(update_fields=['message'])
            print(f"[gpr-worker-msg] {msg}")
    
    # Path to prebuilt DIAMOND database for complexes
    complexes_db_path = os.path.join(
        settings.BASE_DIR,
        'data', 'complexes_blast_db', 'complexes.dmnd'
    )
    
    if complexes_fasta_path is None:
        complexes_fasta_path = os.path.join(
            settings.BASE_DIR,
            'data', 'complexes_blast_db', 'complexes.fasta'
        )
    
    # Step 1: Extract ONLY model genes from the genome FASTA file
    update_job_message("Extracting model genes from genome FASTA...")
    try:
        model_genes_fasta, genes_found = extract_model_genes_from_fasta(
            query_genome_fasta,
            existing_model_genes
        )
        update_job_message(f"Found {genes_found}/{len(existing_model_genes)} model genes in genome FASTA file")
    except Exception as e:
        raise Exception(f"Failed to extract model genes from genome: {str(e)}")
    
    temp_genome_fasta_path = model_genes_fasta
    
    # Step 2: Search model genes against PREBUILT complexes database using DIAMOND (much faster!)
    update_job_message("Searching model genes against complexes database (DIAMOND)...")
    blast_results = run_diamond_search(temp_genome_fasta_path, complexes_db_path, evalue_threshold)
    update_job_message(f"Search complete: {len(blast_results)} matches found")
    
    # Step 2b: Create genome database for searching complex members (using original full genome)
    update_job_message("Creating genome database for complex member search...")
    genome_db_path, genome_db_dir = create_diamond_db(query_genome_fasta, 'genome_db')
    
    # Clean up temp genome FASTA since we have the database now
    if os.path.exists(temp_genome_fasta_path):
        os.unlink(temp_genome_fasta_path)
    
    try:
        # Step 3: Filter results to get best hits per query gene
        update_job_message("Filtering results by E-value...")
        best_hits = get_best_blast_hit(blast_results, evalue_threshold)
        update_job_message(f"Found {len(best_hits)} best hits after filtering")
        
        # Step 4: For each best hit, find ALL complexes containing that UniProt
        update_job_message("Analyzing complex assignments...")
        gpr_rules = {}
        all_new_genes = set()
        
        for idx, (genome_gene, best_hit) in enumerate(best_hits.items(), 1):
            uniprot_id = best_hit['uniprot_id']  # Complex protein that matched this genome gene
            update_job_message(f"Processing gene {idx}/{len(best_hits)}: {genome_gene}...")
            
            # Find ALL complexes that contain this UniProt ID
            complexes_containing_uniprot = find_all_complexes_for_uniprot(
                uniprot_id, 
                complex_to_stoichiometry
            )
            
            # Step 5: For each complex, build GPR rule
            for complex_id in complexes_containing_uniprot:
                try:
                    # Create a hit object with the genome gene and matched complex protein
                    hit_for_complex = {
                        'gene_id': genome_gene,  # The genome gene we found
                        'uniprot_id': uniprot_id,  # The complex protein it matched
                        'evalue': best_hit['evalue'],
                        'bitscore': best_hit['bitscore'],
                        'identity': best_hit['identity'],
                        'query_coverage': best_hit['query_coverage']
                    }
                    
                    gpr_result = build_gpr_rule_from_complex(
                        hit_for_complex,
                        complex_to_stoichiometry,
                        genome_db_path,
                        existing_model_genes,
                        complexes_fasta_path,
                        complex_id,
                        evalue_threshold
                    )
                    
                    if gpr_result:
                        # Use (gene_id, complex_id) as key to handle multiple assignments
                        key = (genome_gene, complex_id)
                        gpr_rules[key] = gpr_result
                        all_new_genes.update(gpr_result['new_genes'])
                
                except Exception as e:
                    import traceback
                    print(f"[gpr-error] Failed to build GPR rule for gene={genome_gene}, complex={complex_id}: {e}")
                    traceback.print_exc()
        
        # Step 6: Generate summary
        update_job_message("Generating analysis summary...")
        summary = generate_gpr_summary(blast_results, existing_model_genes, gpr_rules)
        update_job_message(f"Analysis complete: {len(gpr_rules)} GPR rules generated")
        
        return {
            'gpr_rules': gpr_rules,
            'new_genes_found': list(all_new_genes),
            'summary': summary,
            'best_hits': best_hits,  # Best hit per gene
            'blast_results': blast_results,  # All BLAST results
            'genome_db_path': genome_db_path,
            'genome_db_dir': genome_db_dir
        }
    
    finally:
        # Note: cleanup_blast_db is called by the worker after populate_blast_results_for_complexes
        pass


def format_gpr_rule_text(genome_genes_dict, stoichiometry_dict):
    """
    Build a GPR rule text from a set of genome genes belonging to a complex.
    
    Args:
        genome_genes_dict: Dict mapping genome_gene → {uniprot_id, stoichiometry_count, ...}
        stoichiometry_dict: Dict mapping uniprot_id → stoichiometry_count
    
    Returns:
        GPR rule string, e.g. "(gene_1 AND gene_2(2))" with stoichiometry
    
    Example:
        Input: {
            'gene_1': {'uniprot_id': 'P12345', 'stoichiometry': 1},
            'gene_2': {'uniprot_id': 'P67890', 'stoichiometry': 2}
        }
        Output: "(gene_1 AND gene_2(2))"
    """
    if not genome_genes_dict:
        return ""
    
    gene_rules = []
    for genome_gene, gene_info in genome_genes_dict.items():
        uniprot_id = gene_info.get('uniprot_id', '')
        stoich_count = stoichiometry_dict.get(uniprot_id, 1)
        
        # Note: DO NOT add stoichiometry notation to gene names
        # COBRApy doesn't support "gene(4)" syntax in GPR rules
        # Just use the gene name - stoichiometry is handled at reaction level
        gene_rules.append(genome_gene)
    
    # Join with AND
    gpr_text = " AND ".join(gene_rules)
    return f"({gpr_text})" if gpr_text else ""


def create_gpr_rules_for_reactions(reactions_with_gene, selected_complex_info):
    """
    Create GPR rules for all reactions containing the selected complex's genes.
    
    Args:
        reactions_with_gene: List of reaction info dicts with 'reaction_id', 'gpr', etc
        selected_complex_info: Dict with 'genome_genes', 'stoichiometry', 'complex_id'
    
    Returns:
        Dict mapping reaction_id → new_gpr_rule
    
    Example output:
        {
            'R_REACTION1': '(gene_1 AND gene_2(2))',
            'R_REACTION2': '(gene_1 AND gene_2(2))'
        }
    """
    gpr_rules = {}
    
    if not reactions_with_gene:
        return gpr_rules
    
    # Get the genome genes and stoichiometry from selected complex
    genome_genes_dict = selected_complex_info.get('genome_genes', {})
    stoichiometry_dict = selected_complex_info.get('stoichiometry', {})
    
    # Build the GPR rule for this complex
    complex_gpr_rule = format_gpr_rule_text(genome_genes_dict, stoichiometry_dict)
    
    if not complex_gpr_rule:
        return gpr_rules
    
    # Apply this rule to all reactions containing genes from this complex
    for reaction in reactions_with_gene:
        reaction_id = reaction['reaction_id']
        # Use the complex's GPR rule for this reaction
        # (In real scenario, you might want to combine with existing GPR)
        gpr_rules[reaction_id] = complex_gpr_rule
    
    return gpr_rules


def update_model_with_gpr_rules(model_path, gpr_rules_dict, output_path, merge=True):
    """
    Update metabolic model with new GPR rules for reactions.
    
    Args:
        model_path: Path to original metabolic model (JSON or SBML)
        gpr_rules_dict: Dict mapping reaction_id → new_gpr_rule_text
        output_path: Where to save the updated model
        merge: If True, merge new rules with existing old rules (preserve
               non-overlapping 'or' branches). If False, set rules verbatim.
    
    Returns:
        {'success': True/False, 'path': output_path, 'updates': count}
    
    Supported formats: JSON (.json) and SBML (.xml, .sbml)
    """
    try:
        import re
        updated_count = 0
        
        # Normalize GPR rules: convert uppercase operators to lowercase for COBRApy compatibility
        def normalize_gpr_rule(rule):
            """Normalize GPR rule to use lowercase and/or operators"""
            if not rule:
                return rule
            # Use regex to replace AND/OR with and/or, handling any amount of whitespace
            # This correctly handles "AND", " AND ", "AND ", etc.
            normalized = re.sub(r'\bAND\b', 'and', rule, flags=re.IGNORECASE)
            normalized = re.sub(r'\bOR\b', 'or', normalized, flags=re.IGNORECASE)
            # Clean up multiple spaces
            normalized = re.sub(r'\s+', ' ', normalized)
            return normalized.strip()
        
        def strip_stoichiometry(rule):
            """Remove (N) stoichiometry notation from gene names.
            COBRApy parses GPR rules as Python AST, so gene(4) is
            interpreted as a function call and crashes."""
            return re.sub(r'(\b[A-Za-z0-9_.]+)\(\d+\)', r'\1', rule) if rule else rule

        def _genes_overlap(old_rule, new_rule):
            """Return True if at least one gene in new_rule also appears in old_rule.
            If the old rule is empty, always allow the new rule."""
            if not old_rule:
                return True
            old_genes = set(parse_gpr_rule(old_rule))
            new_genes = set(parse_gpr_rule(new_rule))
            return bool(old_genes & new_genes)

        def _split_top_level_or(rule):
            """Split a GPR rule at top-level 'or' operators, respecting parentheses."""
            if not rule:
                return []
            branches = []
            depth = 0
            current = []
            for token in rule.split():
                depth += token.count('(') - token.count(')')
                if token.lower() == 'or' and depth == 0:
                    branches.append(' '.join(current))
                    current = []
                else:
                    current.append(token)
            if current:
                branches.append(' '.join(current))
            return [b.strip() for b in branches if b.strip()]

        def _merge_with_old_rule(old_rule, new_rule):
            """Merge a new GPR rule into the old rule.
            Branches (top-level 'or' parts) of the old rule whose genes
            overlap with the new rule are replaced; non-overlapping
            branches are preserved and combined with 'or'.
            Example: old='A or B', new='(A and C)' → '(A and C) or B'"""
            if not old_rule:
                return new_rule
            new_genes = set(parse_gpr_rule(new_rule))
            old_branches = _split_top_level_or(old_rule)
            kept = []
            for branch in old_branches:
                branch_genes = set(parse_gpr_rule(branch))
                if not (branch_genes & new_genes):
                    kept.append(branch)
            all_parts = [new_rule] + kept
            merged = ' or '.join(all_parts)
            print(f"[update_model] merge: old='{old_rule}' + new='{new_rule}' → '{merged}'")
            return merged

        # Normalize all rules first
        normalized_rules = {rxn_id: strip_stoichiometry(normalize_gpr_rule(rule))
                           for rxn_id, rule in gpr_rules_dict.items()}
        
        # Track the actual merged rules that end up in the model
        merged_rules = {}
        
        # Debug: log all rules before processing
        for rxn_id, rule in normalized_rules.items():
            if rule:
                print(f"[update_model] {rxn_id}: {rule}")
        
        if model_path.endswith('.json'):
            # Update JSON format model
            with open(model_path, 'r') as f:
                model_data = json.load(f)
            
            if 'reactions' in model_data:
                for reaction in model_data['reactions']:
                    reaction_id = reaction.get('id')
                    if reaction_id in normalized_rules:
                        old_gpr = reaction.get('gpr', '')
                        new_gpr = normalized_rules[reaction_id]
                        if merge:
                            if not _genes_overlap(old_gpr, new_gpr):
                                print(f"[update_model] SKIP {reaction_id}: no gene overlap between old and new rule")
                                continue
                            final_rule = _merge_with_old_rule(old_gpr, new_gpr)
                        else:
                            final_rule = new_gpr
                        reaction['gpr'] = final_rule
                        merged_rules[reaction_id] = final_rule
                        updated_count += 1
            
            # Write updated model
            with open(output_path, 'w') as f:
                json.dump(model_data, f, indent=2)
        
        elif model_path.endswith('.xml') or model_path.endswith('.sbml'):
            # Update SBML format model using COBRApy
            try:
                import cobra
                model = cobra.io.read_sbml_model(model_path)
                
                for reaction in model.reactions:
                    if reaction.id in normalized_rules:
                        rule_to_set = normalized_rules[reaction.id]
                        old_gpr = str(reaction.gpr) if reaction.gpr else ''
                        if merge:
                            if not _genes_overlap(old_gpr, rule_to_set):
                                print(f"[update_model] SKIP {reaction.id}: no gene overlap between old and new rule")
                                continue
                            final_rule = _merge_with_old_rule(old_gpr, rule_to_set)
                        else:
                            final_rule = rule_to_set
                        try:
                            print(f"[update_model] Setting {reaction.id} to: {final_rule}")
                            reaction.gene_reaction_rule = final_rule
                            merged_rules[reaction.id] = final_rule
                            updated_count += 1
                        except Exception as e:
                            print(f"[update_model] ERROR setting {reaction.id}: {type(e).__name__}: {str(e)}")
                            print(f"[update_model] Rule was: {final_rule}")
                            raise
                
                cobra.io.write_sbml_model(model, output_path)
            except ImportError:
                # Fallback: XML-based update without COBRApy
                with open(model_path, 'r') as f:
                    content = f.read()
                
                # This is a simplified approach - proper SBML parsing recommended
                import re
                for reaction_id, gpr_rule in gpr_rules_dict.items():
                    # Find reaction and update its gpr attribute
                    pattern = f'<reaction id="{reaction_id}"[^>]*>'
                    # This is limited - COBRApy is recommended
                    updated_count += 1
                
                with open(output_path, 'w') as f:
                    f.write(content)
        
        return {
            'success': True,
            'path': output_path,
            'updates': updated_count,
            'merged_rules': merged_rules
        }
    
    except Exception as e:
        return {
            'success': False,
            'error': str(e),
            'updates': 0
        }


def build_gpr_rules_from_unique_cases(unique_cases, complex_to_stoichiometry, gpr_blast_results=None, best_hits=None, prebuilt_gpr_rules=None):
    """
    Build complete GPR rules from unique cases grouped by complex.

    Uses pre-built GPR rules from analyze_gpr_with_complexes() which already
    DIAMOND-searched each complex member against the genome.  Falls back to
    best-hits cross-referencing only when pre-built rules are unavailable.

    Args:
        unique_cases: List of unique case dicts (from collect_ambiguous_cases)
        complex_to_stoichiometry: Dict {complex_id: stoichiometry_string}
        gpr_blast_results: (unused, kept for API compat)
        best_hits: Optional dict {gene_id: hit_dict} — fallback cross-reference
        prebuilt_gpr_rules: Dict {(genome_gene, complex_id): gpr_result_dict}
                            from analyze_gpr_with_complexes().  Each value has
                            key 'gpr_rule' with the AND-rule string.

    Returns:
        Dict mapping reaction_id → gpr_rule_text
    """
    gpr_rules_by_reaction = {}
    if prebuilt_gpr_rules is None:
        prebuilt_gpr_rules = {}

    # ---- Phase 1: collect which complex applies per gene ----
    gene_to_complex = {}   # {genome_gene: complex_id}
    complex_ids_seen = set()

    for unique_case in unique_cases:
        genome_gene = unique_case.get('genome_gene')
        if not unique_case.get('possible_complexes'):
            continue
        complex_id = unique_case['possible_complexes'][0].get('complex_id')
        gene_to_complex[genome_gene] = complex_id
        complex_ids_seen.add(complex_id)

    # ---- Phase 2: build per-complex GPR rule ----
    complex_to_gpr_rule = {}
    for complex_id in complex_ids_seen:
        # Try to find a pre-built rule for this complex (any gene that matched it)
        prebuilt = None
        for (g, cid), rule_dict in prebuilt_gpr_rules.items():
            if cid == complex_id and rule_dict:
                prebuilt = rule_dict
                break

        if prebuilt and prebuilt.get('gpr_rule'):
            rule_text = re.sub(r'\bAND\b', 'and', prebuilt['gpr_rule'], flags=re.IGNORECASE)
            rule_text = re.sub(r'\bOR\b', 'or', rule_text, flags=re.IGNORECASE)
            # Wrap in parens if it contains 'and' and isn't already wrapped
            if ' and ' in rule_text and not rule_text.startswith('('):
                rule_text = f"({rule_text})"
            complex_to_gpr_rule[complex_id] = rule_text
            print(f"[build_gpr_rules] Complex {complex_id}: using pre-built rule: {rule_text}")
        else:
            # Fallback: reconstruct from best_hits cross-reference
            stoich_str = complex_to_stoichiometry.get(complex_id, '')
            stoich_dict = parse_stoichiometry_string(stoich_str) if stoich_str else {}

            # Build uniprot→gene map from best_hits
            uniprot_to_gene = {}
            if best_hits:
                for gene_id, hit in best_hits.items():
                    uid = hit.get('uniprot_id', '')
                    if uid:
                        uniprot_to_gene[normalize_uniprot_id(uid).upper()] = gene_id
            # Also from cases themselves
            for case in unique_cases:
                uid = case.get('uniprot_id', '')
                gene = case.get('genome_gene', '')
                if uid and gene:
                    uniprot_to_gene[normalize_uniprot_id(uid).upper()] = gene

            gene_parts = []
            for uid, count in stoich_dict.items():
                # Skip non-protein entries (CHEBI compounds, etc.)
                if uid.startswith('CHEBI:') or uid.startswith('chebi:'):
                    continue
                mapped = uniprot_to_gene.get(uid.upper(), '')
                if mapped:
                    # Note: DO NOT add stoichiometry notation (count) to gene names
                    # COBRApy doesn't support "gene(4)" syntax in GPR rules
                    # Stoichiometry is handled at the reaction level, not in GPR
                    gene_parts.append(mapped)
            if gene_parts:
                rule = ' and '.join(gene_parts)
                rule = f"({rule})"
                complex_to_gpr_rule[complex_id] = rule
                print(f"[build_gpr_rules] Complex {complex_id}: fallback rule: {rule}")

    # ---- Phase 3: map reactions → complex rules (OR when multiple) ----
    # Only assign a complex rule to a reaction if the rule still contains the
    # original gene that linked this complex to the reaction.  When a different
    # gene in the genome scores better for the same complex member, the
    # original gene may be absent from the rule → skip that assignment.
    reaction_to_complex_rules = {}  # {reaction_id: {complex_id: gpr_rule}}
    for unique_case in unique_cases:
        reactions = unique_case.get('reactions_with_gene', [])
        genome_gene = unique_case.get('genome_gene', '')
        complex_id = unique_case['possible_complexes'][0]['complex_id'] if unique_case.get('possible_complexes') else None
        if complex_id and complex_id in complex_to_gpr_rule:
            gpr_rule = complex_to_gpr_rule[complex_id]
            # Check that the rule contains the gene that linked the complex to
            # this reaction; if not, a different gene replaced it → skip
            if genome_gene and genome_gene not in gpr_rule:
                print(f"[build_gpr_rules] Skipping complex {complex_id} for gene {genome_gene}: "
                      f"rule '{gpr_rule}' no longer contains the original gene")
                continue
            for reaction in reactions:
                reaction_id = reaction['reaction_id']
                if reaction_id not in reaction_to_complex_rules:
                    reaction_to_complex_rules[reaction_id] = {}
                reaction_to_complex_rules[reaction_id][complex_id] = gpr_rule

    for reaction_id, complex_rules in reaction_to_complex_rules.items():
        distinct_rules = list(dict.fromkeys(complex_rules.values()))
        gpr_rules_by_reaction[reaction_id] = ' or '.join(distinct_rules)

    return gpr_rules_by_reaction


def build_complex_summary(cases, complex_to_stoichiometry, gpr_rules, best_hits=None, complex_descriptions=None, prebuilt_gpr_rules=None, blast_results=None):
    """
    Build a summary table GROUPED BY REACTION showing all genes and complexes involved.

    Args:
        cases: List of case dicts (from collect_ambiguous_cases, unique or resolved)
        complex_to_stoichiometry: Dict {complex_id: stoichiometry_string}
        gpr_rules: Dict {reaction_id: gpr_rule_text}
        best_hits: Optional dict {gene_id: hit_dict} for looking up other gene mappings
        complex_descriptions: Optional dict {complex_id: description_text}
        prebuilt_gpr_rules: Optional dict {(gene, complex_id): rule_dict} from
                            analyze_gpr_with_complexes(), used to extract the
                            genome gene assigned to each complex member.
        blast_results: Optional list of BLAST result objects with gene_id, complex_id, evalue

    Returns:
        List of dicts grouped by reaction, with keys:
        - reaction_id
        - reaction_name
        - gpr_rule
        - complexes_in_reaction: [
            {
              'complex_id': cpx_id,
              'complex_description': description,
              'complex_members': [{'uniprot': uid, 'genome_gene': gene, 'evalue': e_value}]
            }
          ]
    """
    if complex_descriptions is None:
        complex_descriptions = {}
    if prebuilt_gpr_rules is None:
        prebuilt_gpr_rules = {}

    # Build e-value map from BLAST results: {(gene_id, complex_id): evalue}
    # Note: in BLAST results, complex_id field actually stores the UniProt ID
    evalue_map = {}
    if blast_results:
        for br in blast_results:
            gene_id = br.gene_id if hasattr(br, 'gene_id') else br.get('gene_id', '')
            uniprot_id = br.complex_id if hasattr(br, 'complex_id') else br.get('complex_id', '')
            evalue = br.evalue if hasattr(br, 'evalue') else br.get('evalue')
            if gene_id and uniprot_id:
                uid_norm = normalize_uniprot_id(uniprot_id).upper()
                evalue_map[(gene_id, uid_norm)] = evalue

    # Map UniProt → genome gene (from best_hits if available)
    uniprot_to_gene = {}
    if best_hits:
        for gene_id, hit in best_hits.items():
            uid = hit.get('uniprot_id', '')
            if uid:
                uniprot_to_gene[normalize_uniprot_id(uid).upper()] = gene_id

    # Also build from the cases themselves
    for case in cases:
        uid = case.get('uniprot_id', '')
        gene = case.get('genome_gene', '')
        if uid and gene:
            uniprot_to_gene[normalize_uniprot_id(uid).upper()] = gene

    # Extract genome-gene mappings from pre-built GPR rules
    complex_member_genes = {}  # {complex_id: {uniprot: genome_gene}}
    for (g, cid), rule_dict in prebuilt_gpr_rules.items():
        if not rule_dict:
            continue
        mgm = rule_dict.get('member_gene_map', {})
        if cid not in complex_member_genes:
            complex_member_genes[cid] = {}
        for uid, genome_gene in mgm.items():
            complex_member_genes[cid][normalize_uniprot_id(uid).upper()] = genome_gene

    # --- Organize by reaction, with complexes deduplicated ---
    # {reaction_id: {
    #   'reaction_name': name,
    #   'complexes': {
    #     complex_id: {'description', 'members': [...]}  # deduplicated by complex_id
    #   }
    # }}
    reactions_data = {}
    
    for case in cases:
        complexes = case.get('possible_complexes', [])
        reactions = case.get('reactions_with_gene', [])
        selected_complex_id = case.get('selected_complex')

        # Determine which complex applies
        if selected_complex_id:
            # User explicitly selected this complex — show only it
            target_complexes = [c for c in complexes if c['complex_id'] == selected_complex_id]
            if not target_complexes:
                target_complexes = [{'complex_id': selected_complex_id, 'stoichiometry': {}}]
        elif 'selected_complex' in case:
            # Key exists but is empty/None — user chose "None", skip entirely
            continue
        else:
            # Key not present — unique case (auto-resolved), use first complex
            target_complexes = complexes[:1]

        for cpx in target_complexes:
            cpx_id = cpx['complex_id']

            # Get stoichiometry
            stoich = cpx.get('stoichiometry', {})
            if not stoich:
                stoich_str = complex_to_stoichiometry.get(cpx_id, '')
                stoich = parse_stoichiometry_string(stoich_str) if stoich_str else {}

            # Build complex members with e-values
            cpx_genes = complex_member_genes.get(cpx_id, {})
            
            # Also extract genome gene mappings from blast_results if available
            cpx_blast_results = cpx.get('blast_results', [])
            blast_result_map = {}  # {uniprot: {genome_gene, evalue, ...}}
            for br in cpx_blast_results:
                uniprot_id = br.get('uniprot_id', '')
                genome_gene = br.get('genome_gene', '')
                evalue = br.get('evalue')
                if uniprot_id and genome_gene:
                    uid_norm = normalize_uniprot_id(uniprot_id).upper()
                    if uid_norm not in blast_result_map:
                        blast_result_map[uid_norm] = {
                            'genome_gene': genome_gene,
                            'evalue': evalue,
                        }
                    # Keep the best (lowest) e-value hit
                    elif evalue and (not blast_result_map[uid_norm].get('evalue') or evalue < blast_result_map[uid_norm]['evalue']):
                        blast_result_map[uid_norm] = {
                            'genome_gene': genome_gene,
                            'evalue': evalue,
                        }
            
            complex_members = []
            for uid, count in stoich.items():
                uid_norm = uid.upper()
                # Priority: blast_results > prebuilt_gpr_rules > uniprot_to_gene
                if uid_norm in blast_result_map:
                    mapped_gene = blast_result_map[uid_norm]['genome_gene']
                    member_evalue = blast_result_map[uid_norm]['evalue']
                else:
                    mapped_gene = cpx_genes.get(uid_norm, '') or uniprot_to_gene.get(uid_norm, '')
                    member_evalue = evalue_map.get((mapped_gene, uid_norm)) if mapped_gene else None
                complex_members.append({
                    'uniprot': uid,
                    'genome_gene': mapped_gene,
                    'evalue': member_evalue,
                })

            # For each reaction, store this complex (deduplicated)
            for reaction in reactions:
                rxn_id = reaction.get('reaction_id', '')
                rxn_name = reaction.get('name', rxn_id)
                
                # Initialize reaction if not yet seen
                if rxn_id not in reactions_data:
                    reactions_data[rxn_id] = {
                        'reaction_name': rxn_name,
                        'original_gpr': reaction.get('gpr', ''),
                        'complexes': {}
                    }
                
                # Add/update complex for this reaction (deduplicated by complex_id)
                if cpx_id not in reactions_data[rxn_id]['complexes']:
                    reactions_data[rxn_id]['complexes'][cpx_id] = {
                        'complex_id': cpx_id,
                        'complex_description': complex_descriptions.get(cpx_id, ''),
                        'complex_members': complex_members,
                    }

    # --- Convert to output format (one row per reaction, sorted by reaction_id) ---
    rows = []
    for rxn_id in sorted(reactions_data.keys()):
        rxn_info = reactions_data[rxn_id]
        gpr_rule = gpr_rules.get(rxn_id, '')
        
        # Build complexes_in_reaction list (deduplicated)
        complexes_in_reaction = sorted(
            rxn_info['complexes'].values(),
            key=lambda c: c['complex_id']
        )
        
        rows.append({
            'reaction_id': rxn_id,
            'reaction_name': rxn_info['reaction_name'],
            'original_gpr': rxn_info.get('original_gpr', ''),
            'gpr_rule': gpr_rule,
            'complexes_in_reaction': complexes_in_reaction,
        })

    return rows

