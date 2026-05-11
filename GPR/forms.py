from django import forms
from .models import GPRAnalysisJob, GPRAmbiguousCaseResolution

class GPRAnalysisForm(forms.ModelForm):
    """Form for uploading metabolic model and genome FASTA"""
    organism_name = forms.CharField(
        max_length=255,
        required=True,
        label="Organism Name",
        widget=forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'e.g., Saccharomyces cerevisiae'})
    )
    
    metabolic_model = forms.FileField(
        label="Metabolic Model (SBML or JSON)",
        widget=forms.FileInput(attrs={'class': 'form-control', 'accept': '.xml,.json,.sbml'})
    )
    
    genome_fasta = forms.FileField(
        label="Genome FASTA File",
        widget=forms.FileInput(attrs={'class': 'form-control', 'accept': '.fasta,.fa,.fna'})
    )
    
    class Meta:
        model = GPRAnalysisJob
        fields = ['organism_name', 'metabolic_model', 'genome_fasta']


class ComplexSelectionForm(forms.Form):
    """Dynamic form for selecting complexes for ambiguous cases"""
    
    def __init__(self, ambiguous_cases, complex_metadata_dict, *args, **kwargs):
        """
        Args:
            ambiguous_cases: List of ambiguous case dicts with structure:
                {
                    'genome_gene': gene_id,
                    'reaction_id': reaction_id,
                    'possible_complexes': [complex_id, ...]
                }
            complex_metadata_dict: Dict {complex_id: {'name': name, 'aliases': [...]}}
        """
        super().__init__(*args, **kwargs)
        
        self.ambiguous_cases = ambiguous_cases
        self.complex_metadata = complex_metadata_dict
        
        # Create a radio field for each ambiguous case
        for idx, case in enumerate(ambiguous_cases):
            gene = case['genome_gene']
            reaction = case['reaction_id']
            complexes = case['possible_complexes']
            
            # Create field name
            field_name = f"case_{idx}"
            
            # Create choices
            choices = []
            for complex_id in complexes:
                metadata = complex_metadata_dict.get(complex_id, {})
                label = f"{complex_id}: {metadata.get('name', 'Unknown')}"
                if metadata.get('aliases'):
                    label += f" ({', '.join(metadata['aliases'])})"
                choices.append((complex_id, label))
            
            # Create radio field
            self.fields[field_name] = forms.ChoiceField(
                choices=choices,
                widget=forms.RadioSelect(attrs={'class': 'form-check-input'}),
                label=f"{gene} in Reaction {reaction}",
                initial=choices[0][0] if choices else None
            )
            
            # Store metadata for template
            self.case_info = {
                field_name: {
                    'gene': gene,
                    'reaction': reaction,
                    'complexes': complexes
                }
            }

