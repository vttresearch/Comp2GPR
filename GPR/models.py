from django.db import models
from django.contrib.auth.models import User

class GPRAnalysisJob(models.Model):
    """Store GPR analysis jobs"""
    STATUS_CHOICES = [
        ('queued', 'Queued'),
        ('processing', 'Processing'),
        ('ambiguity_review', 'Awaiting User Resolution'),
        ('completed', 'Completed'),
        ('error', 'Error'),
    ]
    
    user = models.ForeignKey(User, on_delete=models.CASCADE, related_name='gpr_jobs', null=True, blank=True)
    session_key = models.CharField(max_length=40, blank=True, null=True, db_index=True)
    job_id = models.CharField(max_length=100, unique=True)
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default='queued')
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    
    # Input files
    metabolic_model = models.FileField(upload_to='gpr_models/')
    genome_fasta = models.FileField(upload_to='gpr_genomes/')
    
    # Results
    results_file = models.FileField(upload_to='gpr_results/', null=True, blank=True)
    error_message = models.TextField(null=True, blank=True)
    message = models.TextField(null=True, blank=True, help_text="Current processing status message")
    
    # Metadata
    organism_name = models.CharField(max_length=255, null=True, blank=True)
    gene_count = models.IntegerField(null=True, blank=True)
    matched_genes = models.IntegerField(default=0)
    
    # Store ambiguous cases as JSON for review
    ambiguous_cases_json = models.TextField(null=True, blank=True)
    
    # Track if ambiguous cases need user resolution (independent of job status)
    has_unresolved_ambiguities = models.BooleanField(default=False, help_text="Indicates if ambiguous cases need user review (independent of job status)")
    
    # GPR rules and updated model
    gpr_rules_json = models.TextField(null=True, blank=True, help_text="Mapping of reaction_id to GPR rule text")
    updated_model_path = models.CharField(max_length=512, null=True, blank=True, help_text="Path to updated model with GPR rules")
    complex_summary_json = models.TextField(null=True, blank=True, help_text="Summary of gene-complex-reaction assignments for display")
    
    @property
    def queue_position(self):
        """Return 1-based position in queue, or None if not queued.
        
        Jobs stuck in 'processing' for over 2 hours are considered stale
        and excluded from the count.
        """
        if self.status != 'queued':
            return None
        from django.utils import timezone
        from datetime import timedelta
        stale_cutoff = timezone.now() - timedelta(hours=2)
        ahead = GPRAnalysisJob.objects.filter(
            status='queued',
            created_at__lt=self.created_at,
        ).count()
        # +1 if there's a recently-active processing job (not stale)
        running = GPRAnalysisJob.objects.filter(
            status='processing',
            updated_at__gte=stale_cutoff,
        ).exists()
        return ahead + (1 if running else 0) + 1

    @classmethod
    def queue_length(cls):
        """Total number of queued + actively processing jobs.
        
        Excludes jobs stuck in 'processing' for over 2 hours.
        """
        from django.utils import timezone
        from datetime import timedelta
        stale_cutoff = timezone.now() - timedelta(hours=2)
        queued = cls.objects.filter(status='queued').count()
        active = cls.objects.filter(
            status='processing',
            updated_at__gte=stale_cutoff,
        ).count()
        return queued + active

    def __str__(self):
        return f"GPR Job {self.job_id} - {self.status}"
    
    class Meta:
        ordering = ['-created_at']


class GPRBlastResult(models.Model):
    """Store individual BLAST results"""
    gpr_job = models.ForeignKey(GPRAnalysisJob, on_delete=models.CASCADE, related_name='blast_results')
    
    gene_id = models.CharField(max_length=255)
    complex_id = models.CharField(max_length=255)
    evalue = models.FloatField()
    bitscore = models.FloatField()
    identity = models.FloatField()  # Percentage identity
    query_coverage = models.FloatField()  # Query coverage percentage
    
    def __str__(self):
        return f"{self.gene_id} -> {self.complex_id}"
    
    class Meta:
        ordering = ['-bitscore']


class ComplexMetadata(models.Model):
    """Store complex names and aliases"""
    complex_id = models.CharField(max_length=255, unique=True, primary_key=True)
    name = models.CharField(max_length=500)
    aliases = models.TextField(blank=True, help_text="Comma-separated list of aliases")
    description = models.TextField(blank=True)
    
    def get_aliases_list(self):
        """Return aliases as list"""
        return [a.strip() for a in self.aliases.split(',') if a.strip()]
    
    def __str__(self):
        return f"{self.complex_id}: {self.name}"
    
    class Meta:
        verbose_name_plural = "Complex Metadata"


class GPRAmbiguousCaseResolution(models.Model):
    """Store user's resolution of ambiguous complex assignments"""
    gpr_job = models.ForeignKey(GPRAnalysisJob, on_delete=models.CASCADE, related_name='complex_resolutions')
    
    genome_gene = models.CharField(max_length=255)
    reaction_id = models.CharField(max_length=255)
    selected_complex = models.CharField(max_length=255)
    timestamp = models.DateTimeField(auto_now_add=True)
    
    class Meta:
        ordering = ['genome_gene', 'reaction_id']
    
    def __str__(self):
        return f"{self.genome_gene} in {self.reaction_id} -> {self.selected_complex}"
