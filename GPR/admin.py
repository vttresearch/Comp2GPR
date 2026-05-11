from django.contrib import admin
from .models import GPRAnalysisJob, GPRBlastResult, ComplexMetadata, GPRAmbiguousCaseResolution


@admin.register(GPRAnalysisJob)
class GPRAnalysisJobAdmin(admin.ModelAdmin):
    list_display = ('job_id', 'user', 'organism_name', 'status', 'created_at', 'gene_count', 'matched_genes')
    list_filter = ('status', 'created_at')
    search_fields = ('job_id', 'organism_name', 'user__username')
    readonly_fields = ('job_id', 'created_at', 'updated_at')


@admin.register(GPRBlastResult)
class GPRBlastResultAdmin(admin.ModelAdmin):
    list_display = ('gene_id', 'complex_id', 'bitscore', 'evalue', 'identity', 'query_coverage')
    list_filter = ('bitscore', 'identity')
    search_fields = ('gene_id', 'complex_id')
    readonly_fields = ('gpr_job', 'gene_id', 'complex_id')


@admin.register(ComplexMetadata)
class ComplexMetadataAdmin(admin.ModelAdmin):
    list_display = ('complex_id', 'name')
    search_fields = ('complex_id', 'name', 'aliases')
    fields = ('complex_id', 'name', 'aliases', 'description')


@admin.register(GPRAmbiguousCaseResolution)
class GPRAmbiguousCaseResolutionAdmin(admin.ModelAdmin):
    list_display = ('gpr_job', 'genome_gene', 'reaction_id', 'selected_complex', 'timestamp')
    list_filter = ('gpr_job', 'timestamp')
    search_fields = ('genome_gene', 'reaction_id', 'selected_complex')
    readonly_fields = ('gpr_job', 'genome_gene', 'reaction_id', 'timestamp')

