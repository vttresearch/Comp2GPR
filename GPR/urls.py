from django.urls import path
from . import views

app_name = 'gpr'

urlpatterns = [
    path('', views.gpr_analysis_index, name='index'),
    path('job/<int:job_id>/', views.gpr_job_detail, name='job_detail'),
    path('job/<int:job_id>/status/', views.gpr_job_status, name='job_status'),
    path('job/<int:job_id>/download/', views.gpr_results_download, name='results_download'),
    path('job/<int:job_id>/resolve-ambiguities/', views.gpr_resolve_ambiguities, name='resolve_ambiguities'),
    path('job/<int:job_id>/user-selections/', views.gpr_user_selections, name='user_selections'),
    path('job/<int:job_id>/edit-gpr-rule/', views.edit_gpr_rule, name='edit_gpr_rule'),
]
