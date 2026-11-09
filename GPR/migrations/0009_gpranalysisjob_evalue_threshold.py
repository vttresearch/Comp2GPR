# Generated migration for evalue_threshold field

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('GPR', '0008_alter_gprambiguouscaseresolution_unique_together'),
    ]

    operations = [
        migrations.AddField(
            model_name='gpranalysisjob',
            name='evalue_threshold',
            field=models.FloatField(default=1e-30, help_text='E-value threshold for BLAST hits (lower = stricter)'),
        ),
    ]
