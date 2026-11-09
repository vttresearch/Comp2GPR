from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('GPR', '0009_gpranalysisjob_evalue_threshold'),
    ]

    operations = [
        migrations.AddField(
            model_name='gpranalysisjob',
            name='needs_finalization',
            field=models.BooleanField(default=False, help_text='Awaiting background worker to apply resolved GPR rules to the model'),
        ),
    ]
