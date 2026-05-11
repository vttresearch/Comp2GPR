# GPR Analysis – Django App

A Django application for **Gene-Protein-Reaction (GPR) rule analysis** of metabolic models. It uses DIAMOND/BLAST sequence alignment to map genome genes to protein complexes and automatically generates GPR rules for SBML metabolic models.

## Features

- Upload a metabolic model (SBML/JSON) and a genome FASTA file
- DIAMOND-accelerated protein sequence search against a curated complexes database
- Automatic GPR rule generation for reactions
- Interactive resolution of ambiguous gene-to-complex assignments
- Download of updated metabolic models with new GPR rules
- Session-based job tracking (works for both authenticated and anonymous users)
- Background queue worker for processing jobs

---

## Prerequisites

### System packages

| Tool | Purpose | Install |
|------|---------|---------|
| **Python** ≥ 3.9 | Runtime | System package manager |
| **DIAMOND** | Fast protein alignment (primary) | `sudo apt-get install diamond-aligner` |
| **NCBI BLAST+** | Protein alignment (fallback) | `sudo apt-get install ncbi-blast+` |

Verify both are on your `PATH`:

```bash
diamond version
blastp -version
```

### Python packages

Install from the included requirements file:

```bash
pip install -r requirements.txt
```

---

## Installation

### 1. Clone the repository and install dependencies

```bash
git clone <repo-url>
cd GPR
pip install -r requirements.txt
```

The project structure:

```
├── manage.py
├── requirements.txt
├── run_unified_worker.py   # background job queue worker
├── config/                 # Django project settings
│   ├── settings.py
│   ├── urls.py
│   ├── wsgi.py
│   └── asgi.py
├── GPR/                    # the GPR analysis app
│   ├── admin.py
│   ├── forms.py
│   ├── models.py
│   ├── services.py
│   ├── views.py
│   ├── urls.py
│   ├── migrations/
│   ├── templates/
│   └── templatetags/
├── data/                   # reference databases
│   ├── complexes_blast_db/
│   │   ├── complexes.dmnd
│   │   └── complexes.fasta
│   └── gpr_reference/
│       ├── complex_stoichiometry.pkl
│       ├── complex_metadata.pkl
│       └── complexes_list.pkl
└── media/                  # upload & output directory
    ├── gpr_genomes/
    ├── gpr_models/
    └── sbml/
```

### 2. Run migrations

```bash
python manage.py migrate
```

### 3. Load complex metadata into the database

The `ComplexMetadata` table should be populated for the ambiguity-resolution UI to show complex names and descriptions:

```bash
python manage.py shell
```

```python
import pickle
from GPR.models import ComplexMetadata

with open('data/gpr_reference/complex_metadata.pkl', 'rb') as f:
    metadata = pickle.load(f)

for complex_id, description in metadata.items():
    ComplexMetadata.objects.update_or_create(
        complex_id=complex_id,
        defaults={'name': complex_id, 'description': description}
    )
```

### 4. Start the background queue worker

Run the worker in a separate terminal (or as a systemd service):

```bash
python run_unified_worker.py
```

The worker:
- Polls the database every 5 seconds for queued jobs
- Runs up to 4 jobs concurrently
- Writes a PID file to `media/queue_worker.pid`
- Handles SIGTERM/SIGINT gracefully

### 5. Start the development server

```bash
DJANGO_DEBUG=true python manage.py runserver
```

> **Note:** `DEBUG` defaults to `False`. You must set `DJANGO_DEBUG=true` for the
> development server to serve media files (uploaded genomes, models, results).

---

## Usage

1. Navigate to `/GPR/` in your browser
2. Enter an organism name, upload a metabolic model (`.xml`, `.sbml`, or `.json`) and a genome FASTA file (`.fasta`, `.fa`, `.fna`)
3. The job is queued and processed by the background worker
4. View results at `/GPR/job/<id>/` — BLAST hits, GPR rules, complex assignments
5. If ambiguous cases exist, resolve them at `/GPR/job/<id>/resolve-ambiguities/`
6. Download the updated model with the new GPR rules

---

## What's Included

This repository ships as a ready-to-run Django project:

| Item | Status |
|---|---|
| `manage.py` | Included |
| `config/settings.py` | Included — pre-configured with `'GPR'` in `INSTALLED_APPS` and `MEDIA_ROOT` |
| `config/urls.py` | Included — GPR routes and media serving pre-configured |
| `run_unified_worker.py` | Included — background queue worker |
| `requirements.txt` | Included |
| Static file collection | Not needed — templates use Bootstrap 5 via CDN |

### Recommended production setup

- Use a process manager (systemd, supervisor) for the queue worker instead of running it in a bare shell
- Use PostgreSQL instead of SQLite for concurrent access (`select_for_update` is used for queue locking)
- Serve `media/` via Nginx/Apache in production rather than Django's `static()` helper
- Set appropriate `FILE_UPLOAD_MAX_MEMORY_SIZE` in settings for large genome files

---

