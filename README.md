# photozpy

A Python pipeline for photometry and photo-z–related analysis (Swift/UVOT + SARA images).

## Recommended workflow

- Use **conda** to manage *all dependencies*
- Use **pip only** for installing `photozpy` in **editable** mode:
  - `python -m pip install -e . --no-deps`

## Installation

### 1) Install conda

Install Anaconda (or Miniconda) following the official guide:
https://docs.anaconda.com/free/anaconda/install/index.html

(Optional) Use conda-forge with strict channel priority to reduce dependency conflicts:

```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
```

### 2) Clone the repository

```bash
git clone https://github.com/Yong2Sheng/photozpy.git
cd photozpy
```

### 3) Create and activate the conda environment

Create the environment from `environment.yml`:

```bash
conda env create -f environment.yml
```

Activate the environment (the env name is defined by `name:` in `environment.yml`):

```bash
conda activate <ENV_NAME>
```

If you edit `environment.yml` later, update the environment with:

```bash
conda env update -f environment.yml --prune
```

### 4) Install `photozpy` in editable mode and install pre-push hooks

From the repository root (same folder as `pyproject.toml`):

```bash
python -m pip install -e . --no-deps
pre-commit install -f

```

## Quick start

Tutorial notebooks are under `docs/`, for example:

- `docs/Tutorial.ipynb`

## Testing

Smoke test (used by the pre-push hook):

```bash
pytest -m smoke
```

Full test + coverage (HTML report will be written to `tests/coverage_report/`):

```bash
pytest --cov=photozpy --cov-report=term --cov-report=html:tests/coverage_report
```

## License

(TODO) Add license information here.
