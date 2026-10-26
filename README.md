# photozpy

A Python pipeline for photometry and photo-z–related analysis (Swift/UVOT + SARA images).

## Recommended workflow

- Use **conda** to manage Python dependencies.
- Install the external **Astrometry.net** solver separately for local plate solving (see step 5).
- Use **pip only** for installing `photozpy` in **editable** mode:
  - `python -m pip install -e . --no-deps`

## Installation

### 1) Install conda

Install Anaconda (or Miniconda/Miniforge3) following the installation guide for your chosen distribution.

### 2) Clone the repository

```bash
git clone https://github.com/Yong2Sheng/photozpy.git
cd photozpy
```

### 3) Create and activate the conda environment

Create the environment from `environment.yml`. The `-n photozpy` option overrides the name defined in the YAML file:

```bash
conda env create -f environment.yml -n photozpy
```

Activate the environment using the name specified above:

```bash
conda activate photozpy
```

If you edit `environment.yml` later, update the environment with:

```bash
conda env update -f environment.yml -n photozpy --prune
```

### 4) Install `photozpy` in editable mode and install pre-push hooks

From the repository root (same folder as `pyproject.toml`):

```bash
python -m pip install -e . --no-deps
pre-commit install -f

```

Manually run the same checks that would run on push:

```bash
pre-commit run --all-files --hook-stage pre-push
```

### 5) Install and configure local plate solving on macOS

Local plate solving requires the external Astrometry.net program `solve-field` and astrometric index files. These are installed separately from the conda environment. They are not required when using the online Astrometry.net API.

#### Install Astrometry.net

Install [Homebrew](https://brew.sh/) first, then run:

```bash
brew install astrometry-net
```

On Apple Silicon, use native ARM Homebrew. After activating the Photozpy conda environment, check that the solver is available:

```bash
command -v solve-field
solve-field --version
```

Photozpy invokes `solve-field` through `PATH`, so it must also be available to the Python or Jupyter process running the pipeline.

#### Download the SARA index files

For the SARA setup, download the 5200-series LIGHT indexes `5203`, `5204`, and `5205`, with 48 tiles per scale (144 files total). Other telescope setups may require different indexes; see the [Astrometry.net index guide](https://astrometry.net/doc/readme.html#getting-index-files).

Install `wget` if it is not already available:

```bash
brew install wget
```

Run the following from the repository root in Bash or zsh. The files are downloaded from the [5200-series LIGHT directory](https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/) into `data/astrometry-indexes/` inside the repository:

```bash
mkdir -p data/astrometry-indexes

for j in 5203 5204 5205; do
  for ((i=0; i<48; i++)); do
    I=$(printf '%02d' "$i")
    wget -c -P data/astrometry-indexes \
      "https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/index-${j}-${I}.fits" || break 2
  done
done
```

If a download fails, resolve the reported error and rerun the loop; `-c` resumes partial downloads.

Only `/data/astrometry-indexes/` is excluded by `.gitignore`; other small files under `data/` can still be tracked in Git. The indexes are local external data and are not included when cloning the repository.

#### Configure the index directory

After `brew install astrometry-net` completes successfully, locate the default configuration file provided by the installation:

```bash
ls "$(brew --prefix astrometry-net)/etc/astrometry.cfg"
```

From the repository root, obtain the absolute index directory path:

```bash
(cd data/astrometry-indexes && pwd)
```

Open `astrometry.cfg` in a text editor. Set an `add_path` entry to the absolute directory path printed above, replacing the default index-directory entry or adding one if none is present. Enable automatic index loading by ensuring the file contains a line consisting of `autoindex` without a leading `#`. If the line is `#autoindex`, remove the `#`; if it is absent, add it. The relevant settings should look like this:

```text
add_path /absolute/path/to/photozpy/data/astrometry-indexes
autoindex
```

Replace `/absolute/path/to/photozpy/data/astrometry-indexes` with the actual absolute directory path; do not use `~` or `$HOME` inside the configuration file. Keep the other configuration settings unchanged and save the file. `autoindex` is a configuration directive, not a separate package to install.

#### Verify local solving

Test with a FITS image and write the outputs to a new directory:

```bash
solve-field /absolute/path/to/test.fits \
  --no-plots \
  --dir /absolute/path/to/new-test-output
```

Confirm that solving succeeds and inspect the resulting WCS. The version check alone does not verify that the index files are configured correctly.

### 6) Install DS9 for interactive workflows

Photozpy's interactive DS9 workflows require a separately installed [SAOImage DS9](https://sites.google.com/cfa.harvard.edu/saoimageds9/download). Install the version appropriate for your operating system and open DS9 before starting an interactive workflow; Photozpy does not launch DS9 automatically.

Photozpy communicates with DS9 through the XPA commands `xpaget` and `xpaset`. The `xpa` package is included in `environment.yml`; both commands must be available on `PATH` in the Python or Jupyter environment running Photozpy. With that environment activated and DS9 running, verify communication with:

```bash
xpaget ds9 version
```

This should return the running DS9 version. DS9 is required for the interactive DS9 features, rather than for every pipeline operation.

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
