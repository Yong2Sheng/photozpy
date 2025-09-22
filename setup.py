from setuptools import setup, find_packages  # find_namespace_packages

setup(
    name="photozpy",
    python_requires="==3.9.19",
    version="0.1.0",
    author="Yong Sheng",
    author_email="sheng2@clemson.edu",
    description="Automatic pipeline for data analysis UVOT and SARA images.",
    # Use the recommended src/ layout.
    package_dir={"": "src"},
    packages=find_packages(where="src"),  # finds src/photozpy and subpackages
    # --- Runtime dependencies (what your code imports at run time) ---
    install_requires=[
        "numpy==1.26.4",
        "pandas==2.2.2",
        "tqdm",
        "astropy==6.0.1",
        "swifttools==3.0.21",
        "ipywidgets",
        "jupyterlab",
        "chardet",
        "ccdproc==2.4.2",
        "photutils==1.11.0",
        "astroalign==2.5.1",
        "scipy==1.13.1",
        "matplotlib==3.9.0",
        "astroquery==0.4.7",
        "regions==0.8",
    ],
    # --- Dev extras (tools needed for testing, linting, typing; not runtime) ---
    # they are not needed by the users at runtime, this makes the module lightweighted.
    # to install with the extra require, run pip install -e ".[dev]"
    extras_require={
        "dev": [
            "pytest>=8",
            "pytest-cov>=5",
            "coverage>=7",
            "ruff>=0.5",
            "black>=24.0",
            "mypy>=1.8",
            "pre-commit>=3.6",
        ]
    },
)
