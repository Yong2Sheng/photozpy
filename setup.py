from setuptools import setup, find_namespace_packages, find_packages

setup(
    name = "photozpy",
    python_requires=">=3.7",
    version = "0.1",
    author = "Yong Sheng",
    author_email = "sheng2@clemson.edu",
    description = "Automatic pipeline for data analysis UVOT and SARA images.",
    
    package_dir = {"": "src"},

    packages = find_packages(where="src"),

    install_requires = [
        "numpy",
        "pandas",
        "tqdm",
        "astropy",
        "swifttools",
        "ipywidgets",
        "jupyterlab",
        "chardet",
        "ccdproc",
        "photutils",
        "astroalign",
        "scipy",
        "matplotlib",
        "astroquery",
        "regions"
        ]

)
