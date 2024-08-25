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
        "regions==0.8"
        ]

)
