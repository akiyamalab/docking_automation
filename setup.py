# type: ignore
from setuptools import find_packages, setup

setup(
    name="docking_automation",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "dask",
        "dask-jobqueue",
        "matplotlib",
        "scipy",
        "sortedcontainers",
        "openbabel-wheel>=3.0.0",
        "numpy<2.0.0",  # rdkit の都合でnumpyは2.0.0未満
        "rdkit==2023.9.1",  # meekoの関係で2023.9.6より新しいものは使えない、2023.9.6はrdkit-stubsのバグで使えない
        "rdkit-stubs",
        "meeko==0.5.0",
        "tqdm", 
        "types-tqdm",
        # 以下は alphacutter の依存関係
        "biopython",
        "biopandas",
        "networkx",
        "pcmap",
    ],
    extras_require={
        "dev": [
            "pytest",
            "pytest-cov",
            "mypy",
            "flake8",
            "black",
            "sphinx",
            "types-setuptools",  # setuptoolsの型定義
        ],
    },
    python_requires=">=3.10",
    author="Your Name",
    author_email="your.email@example.com",
    description="A domain-driven design approach to protein-compound docking automation",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
    ],
)
