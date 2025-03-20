# type: ignore
from setuptools import setup, find_packages

setup(
    name="docking_automation",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.22.0",       # 数値計算ライブラリ
        "pandas>=1.4.0",       # データ解析ライブラリ
        "sqlalchemy>=1.4.0",   # ORM/データベース操作
        "pydantic>=1.9.0",     # データバリデーション
        "typer>=0.6.0",        # CLIツール
        "dask>=2022.5.0",      # 並列計算
        "dask-jobqueue>=0.8.0", # ジョブスケジューラ連携
        "psutil>=5.9.0",       # システムリソース情報
        "matplotlib>=3.5.0",   # グラフ描画
        "scipy>=1.8.0",        # 科学計算
        "sortedcontainers>=2.4.0",  # ソート済みコンテナ
        # 以下のパッケージは以前はオプショナルでしたが、常にインストールされるようになりました
        "rdkit>=2022.3.1",        # 化学情報学ライブラリ（conda経由でインストール推奨）
        "openbabel-wheel>=3.1.0",  # 分子ファイル形式変換
        "meeko>=0.4.0",         # リガンド準備（AutoDock Vina用）
    ],
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.1.0",
            "mypy>=0.960",
            "flake8>=4.0.0",
            "black>=22.3.0",
            "sphinx>=5.0.0",
            "types-setuptools",  # setuptoolsの型定義
        ],
    },
    python_requires=">=3.8",
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
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
)