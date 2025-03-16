# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'ドッキング自動化システム'
copyright = '2025, 開発チーム'
author = '開発チーム'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'sphinx.ext.viewcode',
    'sphinx_autodoc_typehints',
    'sphinx.ext.graphviz',
    'sphinx.ext.inheritance_diagram',
]

templates_path = ['_templates']
exclude_patterns = []

language = 'ja'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# -- autodoc の設定 ---------------------------------------------------------
autodoc_member_order = 'bysource'
autoclass_content = 'both'
autodoc_typehints = 'description'

# -- napoleon の設定 --------------------------------------------------------
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = True
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
# LaTeX関連の設定を更新 (XeTeX用)
latex_engine = 'xelatex'
latex_elements = {
    'papersize': 'a4paper',
    'pointsize': '10pt',
    'preamble': r'''
\usepackage{xeCJK}
\setCJKmainfont{IPAMincho}
\setCJKsansfont{IPAGothic}
\setCJKmonofont{IPAGothic}
\XeTeXlinebreaklocale "ja"
\XeTeXlinebreakskip = 0pt plus 1pt
''',
    'figure_align': 'htbp',
    'fncychap': r'\usepackage[Bjarne]{fncychap}',
    'extraclassoptions': ',openany,oneside',
    'maketitle': r'\newcommand\sphinxbackoftitlepage{}\sphinxmaketitle',
    'sphinxsetup': r'warningBgColor={rgb}{1,1,1}',
}
