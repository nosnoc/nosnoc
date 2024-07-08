# Configuration file for the Sphinx documentation builder.

import os

# -- Project information

project = 'nosnoc'
copyright = '2024, Armin Nurkanovic, Jonathan Frey, Anton Pozharskiy, Moritz Diehl'
author = 'Armin Nurkanovic'

release = '0.5'
version = '0.5.0'

# -- General configuration

extensions = [
    'sphinx.ext.napoleon',
    'sphinx.ext.autodoc',
    'sphinx.ext.todo',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinxcontrib.matlab',
    'sphinx.ext.mathjax',
    'sphinx_math_dollar',
    'sphinxcontrib.bibtex',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
    'vdx': ('https://vdx.readthedocs.io/en/latest/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'furo'
pygments_style = "sphinx"
pygments_dark_style = "monokai"

# -- Options for EPUB output
epub_show_urls = 'footnote'

# -- Options for bibtex
bibtex_bibfiles = ['refs.bib']
bibtex_reference_style = 'label'
bibtex_default_style = 'plain'

# -- Other Options TODO(@anton) organize and document.

matlab_src_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
matlab_keep_package_prefix = False
matlab_show_property_default_value = True
primary_domain = "mat"
default_role = "obj"
autodoc_member_order = 'bysource'

autodoc_default_options = {'members': True, 'show-inheritance': True}
autosummary_generate = True

# -- magic to prevent mathjax from overriding sphinx_math_dollar.

mathjax3_config = {
    "tex": {
        "inlineMath": [['\\(', '\\)']],
        "displayMath": [["\\[", "\\]"]],
        "macros": {
            # Math notation
            "R": "\\mathbb{R}",
            "C": "\\mathcal{C}",
            "xdot": "\\dot{x}",
            "Set": ["\\left\\{\\, #1\\ \\middle|\\ #2\\,\\right\\}", 2],
            "argmin": "\\operatorname{arg\\,min}",
            "argmax": "\\operatorname{arg\\,max}",
            "conv": "\\operatorname{conv}",
        }
    }
}
