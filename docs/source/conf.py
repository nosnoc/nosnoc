# Configuration file for the Sphinx documentation builder.

import os

# -- Project information

project = 'nosnoc'
copyright = '2024, Armin Nurkanovic, Jonathan Frey, Anton Pozharskiy, Moritz Diehl'
author = 'Anton Pozharskiy'

release = '0.5'
version = '0.5.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinxcontrib.matlab',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'

matlab_src_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
matlab_short_links = True
matlab_auto_link = "basic"
primary_domain = "mat"
autodoc_member_order = 'bysource'
