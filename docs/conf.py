# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
from docutils import nodes, utils

# sys.path.insert(0, os.path.abspath('.'))

# sys.path.insert(0, os.path.abspath('..'))  # Adjust this as necessary
#
# sys.path.insert(0, os.path.abspath('../pyranges'))  # Adjust this as necessary
#
# sys.path.insert(0, os.path.abspath('.'))  # Adjust this as necessary

# -- Project information -----------------------------------------------------

project = "pyranges"
copyright = "2024, Endre Bakken Stovner, Marco Mariotti"
author = "Endre Bakken Stovner, Marco Mariotti"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.


extensions = [
    "sphinx.ext.napoleon",
    # "sphinxcontrib.napoleon",
    # "autoapi.extension",
    "sphinx.ext.autodoc",
    # "autoapi.extension",
    "sphinx.ext.autosummary",
    "sphinx.ext.doctest",
    # "sphinx.ext.doctest"
]

autosummary_generate = True  # Enable summary table generation


autodoc_default_options = {
    "members": True,
    "imported-members": True,
    # other options...
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_generated_hidden*"]

master_doc = "index"

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = "alabaster"
html_theme = "sphinx_rtd_theme"


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]


def monospaced_link(name, rawtext, text, lineno, inliner, options={}, content=[]):
    url = text.split(" ")[-1].strip("<>")
    clickable_text = " ".join(text.split(" ")[:-1])
    # Create a reference node, which is the docutils node for hyperlinks
    unescaped_text = utils.unescape(text)

    node = nodes.reference(rawtext, clickable_text, refuri=url, **options)

    # Add a special class to this node
    node["classes"].append("monospaced-link")
    return [node], []


def setup(app):
    app.add_role("mslink", monospaced_link)
    app.add_css_file("custom.css")
