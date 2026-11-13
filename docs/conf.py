# Sphinx configuration for the GUST documentation.
#
#     make -C docs html        (warnings are errors; see docs/Makefile)

project = "GUST"
author = "GUST contributors"
copyright = "2026, GUST contributors"

extensions = ["sphinx.ext.mathjax"]
templates_path = []
exclude_patterns = ["_build"]

html_theme = "sphinx_rtd_theme"
html_static_path = []

master_doc = "index"

# Labelled equations are numbered and can be cited with :eq:`label`.
numfig = True

# PDF (make -C docs latexpdf): xelatex, because the text carries Unicode
# (primes, Greek letters, dashes) that pdflatex refuses. Sphinx's default
# xelatex fonts (GNU FreeFont, by file name) come with TeX Live and, on
# Ubuntu, with fonts-freefont-otf.
latex_engine = "xelatex"
latex_documents = [(master_doc, "gust.tex", "GUST documentation", author, "manual")]
latex_elements = {
    "papersize": "letterpaper",
    "pointsize": "10pt",
}
