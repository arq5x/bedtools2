# Makefile template for Sphinx documentation

# You can set these variables from the command line.
SPHINXOPTS      =
PAPER           =
SPHINX_BUILDDIR = _build

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) $(top_srcdir)/docs

.PHONY: default all help clean html dirhtml singlehtml pickle json htmlhelp qthelp devhelp epub latex latexpdf text man changes linkcheck doctest

if ENABLE_DOCUMENTATION
all: $(pdf_DATA) .html-stamp

# pdf target
pdf_DATA = $(SPHINX_BUILDDIR)/latex/bedtools.pdf
$(SPHINX_BUILDDIR)/latex/bedtools.pdf:
	$(MAKE) latexpdf

# html target
.html-stamp:
	$(MAKE) html

install-data-local:
	$(MKDIR_P) "$(DESTDIR)$(htmldir)" || exit 1
	@cp -R $(SPHINX_BUILDDIR)/html/. "$(DESTDIR)$(htmldir)" || exit 1
endif

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make standalone HTML files"
	@echo "  dirhtml    to make HTML files named index.html in directories"
	@echo "  singlehtml to make a single large HTML file"
	@echo "  pickle     to make pickle files"
	@echo "  json       to make JSON files"
	@echo "  htmlhelp   to make HTML files and a HTML help project"
	@echo "  qthelp     to make HTML files and a qthelp project"
	@echo "  devhelp    to make HTML files and a Devhelp project"
	@echo "  epub       to make an epub"
	@echo "  latex      to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  latexpdf   to make LaTeX files and run them through pdflatex"
	@echo "  text       to make text files"
	@echo "  man        to make manual pages"
	@echo "  changes    to make an overview of all changed/added/deprecated items"
	@echo "  linkcheck  to check all external links for integrity"
	@echo "  doctest    to run all doctests embedded in the documentation (if enabled)"

clean:
	-rm -rf $(SPHINX_BUILDDIR)/*
	-rm -f .html-stamp .pdf-stamp

html:
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(SPHINX_BUILDDIR)/html && \
		touch .html-stamp
	@echo
	@echo "Build finished. The HTML pages are in $(SPHINX_BUILDDIR)/html."

dirhtml:
	$(SPHINXBUILD) -b dirhtml $(ALLSPHINXOPTS) $(SPHINX_BUILDDIR)/dirhtml
	@echo
	@echo "Build finished. The HTML pages are in $(SPHINX_BUILDDIR)/dirhtml."

singlehtml:
	$(SPHINXBUILD) -b singlehtml $(ALLSPHINXOPTS) $(SPHINX_BUILDDIR)/singlehtml
	@echo
	@echo "Build finished. The HTML page is in $(SPHINX_BUILDDIR)/singlehtml."

pickle:
	$(SPHINXBUILD) -b pickle $(ALLSPHINXOPTS) $(SPHINX_BUILDDIR)/pickle
	@echo
	@echo "Build finished; now you can process the pickle files."

json:
	$(SPHINXBUILD) -b json $(ALLSPHINXOPTS) $(SPHINX_BUILDDIR)/json
	@echo
	@echo "Build finished; now you can process the JSON files."

htmlhelp:
	$(SPHINXBUILD) -b htmlhelp $(ALLSPHINXOPTS) $(SPHINX_BUILDDIR)/htmlhelp
	@echo
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in $(SPHINX_BUILDDIR)/htmlhelp."

qthelp:
	$(SPHINXBUILD) -b qthelp $(ALLSPHINXOPTS) $(SPHINX_BUILDDIR)/qthelp
	@echo
	@echo "Build finished; now you can run "qcollectiongenerator" with the" \
	      ".qhcp project file in $(SPHINX_BUILDDIR)/qthelp, like this:"
	@echo "# qcollectiongenerator $(SPHINX_BUILDDIR)/qthelp/track.qhcp"
	@echo "To view the help file:"
	@echo "# assistant -collectionFile $(SPHINX_BUILDDIR)/qthelp/track.qhc"

devhelp:
	$(SPHINXBUILD) -b devhelp $(ALLSPHINXOPTS) $(SPHINX_BUILDDIR)/devhelp
	@echo
	@echo "Build finished."
	@echo "To view the help file:"
	@echo "# mkdir -p $$HOME/.local/share/devhelp/track"
	@echo "# ln -s $(SPHINX_BUILDDIR)/devhelp $$HOME/.local/share/devhelp/track"
	@echo "# devhelp"

epub:
	$(SPHINXBUILD) -b epub $(ALLSPHINXOPTS) $(SPHINX_BUILDDIR)/epub
	@echo
	@echo "Build finished. The epub file is in $(SPHINX_BUILDDIR)/epub."

latex:
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(SPHINX_BUILDDIR)/latex
	@echo
	@echo "Build finished; the LaTeX files are in $(SPHINX_BUILDDIR)/latex."
	@echo "Run \`make' in that directory to run these through (pdf)latex" \
	      "(use \`make latexpdf' here to do that automatically)."

latexpdf:
	@echo "Running LaTeX files through pdflatex..."
	$(MAKE) latex
	$(MAKE) -C $(SPHINX_BUILDDIR)/latex all-pdf
	@echo "pdflatex finished; the PDF files are in $(SPHINX_BUILDDIR)/latex."

text:
	$(SPHINXBUILD) -b text $(ALLSPHINXOPTS) $(SPHINX_BUILDDIR)/text
	@echo
	@echo "Build finished. The text files are in $(SPHINX_BUILDDIR)/text."

man:
	$(SPHINXBUILD) -b man $(ALLSPHINXOPTS) $(SPHINX_BUILDDIR)/man
	@echo
	@echo "Build finished. The manual pages are in $(SPHINX_BUILDDIR)/man."

changes:
	$(SPHINXBUILD) -b changes $(ALLSPHINXOPTS) $(SPHINX_BUILDDIR)/changes
	@echo
	@echo "The overview file is in $(SPHINX_BUILDDIR)/changes."

linkcheck:
	$(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) $(SPHINX_BUILDDIR)/linkcheck
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in $(SPHINX_BUILDDIR)/linkcheck/output.txt."

doctest:
	$(SPHINXBUILD) -b doctest $(ALLSPHINXOPTS) $(SPHINX_BUILDDIR)/doctest
	@echo "Testing of doctests in the sources finished, look at the " \
	      "results in $(SPHINX_BUILDDIR)/doctest/output.txt."
