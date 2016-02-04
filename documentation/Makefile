PHYS=physics_doc


physics:
	pdflatex $(PHYS) ; bibtex $(PHYS) ; pdflatex $(PHYS) ; pdflatex $(PHYS)

all:	physics

clean:
	rm -rf *dvi *idx *log *toc *bbl *aux *blg *out *~ *spl  $(PHYS)Notes.bib

deepclean: clean
	rm $(PHYS).pdf
