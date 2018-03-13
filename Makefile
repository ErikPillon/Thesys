title = thesys

chapters = ch/*.tex
bibliography = bibliography.bib

all: $(title).tex

clean:
	rm -f *.pdf *.toc *.aux *.out *.log *.dvi *.bbl

$(title).pdf: $(title).tex
	pdflatex $(title).tex

$(title).tex: $(chapters) $(bibliography)
	touch $(title).tex


full:
	pdflatex $(title).tex
	bibtex $(title).aux
	pdflatex $(title).tex
	pdflatex $(title).tex
	




foo :
	for i in $(chapters); do \
		echo $$i ; \
		done


