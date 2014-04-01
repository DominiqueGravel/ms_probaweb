text = ms_probaweb
pdf = ms_probaweb.pdf

all: $(pdf)

$(pdf): $(text).tex
	pdflatex $(text).tex
	bibtex $(text)
	pdflatex $(text).tex
