text = ms_probaweb.tex
pdf = ms_probaweb.pdf

all: $(pdf)

$(pdf): $(text)
	pdflatex $(text)
