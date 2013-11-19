text = ms_probaweb.md
pdf = probaweb_PDF.pdf

all: $(pdf)

$(pdf): $(text)
	pandoc $(text) -o $(pdf)
