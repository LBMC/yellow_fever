all: mat_meth.pdf mat_meth.docx

mat_meth.pdf: mat_meth.tex
	xelatex mat_meth.tex

mat_meth.docx: mat_meth.tex bibliography.bib
	pandoc mat_meth.tex -o mat_meth.docx  --pdf-engine=xelatex --bibliography bibliography.bib

