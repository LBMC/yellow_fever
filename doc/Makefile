all: 1_QC.pdf 2_cell_type.pdf

1_QC.pdf: 1_QC.md
	R -e 'require("rmarkdown"); rmarkdown::render("1_QC.md")'
2_cell_type.pdf: 2_cell_type.md
	R -e 'require("rmarkdown"); rmarkdown::render("2_cell_type.md")'
