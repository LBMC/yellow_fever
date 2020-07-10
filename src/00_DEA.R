#!Rscript
source("src/00_functions.R")
# /src/00_DEA.R --counts_matrix results/test_counts_matrix.csv --cells_matrix results/test_cells_matrix.csv --output results/test_DEA.csv --formula "p_PLS_DEA_cell_type + clonality" --test "p_PLS_DEA_cell_type" --cpus 6
option_list = list(
  make_option(
    c("-x", "--counts_matrix"),
    type = "character",
    default = NULL,
    help = "csv files of the counts (genes are row and cells as columns)",
    metavar = "character"
  ),
  make_option(
    c("-c", "--cells_matrix"),
    type = "character",
    default = NULL,
    help = "csv files of the cells (cells are row and cell infos as columns)",
    metavar = "character"
  ),
  make_option(
    c("-f", "--formula"),
    type = "character",
    default = NULL,
    help = "formula to model the --cells_matrix columns (e.g: \"DEA_cell_type + clonality\")",
    metavar = "character"
  ),
  make_option(
    c("-t", "--test"),
    type = "character",
    default = NULL,
    help = "name of the variable to test for in the formula (e.g: \"DEA_cell_type\")",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = NULL,
    help = "name of csv file to write for the results",
    metavar = "character"
  ),
  make_option(
    c("-m", "--cpus"),
    type = "character",
    default = NULL,
    help = "number of cpus to use for the computation",
    metavar = "character"
  )
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

print(opt)

sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(
    counts = data.table::fread(opt$counts_matrix) %>%
      as.matrix() %>%
      Matrix::Matrix(sparse = T)
  ),
  colData = read.csv(opt$cells_matrix),
  rowData = data.frame(
    id = data.table::fread(opt$counts_matrix) %>% rownames()
  )
)
colnames(sce) <- colData(sce) %>% rownames()
rownames(sce) <- rowData(sce)$id


DEA_results <- DEA(
  sce = sce,
  test = opt$test,
  formula = str_c("count ~ 1 + ", opt$formula),
  assay_name = "counts",
  cpus = as.numeric(opt$cpus)
)


rowData(sce)$pval_DEA <- NA
rowData(sce)$pval_DEA <- 
  get_genes_pval(DEA_results, sce)

rowData(sce)$pval_DEA_adj <- p.adjust(
  rowData(sce)$pval_DEA,
  method = "BH"
)

rowData(sce) %>% 
  as.data.frame() %>% 
  write_csv(, path = opt$output)