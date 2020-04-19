library(tidyverse)

data_dir <- "data/salmon_output/"

# build the transcript to gene table
tx2g <- list.files(
  path = data_dir,
  pattern = ".*quant\\.sf",
  recursive = T) %>%
  .[1] %>%
  stringr::str_c(data_dir, .) %>%
  readr::read_tsv() %>%
  dplyr::select(Name) %>%
  dplyr::rename(Transcript = Name) %>%
  dplyr::mutate(Gene = Transcript %>% sapply(
    function(x) {
      ifelse(
        stringr::str_detect(x, "\\|"),
        stringr::str_split(x, "\\|", simplify = T)[1, 2],
        x
      )
    }
  ))

# import counts
counts <- list.files(
  path = data_dir,
  pattern = ".*quant\\.sf",
  recursive = T) %>%
  stringr::str_c(data_dir, .) %>%
  tximport::tximport(
    files = .,
    type = "salmon",
    tx2gene = tx2g,
    importer = data.table::fread
    )

# read cell id from file paths
cell_id <- list.files(
  path = data_dir,
  pattern = ".*quant\\.sf",
  recursive = T) %>%
  stringr::str_replace(
    pattern = ".*[_/]P(\\d+)_(\\w{0,1}\\d+)[_/].*$",
    replacement = "P\\1_\\2"
  ) %>%
  stringr::str_replace(
    pattern = "output/(.).*(\\d+)_S.*$",
    replacement = "P999\\1_\\2"
  )

# build singleCellExperiment object
sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(
    counts_raw = counts$counts %>% Matrix::Matrix(sparse = T),
    abundance_raw = counts$abundance %>% Matrix::Matrix(sparse = T)
  ),
  colData = data.frame(
    id = cell_id
  ),
  rowData = data.frame(
    id = tx2g$Gene[match(rownames(counts$counts), tx2g$Gene)],
    infos = tx2g$Transcript[match(rownames(counts$counts), tx2g$Gene)]
  )
)
save(sce, file = "results/sce_raw.Rdata")

# add cells annotation from annotation file
library(SingleCellExperiment)
library(SummarizedExperiment)
load("results/sce_raw.Rdata", verbose = T)

colData(sce) <- colData(sce) %>%
  as_tibble() %>%
  mutate(id = as.character(id)) %>% 
  dplyr::left_join(read_csv("data/Summary_SSEQ.csv")) %>% 
  dplyr::left_join(read_csv("data/2020_01_22_P1856_Meta.csv")) %>%
  mutate(
    experiment = stringr::str_replace(id, "(P\\d+)_.*", "\\1"),
    sequencing = "paired",
    sequencing = ifelse(stringr::str_detect(experiment, "P1373"),
                        "single",
                        sequencing),
    to_QC = F
  ) %>% 
  DataFrame()

# add gene name info
rowData(sce) <- rowData(sce) %>% 
  as_tibble() %>% 
  dplyr::mutate(gene = infos %>% sapply(
    function(x) {
      ifelse(
        stringr::str_detect(x, "\\|"),
        stringr::str_split(x, "\\|", simplify = T)[1, 6],
        x
      )
    }
  )) %>% 
  DataFrame

save(sce, file = "results/sce_annot.Rdata")
