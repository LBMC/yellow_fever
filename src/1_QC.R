setwd("~/projects/yellow_fever/")
library(scRNAtools)
devtools::load_all("../scRNAtools/", reset = T)

################################################################################
infos <- "data/Summary_SSEQ.csv"
counts = paste0("data/salmon_output/paired_end/")
id_regexp = ".*_(P[0-9]{4}_[0-9]{1,4})_.*"
id_regexp_b = "P[0-9]{4}_[0-9]{1,4}"
ERCC_regexp = "^ERCC.*"
feature = "counts"
infos_sep = ","
grouping_FUN = colSums

print("loading quant.sf files...")
dir_list <- list.dirs(counts)
dir_list <- paste0(dir_list, "/quant.sf")
dir_list <- dir_list[file.exists(dir_list)]
names(dir_list) <- gsub(id_regexp, "\\1", dir_list, perl = T)
dir_list <- dir_list[grepl(id_regexp_b, names(dir_list))]
scd_paired <- tximport(dir_list, type = "none", txOut = TRUE,
  txIdCol = "Name", abundanceCol = "TPM", countsCol = "NumReads",
  lengthCol = "EffectiveLength")
print(names(scd_paired))
Warning: 1237 parsing failures.
row # A tibble: 5 x 5 col     row      col               expected   actual expected   <int>    <chr>                  <chr>    <chr> actual 1  1261 NumReads no trailing characters  .500152 file 2  1262 NumReads no trailing characters  .499848 row 3  1329 NumReads no trailing characters    .7657 col 4  1330 NumReads no trailing characters  .141882 expected 5  1331 NumReads no trailing characters .0924225 actual # ... with 1 more variables: file <chr>
Warning: 335 parsing failures.
row # A tibble: 5 x 5 col     row      col               expected     actual expected   <int>    <chr>                  <chr>      <chr> actual 1  1521 NumReads no trailing characters    .118747 file 2  1522 NumReads no trailing characters    .881253 row 3  1570 NumReads no trailing characters    .999967 col 4  1571 NumReads no trailing characters .34432e-05 expected 5  2169 NumReads no trailing characters    .367279 actual # ... with 1 more variables: file <chr>
Warning: 735 parsing failures.
row # A tibble: 5 x 5 col     row      col               expected  actual expected   <int>    <chr>                  <chr>   <chr> actual 1  1261 NumReads no trailing characters .500519 file 2  1262 NumReads no trailing characters .499481 row 3  1609 NumReads no trailing characters .457291 col 4  1610 NumReads no trailing characters .542709 expected 5  3028 NumReads no trailing characters  .23259 actual # ... with 1 more variables: file <chr>
dir_list[1237]
dir_list[335]
dir_list[735]



for (type in c("paired_end", "single_end")) {
  for (feature in c("counts", "length", "abundance")) {
    print(paste0(type, "_", feature))
    scd <- scRNAtools::load_data_salmon(
      infos = "data/Summary_SSEQ.csv",
      counts = paste0("data/salmon_output/", type, "/"),
      feature = feature,
      id_regexp = ifelse(
        type %in% "paired_end",
        ".*_(P[0-9]{4}_[0-9]{1,4})_.*",
        ".*Run2.*_(P[0-9]{4}_[0-9]{1,4})_.*"
      ),
      tximport_obj = paste0("results/tmp/tmp_tximport_", type),
      infos_sep = ",",
      grouping_FUN = ifelse(
        feature %in% "length",
        colSums,
        max
      )
    )
    scd$setfeature(
      "to_QC",
      scd$getfeature("sex") %in% "M" & scd$getfeature("day") %in% c("D15", "D136", "D593")
    )
    save(scd, file = paste0("results/", type, "_", feature, ".Rdata"))
  }
}

load("results/single_end_abundance.Rdata")
scd_single <- scd
scd_single$setfeature("sequencing", "single")

load("results/paired_end_abundance.Rdata")
scd$setfeature("sequencing", "paired")
scd$add(
  infos = scd_single$getfeatures,
  counts = scd_single$getcounts
)
scd$setfeature(
  "to_QC",
  scd$getfeature("sex") %in% "M" & scd$getfeature("day") %in% c("D15", "D136", "D593") & scd$getfeature("sequencing") %in% "paired"
)
save(scd, file = "results/abundance.Rdata")



system("mkdir -p results/QC/QC_paraload")

scRNAtools::QC_paraload_parameters(
  paraload_file = "results/QC/paraload.csv",
  bootstraps = 100000,
  job_boot_number = 50
)

# launch paraload server
system("
bin/paraload --server \
--port 13469 \
--input results/QC/paraload.csv \
--output results/QC/paraload_run.txt \
--log results/QC/paraload.log \
--report results/QC/paraload_report.txt \
--conf src/pbs/QC/QC.conf
")

# launch paralod clients
system("
bin/paraload --client --port 13469 --host pbil-deb
")
system("
while [ $(ps -u modolo | grep paraload | wc -l) -gt 0 ]
do
stat | wc -l
iter=$(echo 200 - $(qstat -u modolo | grep -e \"[RQ]\" | wc -l) | bc)
for ((i = 1;i <= $iter;i += 1))
do
qsub src/pbs/QC/paired_end_QC.pbs &
/bin/sleep 0.5
done
/bin/sleep 3600
done
")

load("results/abundance.Rdata")
scRNAtools::QC_load_bootstraps(
  scd = scd,
  paraload_folder = "results/QC/QC_paraload/"
)
hist(scd$getfeature("QC_score"), breaks = sqrt(scd$getncells))
scRNAtools::QC_classification(scd)

save(scd, file = "results/QC/paired_end_abundance_QC.Rdata")
load("results/QC/paired_end_abundance_QC.Rdata")

check_gene(scd, "CCR7", "sex")
scd_QC <- scd$select(b_cells = scd$getfeature("QC_good") %in% T)

check_gene(scd_QC, "CCR7", "sex")

table(scd$getfeature("day"), scd$getfeature("cell_number"))

devtools::load_all("../scRNAtools/", reset = T)
b_cells = scd$getfeature('cell_number') %in% 1 & scd$getfeature("day") %in% c("D15", "D136", "D593")
system("rm results/tmp/pca_tmp.Rdata")

table(scd$getfeature("QC_good"), scd$getfeature("cell_number"))

scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "QC_good", color_name = "sex",
  tmp_file = "results/tmp/pca_tmp.Rdata")


  scRNAtools::pca_plot(
    scd$select(b_cells = b_cells), color = "day", color_name = "day",
    tmp_file = "results/tmp/pca_tmp.Rdata")

scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "sex", color_name = "sex",
  tmp_file = "results/tmp/pca_tmp.Rdata")

scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "day", color_name = "day",
  tmp_file = "results/tmp/pca_tmp.Rdata")

scRNAtools::pca_plot(
  scd$select(b_cells = b_cells)$select(), color = "cell_number", color_name = "clonality",
  tmp_file = "results/tmp/pca_tmp.Rdata")
