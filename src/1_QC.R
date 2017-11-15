setwd("~/projects/yellow_fever/")
library(scRNAtools)
devtools::load_all("../scRNAtools/", reset = T)

system("perl -pi -e 's/[pP](\\d*_\\d*)/P\\1/g' data/Summary_SSEQ.csv")
<<<<<<< HEAD
################################################################################
Warning: 1237 parsing failures.
row # A tibble: 5 x 5 col     row      col               expected   actual expected   <int>    <chr>                  <chr>    <chr> actual 1  1261 NumReads no trailing characters  .500152 file 2  1262 NumReads no trailing characters  .499848 row 3  1329 NumReads no trailing characters    .7657 col 4  1330 NumReads no trailing characters  .141882 expected 5  1331 NumReads no trailing characters .0924225 actual # ... with 1 more variables: file <chr>
Warning: 335 parsing failures.
row # A tibble: 5 x 5 col     row      col               expected     actual expected   <int>    <chr>                  <chr>      <chr> actual 1  1521 NumReads no trailing characters    .118747 file 2  1522 NumReads no trailing characters    .881253 row 3  1570 NumReads no trailing characters    .999967 col 4  1571 NumReads no trailing characters .34432e-05 expected 5  2169 NumReads no trailing characters    .367279 actual # ... with 1 more variables: file <chr>
Warning: 735 parsing failures.
row # A tibble: 5 x 5 col     row      col               expected  actual expected   <int>    <chr>                  <chr>   <chr> actual 1  1261 NumReads no trailing characters .500519 file 2  1262 NumReads no trailing characters .499481 row 3  1609 NumReads no trailing characters .457291 col 4  1610 NumReads no trailing characters .542709 expected 5  3028 NumReads no trailing characters  .23259 actual # ... with 1 more variables: file <chr>

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
=======
system("perl -pi -e 's/P1306/P1316/g' data/Summary_SSEQ.csv")
################################################################################
for (feature in c("counts", "length", "abundance")) {
  print(feature)
  scd <- scRNAtools::load_data_salmon(
    infos = "data/Summary_SSEQ.csv",
    counts = "data/salmon_output/",
    feature = feature,
    tximport_obj = "results/tmp/tmp_tximport",
    infos_sep = ","
  )
  scd$setfeature(
    "sequencing",
    ifelse(
      grepl("P1373", scd$getfeature("id")),
      "single",
      "paired"
>>>>>>> refs/remotes/origin/master
    )
  )
  scd$setfeature(
    "to_QC",
    scd$getfeature("sex") %in% "M" &
    scd$getfeature("day") %in% c("D15", "D136", "D593") &
    scd$getfeature("sequencing") %in% "paired" &
    !(scd$getfeature("batch") %in% c(6:8)) &
    !(scd$getcells %in% c("P1299_1797", "P1299_1896"))
  )
  save(scd, file = paste0("results/", feature, ".Rdata"))
  print(scd$getcells[rowSums(is.na(scd$getcounts)) != 0])
}

load("results/abundance.Rdata")
b_cells <- scd$getfeature("to_QC")
  | scd$getcells %in% c("P1299_1797", "P1299_1896")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "batch", color_name = "clonality",
  tmp_file = "results/tmp/pca_abundance_outliers_tmp.Rdata")
load("results/tmp/pca_abundance_outliers_tmp.Rdata")
scd$select(b_cells = b_cells)$getfeatures[order(pca_data$x)[1:2],]

load("results/counts.Rdata")
b_cells <- scd$getfeature("to_QC")
  | scd$getcells %in% c("P1299_1797", "P1299_1896")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "batch", color_name = "clonality",
  tmp_file = "results/tmp/pca_counts_outliers_tmp.Rdata")
load("results/tmp/pca_counts_tmp.Rdata", v = T)
scd$select(b_cells = b_cells)$getfeatures[order(pca_data$y, decreasing = T)[1:2],]

################################################################################

load("results/counts.Rdata")
system("mkdir -p results/QC/QC_paraload/abundance/")
system("mkdir -p results/QC/QC_paraload/counts/")

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
--output results/QC/paraload_counts_run.txt \
--log results/QC/paraload_counts.log \
--report results/QC/paraload_counts_report.txt \
--conf src/pbs/QC/QC_counts.conf
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
qsub src/pbs/QC/QC_counts.pbs &
/bin/sleep 0.5
done
/bin/sleep 3600
done
")

load("results/counts.Rdata")
scRNAtools::QC_load_bootstraps(
  scd = scd,
  paraload_folder = "results/QC/QC_paraload/counts"
)
hist(scd$getfeature("QC_score"), breaks = sqrt(scd$getncells))
scRNAtools::QC_classification(scd)
table(scd$getfeature("QC_good"))
save(scd, file = "results/QC/counts_QC.Rdata")

counts_features <- scd$getfeatures

load("results/abundance.Rdata")
scd <- scdata$new(
  infos = counts_features,
  counts = scd$getcounts,
  v = T
)
save(scd, file = "results/QC/abundance_QC.Rdata")


load("results/QC/counts_QC.Rdata")

table(scd$getfeature("day"), scd$getfeature("cell_number"))
table(scd$getfeature("day"), scd$getfeature("QC_good"))
b_cells = scd$getfeature('to_QC')

system("rm results/tmp/pca_*_QC_tmp.Rdata")

load("results/QC/counts_QC.Rdata")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "day", color_name = "day", alpha = "QC_good",
  tmp_file = "results/tmp/pca_counts_QC_tmp.Rdata")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "batch", color_name = "clonality", alpha = "QC_good",
  tmp_file = "results/tmp/pca_counts_QC_tmp.Rdata")


load("results/QC/abundance_QC.Rdata")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "day", color_name = "day", alpha = "QC_good",
  tmp_file = "results/tmp/pca_abundance_QC_tmp.Rdata")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "batch", color_name = "clonality", alpha = "QC_good",
  tmp_file = "results/tmp/pca_abundance_QC_tmp.Rdata")


devtools::load_all("../scRNAtools/", reset = T)
load("results/QC/counts_QC.Rdata")

scd <- normalize(
  scd = scd,
  method = "SCnorm",
  cpus = 4,
  tmp_file = "results/tmp/normalization_tmp.Rdata",
)

save(scd, file = "results/QC/norm_counts_QC.Rdata")
system("rm results/tmp/pca_norm_counts_QC_tmp.Rdata")
b_cells <- scd$getfeature("QC_good") %in% T
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "day", color_name = "day",
  tmp_file = "results/tmp/pca_norm_counts_QC_tmp.Rdata")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "batch", color_name = "clonality",
  tmp_file = "results/tmp/pca_norm_counts_QC_tmp.Rdata")

