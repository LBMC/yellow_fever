setwd("~/projects/yellow_fever/")
library(scRNAtools)
devtools::load_all("../scRNAtools/", reset = T)

system("perl -pi -e 's/[pP](\\d*_\\d*)/P\\1/g' data/Summary_SSEQ.csv")
system("perl -pi -e 's/P1306/P1316/g' data/Summary_SSEQ.csv")
################################################################################
for (feature in c("counts", "length", "abundance")) {
  print(feature)
  scd <- scRNAtools::load_data_salmon(
    infos = "data/Summary_SSEQ.csv",
    counts = "data/salmon_output/",
    feature = feature,
    tximport_obj = "results/tmp/tmp_tximport",
    infos_sep = ",",
    grouping_FUN = ifelse(
      feature %in% "length",
      colSums,
      max
    )
  )
  scd$setfeature(
    "sequencing",
    ifelse(
      grepl("P1373", scd$getfeature("id")),
      "single",
      "paired"
    )
  )
  scd$setfeature(
    "to_QC",
    scd$getfeature("sex") %in% "M" &
    scd$getfeature("day") %in% c("D15", "D136", "D593") &
    scd$getfeature("sequencing") %in% "paired"
  )
  save(scd, file = paste0("results/", feature, ".Rdata"))
  print(scd$getcells[rowSums(is.na(scd$getcounts)) != 0])
}

load("results/abundance.Rdata")
system("mkdir -p results/QC/QC_paraload/abundance/")
system("mkdir -p results/QC/QC_paraload/countximportts/")

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
--output results/QC/paraload_abundance_run.txt \
--log results/QC/paraload_abundance.log \
--report results/QC/paraload_abundance_report.txt \
--conf src/pbs/QC/QC_abundance.conf
")

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
qsub src/pbs/QC/QC_abundance.pbs &
/bin/sleep 0.5
done
/bin/sleep 3600
done
")

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

load("results/abundance.Rdata")
scRNAtools::QC_load_bootstraps(
  scd = scd,
  paraload_folder = "results/QC/QC_paraload/counts"
)
hist(scd$getfeature("QC_score"), breaks = sqrt(scd$getncells))
scRNAtools::QC_classification(scd)
table(scd$getfeature("QC_good"))
save(scd, file = "results/QC/abundance_QC.Rdata")

load("results/counts.Rdata")
scRNAtools::QC_load_bootstraps(
  scd = scd,
  paraload_folder = "results/QC/QC_paraload/counts"
)
hist(scd$getfeature("QC_score"), breaks = sqrt(scd$getncells))
scRNAtools::QC_classification(scd)
table(scd$getfeature("QC_good"))
save(scd, file = "results/QC/counts_QC.Rdata")




load("results/QC/abundance_QC.Rdata")
load("results/QC/counts_QC.Rdata")

check_gene(scd, "CCR7", "sex")
scd_QC <- scd$select(b_cells = scd$getfeature("QC_good") %in% T)

check_gene(scd_QC, "CCR7", "sex")

table(scd$getfeature("day"), scd$getfeature("cell_number"))

devtools::load_all("../scRNAtools/", reset = T)
b_cells = scd$getfeature('to_QC')
system("rm results/tmp/pca_tmp.Rdata")

table(scd$getfeature("QC_good"), scd$getfeature("day"))

devtools::load_all("../scRNAtools/", reset = T)

scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "day", color_name = "day",
  tmp_file = "results/tmp/pca_tmp.Rdata")

summary(scd$select(b_cells = b_cells)$getfeature("QC_good"))
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "QC_good", alpha = "cell_number", color_name = "antigen",
  tmp_file = "results/tmp/pca_tmp.Rdata")

scRNAtools::pca_plot(
  scd$select(b_cells = b_cells)$select(), color = "cell_number", color_name = "clonality",
  tmp_file = "results/tmp/pca_tmp.Rdata")
