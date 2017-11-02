setwd("~/projects/yellow_fever/")
library(scRNAtools)
devtools::load_all("../scRNAtools/", reset = T)

################################################################################
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
      tximport_obj = paste0("results/tmp/tmp_tximport_", type, ".Rdata"),
      infos_sep = ",",
      grouping_FUN = ifelse(
        feature %in% "length",
        colSums,
        max
      )
    )
    save(scd, file = paste0("results/", type, "_", feature, ".Rdata"))
  }
}

load("results/paired_end_abundance.Rdata")
scd$setfeature(
  "to_QC",
  scd$getfeature("sex") %in% "M" & scd$getfeature("day") %in% c("D15", "D136", "D593")
)

system("mkdir -p results/QC/QC_paraload/paired_end")

scRNAtools::QC_paraload_parameters(
  paraload_file = "results/QC/paired_end_paraload.csv",
  bootstraps = 100000,
  job_boot_number = 50
)

# launch paraload server
system("
bin/paraload --server \
--port 13469 \
--input results/QC/paired_end_paraload.csv \
--output results/QC/paired_end_paraload_run.txt \
--log results/QC/paired_end_paraload.log \
--report results/QC/paired_end_paraload_report.txt \
--conf src/pbs/QC/paired_end_QC.conf
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

load("results/paired_end_abundance.Rdata")
scRNAtools::QC_load_bootstraps(
  scd = scd,
  paraload_folder = "results/QC/QC_paraload/paired_end/"
)
hist(scd$getfeature("QC_score"))
scRNAtools::QC_classification(scd)

hist(scd$getfeature("QC_score"), breaks = scd$getncells)


save(scd, file = "results/QC/paired_end_abundance_QC.Rdata")
load("results/QC/paired_end_abundance_QC.Rdata")

check_gene(scd, "CCR7", "sex")
scd_QC <- scd$select(b_cells = scd$getfeature("QC_good") %in% T)

check_gene(scd_QC, "CCR7", "sex")

devtools::load_all("../scRNAtools/", reset = T)
b_cells = scd$getfeature('cell_number') %in% c(1)
system("rm results/tmp/pca_tmp.Rdata")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "QC_good", color_name = "sex",
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
