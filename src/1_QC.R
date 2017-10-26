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
        feature %in% "paired_end",
        ".*_(P[0-9]{4}_[0-9]{1,4})_.*",
        ".*Run2.*_(P[0-9]{4}_[0-9]{1,4})_.*"
      ),
      tximport_obj = paste0("results/tmp/tmp_tximport_", type, ".Rdata"),
      infos_sep = ",",
      grouping_FUN = ifelse(
        feature %in% "length",
        colSums,
        function(x){max[1]}
      )
    )
    save(scd, file = paste0("results/", type, "_", feature, ".Rdata"))
  }
}

system("file_handle.py -f data/matrix_output/male/* data/matrix_output/female/paired_end/* data/matrix_output/female/single_end/*")
# fix cells id
system("perl -pi -e 's/\\S*_(P\\d*_\\d*)_\\S*/\\1/g' data/matrix_output/male/*")
system("perl -pi -e 's/\\S*_(P\\d*_\\d*)_\\S*/\\1/g' data/matrix_output/female/paired_end/*")
system("perl -pi -e 's/\\S*_(P\\d*_\\d*)_\\S*/\\1/g' data/matrix_output/female/single_end/*")
# fix cells id
system("perl -pi -e 's/P1306/P1316/g' data/2017_09_14_Summary_SSEQ_170825.csv")
# combine  2017_09_15_count_genes_D15_P1373_YFV2003_Run2_new.tsv and
# 2017_09_15_count_genes_D15_P1373_YFV2003_new.tsv
count_a <- read.table(
  "data/matrix_output/female/single_end/2017_09_15_count_genes_D15_P1373_YFV2003_Run2_new.tsv",
  h = T)
count_b <- read.table(
  "data/matrix_output/female/single_end/2017_09_15_count_genes_D15_P1373_YFV2003_new.tsv",
  h = T)
count <- count_a + count_b
write.table(count,
  "data/matrix_output/female/single_end/2017_09_15_count_genes_D15_P1373_YFV2003.tsv",
)
system("rm data/matrix_output/female/single_end/2017_09_15_count_genes_D15_P1373_YFV2003_Run2_new.tsv")
system("rm data/matrix_output/female/single_end/2017_09_15_count_genes_D15_P1373_YFV2003_new.tsv")


scd <- scRNAtools::load_data(
  infos = "data/2017_09_14_Summary_SSEQ_170825.csv",
  counts = "data/matrix_output",
  regexp = ".*count_genes.*",
  infos_sep = ","
)
save(scd, file = "results/raw_counts.Rdata")

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
qsub src/pbs/QC/QC.pbs &
/bin/sleep 0.5
done
/bin/sleep 3600
done
")

scRNAtools::QC_load_bootstraps(
  scd = scd,
  paraload_folder = "results/QC/QC_paraload/"
)
hist(scd$getfeature("QC_score"))

scRNAtools::QC_classification(scd)

table(scd$getfeature("QC_good"))

devtools::load_all("../scRNAtools/", reset = T)
check_gene(scd, "CCR7", "sex")
scd_QC <- scd$select(b_cells = scd$getfeature("QC_good") == T)

check_gene(scd_QC, "CCR7", "sex")
