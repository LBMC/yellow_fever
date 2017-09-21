setwd("~/projects/yellow_fever/")
library(scRNAtools)

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


devtools::load_all("../scRNAtools/", reset = T)
scd <- scRNAtools::load_data(
  infos = "data/2017_09_14_Summary_SSEQ_170825.csv",
  counts = "data/matrix_output",
  regexp = ".*count_genes.*"
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
--conf src/pbs/QC/paraload.csv
")

# launch paralod clients
system("
bin/paraload --client --port 13469 --host pbil-deb
")

system("
while [ $(ps -u modolo | grep paraload | wc -l) -gt 0 ]
do
stat | wc -l
iter=$(echo 1000 - $(stat | wc -l) | bc)
for ((i = 1;i <= $iter;i += 1))
do
qsub src/pbs/QC/QC.pbs &
sleep 0.5
done
sleep 3600
done
")
