rm(list=ls())
setwd("~/projects/mold/yellow_fever/")
library(scRNAtools)
devtools::load_all("pkg/", reset = T)

system("perl -pi -e 's/[pP](\\d*_\\d*)/P\\1/g' data/Summary_SSEQ.csv")
system("perl -pi -e 's/P1306/P1316/g' data/Summary_SSEQ.csv")
############################## load data ######################################
system("mkdir -p results/tmp")
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

system("mkdir -p results/QC/pca")
load("results/abundance.Rdata")
b_cells <- scd$getfeature("to_QC") |
  scd$getcells %in% c("P1299_1797", "P1299_1896")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "batch", color_name = "clonality",
  tmp_file = "results/tmp/pca_abundance_outliers_tmp.Rdata")
load("results/tmp/pca_abundance_outliers_tmp.Rdata")
scd$select(b_cells = b_cells)$getfeatures[order(pca_data$x)[1:2],]
ggsave(file = "results/QC/pca/pca_abundance_outliers.pdf")

load("results/counts.Rdata")
b_cells <- scd$getfeature("to_QC") |
  scd$getcells %in% c("P1299_1797", "P1299_1896")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "batch", color_name = "clonality",
  tmp_file = "results/tmp/pca_counts_outliers_tmp.Rdata")
load("results/tmp/pca_counts_tmp.Rdata", v = T)
scd$select(b_cells = b_cells)$getfeatures[order(pca_data$y, decreasing = T)[1:2],]
ggsave(file = "results/QC/pca/pca_counts_outliers.pdf")

#### For P1856 data ####
scd <- scRNAtools::load_data_salmon(
  infos = "data/2020_01_22_P1856_Meta.csv",
  counts = "data/salmon_output/P1856/",
  feature = "counts",
  tximport_obj = "results/tmp/tmp_tximport_P1856",
  infos_sep = ","
)
scd$setfeature(
  "sequencing",
  "paired"
)
scd$setfeature(
  "to_QC",
  T
)
save(scd, file = paste0("results/count_P1856.Rdata"))

############################# QC from M #######################################
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
devtools::load_all("pkg", reset = T)
scRNAtools::QC_load_bootstraps(
  scd = scd,
  paraload_folder = "results/QC/QC_paraload/counts",
)
scRNAtools::QC_classification(
  scd = scd,
  quant = 0.80
)
print(table(scd$getfeature("QC_good")))
null_cells <- scd$getcells[rowSums(scd$getcounts) != 0]
null_genes <- scd$getgenes[colSums(scd$getcounts) != 0]
scd <- scd$select(cells = null_cells, genes = null_genes)
save(scd, file = "results/QC/counts_QC_M.Rdata")

counts_features <- scd$getfeatures
load("results/abundance.Rdata")
scd <- scd$select(cells = null_cells, genes = null_genes)
scd <- scdata$new(
  infos = counts_features,
  counts = scd$getcounts,
  v = T
)
save(scd, file = "results/QC/abundance_QC.Rdata")

############################# QC plots for M ##################################

load("results/QC/counts_QC_M.Rdata")

table(scd$getfeature("batch"), scd$getfeature("QC_good"))
table(scd$getfeature("day"), scd$getfeature("QC_good"))

system("rm results/tmp/pca_*_QC_tmp.Rdata")

# PCA on cell quality
load("results/QC/counts_QC_M.Rdata")

b_cells <- scd$getfeature("to_QC")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "QC_good", color_name = "antigen",
  tmp_file = "results/tmp/pca_counts_QC_tmp.Rdata",
  main = "all day"
)

ggsave(file = "results/QC/pca/pca_counts_QC.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pca_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "batch", color_name = "clonality",
    alpha = "QC_good",
    tmp_file = paste0("results/tmp/pca_counts_", day, "QC_tmp.Rdata"),
    main = day
  )
  ggsave(file = paste0("results/QC/pca/pca_counts_QC_", day, ".pdf"))
}

load("results/QC/counts_QC_M.Rdata")
b_cells = scd$getfeature('QC_good') %in% T
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "batch", color_name = "clonality",
  tmp_file = "results/tmp/pca_counts_QC_good_tmp.Rdata",
  main = "all day"
)
ggsave(file = "results/QC/pca/pca_counts_QC_good.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pca_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "batch", color_name = "clonality",
    tmp_file = paste0("results/tmp/pca_counts_", day, "QC_good_tmp.Rdata"),
    main = day
  )
  ggsave(file = paste0("results/QC/pca/pca_counts_QC_good_", day, ".pdf"))
}


system("rm results/tmp/pCMF_counts_QC_*_tmp.Rdata")
system("rm results/tmp/pCMF_counts_QC_*_good_tmp.Rdata")

system("mkdir -p results/QC/pcmf")
# pCMF on cell quality
devtools::load_all("../scRNAtools/", reset = T)
load("results/QC/counts_QC.Rdata")
b_cells = scd$getfeature('to_QC')
scRNAtools::pCMF_plot(
  scd$select(b_cells = b_cells), color = "batch", color_name = "clonality",
  alpha = "QC_good",
  tmp_file = "results/tmp/pCMF_counts_QC_tmp.Rdata",
  main = "all day",
  ncores = 11
)
ggsave(file = "results/QC/pcmf/pcmf_counts_QC.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pCMF_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "batch", color_name = "clonality",
    alpha = "QC_good",
    tmp_file = paste0("results/tmp/pCMF_counts_", day, "QC_tmp.Rdata"),
    main = day,
    ncores = 11
  )
  ggsave(file = paste0("results/QC/pcmf/pcmf_counts_QC_", day, ".pdf"))
}

load("results/QC/counts_QC.Rdata")
b_cells = scd$getfeature('QC_good') %in% T
scRNAtools::pCMF_plot(
  scd$select(b_cells = b_cells), color = "batch", color_name = "clonality",
  tmp_file = "results/tmp/pCMF_counts_QC_good_tmp.Rdata",
  main = "all day",,
  ncores = 11
)
ggsave(file = "results/QC/pcmf/pcmf_counts_QC_good.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pCMF_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "batch", color_name = "clonality",
    tmp_file = paste0("results/tmp/pCMF_counts_", day, "QC_good_tmp.Rdata"),
    main = day,
    ncores = 11
  )
  ggsave(file = paste0("results/QC/pcmf/pcmf_counts_QC_good_", day, ".pdf"))
}

# cells effect normalization
devtools::load_all("../scRNAtools/", reset = T)
# system("rm results/tmp/normalization_tmp.Rdata")
load("results/QC/counts_QC_M.Rdata")

b_cells = scd$getfeature('QC_good') %in% T
system("rm results/tmp/normalization_tmp.Rdata")
scd <- normalize(
  scd = scd,
  b_cells = scd$getfeature("QC_good") %in% T,
  method = "SCnorm",
  cpus = 4,
  tmp_file = "results/tmp/normalization_tmp.Rdata"
)
save(scd, file = "results/QC/cells_counts_QC_M.Rdata")
system("~/scripts/sms.sh \" cell normalization done\"")

system("rm results/tmp/pca_cells_counts_QC_good_tmp.Rdata")
for (day in c("D15", "D136", "D593"))
  system(paste0("rm results/tmp/pca_cells_counts_", day, "QC_good_tmp.Rdata"))

load("results/QC/cells_counts_QC_M.Rdata")
b_cells = scd$getfeature('QC_good') %in% T
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "batch", color_name = "clonality",
  tmp_file = "results/tmp/pca_cells_counts_QC_good_tmp.Rdata",
  main = "all day"
)
ggsave(file = "results/QC/pca/pca_cells_counts_QC_good.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pca_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "batch", color_name = "clonality",
    tmp_file = paste0("results/tmp/pca_cells_counts_", day, "QC_good_tmp.Rdata"),
    main = day
  )
  ggsave(file = paste0("results/QC/pca/pca_cells_counts_QC_good_", day, ".pdf"))
}

system("rm results/tmp/pCMF_cells_counts_QC_good_tmp.Rdata")
for (day in c("D15", "D136", "D593"))
  system(paste0("rm results/tmp/pCMF_cells_counts_", day, "QC_good_tmp.Rdata"))

load("results/QC/cells_counts_QC.Rdata")
b_cells = scd$getfeature('QC_good') %in% T
scRNAtools::pCMF_plot(
  scd$select(b_cells = b_cells), color = "batch", color_name = "clonality",
  tmp_file = "results/tmp/pCMF_cells_counts_QC_good_tmp.Rdata",
  main = "all day",,
  ncores = 11
)
ggsave(file = "results/QC/pcmf/pcmf_cells_counts_QC_good.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pCMF_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "batch", color_name = "clonality",
    tmp_file = paste0("results/tmp/pCMF_cells_counts_", day, "QC_good_tmp.Rdata"),
    main = day,
    ncores = 11
  )
  ggsave(file = paste0("results/QC/pcmf/pcmf_cells_counts_QC_good_", day, ".pdf"))
}

# batch & cells effect normalization
devtools::load_all("../scRNAtools/", reset = T)
load("results/QC/cells_counts_QC_M.Rdata")
system("rm results/tmp/normalization_cells_mnnCorrect_tmp.Rdata")
b_cells = scd$getfeature('QC_good') %in% T
scd <- normalize(
  scd = scd,
  b_cells = scd$getfeature("QC_good") %in% T,
  method = "mnnCorrect",
  cpus = 5,
  tmp_file = paste0("results/tmp/normalization_cells_mnnCorrect_tmp.Rdata")
)
save(scd, file = "results/QC/CB_counts_QC_M.Rdata")
system("~/scripts/sms.sh \" batch normalization done\"")


system("rm results/tmp/pca_CB_counts_QC_good_tmp.Rdata")
for (day in c("D15", "D136", "D593"))
  system(paste0("rm results/tmp/pca_CB_counts_", day, "QC_good_tmp.Rdata"))

load("results/QC/CB_counts_QC_M.Rdata")
b_cells = scd$getfeature('QC_good') %in% T
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "batch", color_name = "clonality",
  tmp_file = "results/tmp/pca_CB_counts_QC_good_tmp.Rdata",
  main = "all day"
)
ggsave(file = "results/QC/pca/pca_CB_counts_QC_good.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pca_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "batch", color_name = "clonality",
    tmp_file = paste0("results/tmp/pca_CB_counts_", day, "QC_good_tmp.Rdata"),
    main = day
  )
  ggsave(file = paste0("results/QC/pca/pca_CB_counts_QC_good_", day, ".pdf"))
}

system("rm results/tmp/pCMF_CB_counts_QC_good_tmp.Rdata")
for (day in c("D15", "D136", "D593"))
  system(paste0("rm results/tmp/pCMF_CB_counts_", day, "QC_good_tmp.Rdata"))

load("results/QC/CB_counts_QC.Rdata")
b_cells = scd$getfeature('QC_good') %in% T
scRNAtools::pCMF_plot(
  scd$select(b_cells = b_cells), color = "batch", color_name = "clonality",
  tmp_file = "results/tmp/pCMF_CB_counts_QC_good_tmp.Rdata",
  main = "all day",,
  ncores = 11
)
ggsave(file = "results/QC/pcmf/pcmf_CB_counts_QC_good.pdf")
for (day in c("D15", "D136", "D593")) {
  scRNAtools::pCMF_plot(
    scd$select(b_cells = b_cells & scd$getfeature("day") %in% day),
    color = "batch", color_name = "clonality",
    tmp_file = paste0("results/tmp/pCMF_CB_counts_", day, "QC_good_tmp.Rdata"),
    main = day,
    ncores = 11
  )
  ggsave(file = paste0("results/QC/pcmf/pcmf_CB_counts_QC_good_", day, ".pdf"))
}


############################## weird_D15 analysis #############################
setwd("~/projects/yellow_fever/")
library(scRNAtools)
devtools::load_all("pkg/", reset = T)
load("results/QC/CB_counts_QC_M.Rdata")

weird_D15 <- c("P1299_1105", "P1299_1106", "P1299_1111", "P1299_1112", "P1299_1117", "P1299_1129", "P1299_1133", "P1299_1150", "P1299_1151", "P1299_1185", "P1299_1222", "P1299_1263", "P1299_1284", "P1299_1297", "P1299_1299", "P1299_1313", "P1299_1328", "P1299_1336", "P1299_1345", "P1299_1356", "P1299_1364", "P1299_1371", "P1299_1390", "P1299_1397", "P1299_1404", "P1299_1416", "P1299_1429", "P1299_1432", "P1299_1437", "P1299_1445", "P1299_1457", "P1299_1465", "P1299_1466", "P1299_1473", "P1299_1478", "P1299_1770", "P1299_1772", "P1299_1781", "P1299_1795", "P1299_1802", "P1299_1803", "P1299_1810", "P1299_1818", "P1299_1819", "P1299_1826", "P1299_1838", "P1299_1843", "P1299_1847", "P1299_1850", "P1299_1861", "P1299_1881", "P1299_1882", "P1299_1884", "P1299_1908", "P1299_1913", "P1299_1921", "P1299_1922", "P1299_1928", "P1299_1949", "P1299_2012", "P1299_2017", "P1299_2035", "P1299_2052", "P1299_2054", "P1299_2056")

table(weird_D15 %in% scd$getcells)
scd$setfeature("is_weird", scd$getcells %in% weird_D15)
b_cells <- scd$getfeature("QC_good") %in% T
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "is_weird", color_name = "antigen",
  tmp_file = "results/tmp/pca_norm_counts_QC_tmp.Rdata")

b_cells <- scd$getfeature("to_QC") & !scd$getfeature("is_weird") &
  scd$getfeature("sex") %in% "M"
pca_out <- ade4::dudi.pca(log1p(scd$select(b_cells = b_cells)$getcounts),
                    scan = F,
                    nf = scd$getncells - 1,
                    scale = TRUE,
                    center = TRUE)
prop_of_var <- 100 * pca_out$eig[1:3] / sum(pca_out$eig)
ggplot(data = data.frame(x = pca_out$l1$RS1,
                         y = pca_out$l1$RS2,
                         color = as.vector(scd$select(b_cells = b_cells)$getfeature("QC_good")) %in% "TRUE" &
                         scd$select(b_cells = b_cells)$getfeature("cell_number") > 0),
  aes(x = x, y = y, color = color)) +
  labs(x = paste0("first axe: ",
                    round(prop_of_var[1], digit = 2),
                    "% of variance"),
       y = paste0("second axe: ",
                    round(prop_of_var[2], digit = 2),
                    "% of variance"),
       color = "QC pass"
       ) +
  scale_color_manual(values=c("#DF1000", "#2E3191")) +
  geom_point() +
  theme_bw()
ggsave(file = paste0("results/QC/pca_Male_QC.pdf"))

system("rm results/tmp/pca_norm_counts_QC_tmp.Rdata")
b_cells <- scd$getfeature("QC_good") %in% T
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "day", color_name = "day",
  tmp_file = "results/tmp/pca_norm_counts_QC_tmp.Rdata")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "batch", color_name = "clonality",
  tmp_file = "results/tmp/pca_norm_counts_QC_tmp.Rdata")


scRNAtools::pCMF_plot(
  scd$select(b_cells = b_cells), color = "is_weird", color_name = "antigen",
  tmp_file = "results/tmp/pCMF_norm_counts_QC_tmp.Rdata")
devtools::load_all("../scRNAtools/", reset = T)
scRNAtools::pCMF_plot(
  scd$select(b_cells = b_cells), color = "day", color_name = "day",
  tmp_file = "results/tmp/pCMF_norm_counts_QC_tmp.Rdata")
scRNAtools::pCMF_plot(
  scd$select(b_cells = b_cells), color = "batch", color_name = "clonality",
  tmp_file = "results/tmp/pCMF_norm_counts_QC_tmp.Rdata")

devtools::install_github('theislab/kBET')
library(kBET)
load(file = "results/QC/counts_QC.Rdata")
b_cells <- scd$getfeature("to_QC") %in% T
batch_estimate <- list()
for (day in c("D15", "D136", "D593")) {
  print(day)
  batch_estimate[[day]] <- kBET(
    scd$select(b_cells = b_cells)$getcounts,
    scd$select(b_cells = b_cells)$getfeature("batch")
  )
  print(batch_estimate[[day]]$summary)
}
b_cells <- scd$getfeature("QC_good") %in% T
batch_estimate_qc <- list()
for (day in c("D15", "D136", "D593")) {
  print(day)
  batch_estimate_qc[[day]] <- kBET(
    scd$select(b_cells = b_cells)$getcounts,
    scd$select(b_cells = b_cells)$getfeature("batch")
  )
  print(batch_estimate_qc[[day]]$summary)
}
load(file = "results/QC/norm_counts_QC.Rdata")
b_cells <- scd$getfeature("QC_good") %in% T
batch_estimate_qc_norm <- list()
for (day in c("D15", "D136", "D593")) {
  print(day)
  batch_estimate_qc_norm[[day]] <- kBET(
    scd$select(b_cells = b_cells)$getcounts,
    scd$select(b_cells = b_cells)$getfeature("batch")
  )
  print(batch_estimate_qc_norm[[day]]$summary)
}
system("~/scripts/sms.sh \"normalization done\"")
################################################################################
# QC of F data

rm(list=ls())
setwd("~/projects/yellow_fever/")
library(scRNAtools)
devtools::load_all("pkg/", reset = T)
load("results/QC/counts_QC_M.Rdata")
bad_F_cells <- paste0("P1292_", 1097:1192)
scd <- scd$select(b_cells = !( scd$getfeature("id") %in% bad_F_cells ))

F_QC_score <- ifelse(scd$getfeature("sex") %in% "F",
  0.5,
  scd$getfeature("QC_score")
)
F_QC_score <- ifelse(scd$getfeature("cell_number") == 0 &
                     scd$getfeature("sex") %in% "F",
  0,
  F_QC_score
)
F_to_QC <- ifelse(scd$getfeature("sex") %in% "F",
  TRUE,
  scd$getfeature("to_QC")
)

scd$setfeature("QC_score",
  F_QC_score
)
scd$setfeature("to_QC",
  F_to_QC
)
scRNAtools::QC_classification(
  scd = scd,
  is_blank = scd$getfeature("QC_score") < 0.5
)
save(scd, file = "results/QC/counts_QC_F.Rdata")

b_cells = scd$getfeature('sex') %in% "F"
summary(scd$select(b_cells = b_cells)$getfeature("QC_good"))


ggplot(data = data.frame(QC_good = (scd$select(b_cells = b_cells)$getfeature("QC_good")),
                         QC_score = (scd$select(b_cells = b_cells)$getfeature("QC_score")))) +
  geom_point(aes(x = QC_good, y = QC_score))


table(scd$getfeature('sex'), scd$getfeature("QC_good"))

scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "QC_good", color_name = "antigen",
  tmp_file = "results/tmp/pca_counts_F_tmp.Rdata",
  main = "all day"
)

b_cells = scd$getfeature('sex') %in% "F" & scd$getfeature("QC_good") %in% T
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "batch", color_name = "clonality",
  tmp_file = "results/tmp/pca_counts_QC_F_tmp.Rdata",
  main = "all day"
)

scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "day", color_name = "day",
  tmp_file = "results/tmp/pca_counts_QC_F_tmp.Rdata",
  main = "all day"
)

# cells effect normalization
load("results/QC/counts_QC.Rdata")
b_cells = scd$getfeature('sex') %in% "F" & scd$getfeature("QC_good") %in% T
for (day in c("D15", "D90")) {
  system(paste0("rm results/tmp/normalization_", day, "_F_tmp.Rdata"))
  scd <- normalize(
    scd = scd,
    b_cells = b_cells & scd$getfeature("day") %in% day,
    method = "SCnorm",
    cpus = 4,
    tmp_file = paste0("results/tmp/normalization_", day, "_F_tmp.Rdata")
  )
}
save(scd, file = "results/QC/cells_counts_QC_F.Rdata")
load(file = "results/QC/cells_counts_QC_F.Rdata")

scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "batch", color_name = "clonality",
  tmp_file = "results/tmp/pca_cells_QC_F_tmp.Rdata",
  main = "all day"
)
ggsave(file = "results/QC/pca/pca_cells_counts_F_QC_good.pdf")

# batch effect normalization
load("results/QC/cells_counts_QC_F.Rdata")
devtools::load_all("../scRNAtools/", reset = T)
b_cells = scd$getfeature('sex') %in% "F" & scd$getfeature("QC_good") %in% T
for (day in c("D15", "D90")) {
  system(paste0("rm results/tmp/normalization_cells_combat_,", day, "_F_tmp.Rdata"))
  table(scd$select(b_cells = b_cells & scd$getfeature("day") %in% day)$getfeature("batch"))
  scd$select(b_cells = b_cells & scd$getfeature("day") %in% day)$getncells
  scd <- normalize(
    scd = scd,
    b_cells = b_cells & scd$getfeature("day") %in% day,
    method = "mnnCorrect",
    cpus = 5,
    tmp_file = paste0("results/tmp/normalization_cells_combat_,", day, "_F_tmp.Rdata")
  )
}

load("results/QC/cells_counts_QC_F.Rdata")
devtools::load_all("../scRNAtools/", reset = T)
b_cells = scd$getfeature('sex') %in% "F" & scd$getfeature("QC_good") %in% T
scd <- normalize(
  scd = scd,
  b_cells = b_cells,
  method = "mnnCorrect",
  cpus = 5,
  tmp_file = paste0("results/tmp/normalization_cells_mnnCorrect_F_tmp.Rdata")
)
traceback()
save(scd, file = "results/QC/CB_counts_QC_F.Rdata")


load("results/QC/CB_counts_QC_F.Rdata")
scd_mnn <- scd
load("results/QC/cells_counts_QC_F.Rdata")
b_cells = scd$getfeature('sex') %in% "F" & scd$getfeature("QC_good") %in% T
scd <- normalize(
  scd = scd,
  b_cells = b_cells,
  method = "ComBat",
  cpus = 5,
  tmp_file = paste0("results/tmp/normalization_cells_ComBat_F_tmp.Rdata")
)
save(scd, file = "results/QC/CB_counts_QC_F.Rdata")


system("rm results/tmp/pca_CB_QC_F_tmp.Rdata")
scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "batch", color_name = "clonality",
  tmp_file = "results/tmp/pca_CB_QC_F_tmp.Rdata",
  main = "all day"
)
ggsave(file = "results/QC/pca/pca_CB_counts_F_QC_good.pdf")
system("~/scripts/sms.sh \"normalization done\"")

################################# F and M data set ############################
rm(list=ls())
setwd("~/projects/yellow_fever/")
devtools::load_all("pkg/", reset = T)
bad_F_cells <- paste0("P1292_", 1097:1192)

# merge count data
load("results/QC/counts_QC_M.Rdata")
scd <- scd$select(b_cells = !( scd$getfeature("id") %in% bad_F_cells ))
infos_M <- scd$getfeatures
counts_M <- scd$getcounts
load("results/QC/counts_QC_F.Rdata")
b_cells_F <- scd$getfeature("sex") %in% "F"
infos_M[b_cells_F, ] <- scd$select(b_cells = b_cells_F)$getfeatures
counts_M[b_cells_F, ] <- scd$select(b_cells = b_cells_F)$getcounts
scd <- scdata$new(
  infos = infos_M,
  counts = counts_M
)
bad_F_cells <- paste0("P1292_", 1097:1192)
scd <- scd$select(b_cells = !( scd$getfeature("id") %in% bad_F_cells ))
save(scd, file = "results/QC/counts_QC.Rdata")

# merge count data normalize for cell effect
load("results/QC/cells_counts_QC_M.Rdata")
scd <- scd$select(b_cells = !( scd$getfeature("id") %in% bad_F_cells ))
infos_M <- scd$getfeatures
counts_M <- scd$getcounts
load("results/QC/cells_counts_QC_F.Rdata")
b_cells_F <- scd$getfeature("sex") %in% "F"
infos_M[b_cells_F, ] <- scd$select(b_cells = b_cells_F)$getfeatures
counts_M[b_cells_F, ] <- scd$select(b_cells = b_cells_F)$getcounts
scd <- scdata$new(
  infos = infos_M,
  counts = counts_M
)
bad_F_cells <- paste0("P1292_", 1097:1192)
scd <- scd$select(b_cells = !( scd$getfeature("id") %in% bad_F_cells ))
save(scd, file = "results/QC/cells_counts_QC.Rdata")

# merge count data normalize for cell effect and batch effect
load("results/QC/CB_counts_QC_M.Rdata")
scd <- scd$select(b_cells = !( scd$getfeature("id") %in% bad_F_cells ))
infos_M <- scd$getfeatures
counts_M <- scd$getcounts
load("results/QC/CB_counts_QC_F.Rdata")
b_cells_F <- scd$getfeature("sex") %in% "F"
infos_M[b_cells_F, ] <- scd$select(b_cells = b_cells_F)$getfeatures
counts_M[b_cells_F, ] <- scd$select(b_cells = b_cells_F)$getcounts
scd <- scdata$new(
  infos = infos_M,
  counts = counts_M
)
save(scd, file = "results/QC/CB_counts_QC.Rdata")

load(file = "results/QC/cells_counts_QC.Rdata")
outlier <- c("P1373_1013", "P1373_1023", "P1373_1043", "P1373_1058", "P1373_1025")
b_cells <- scd$getfeature("QC_good") %in% T & scd$getfeature("sex") %in% "F"
data <- scd$select(b_cells = b_cells)$getfeatures
ggplot(data = data, aes(x = QC_score, y = QC_score, color = id %in% outlier)) +
  geom_point()

summary(data[data$id %in% outlier,])


###############################################################################
############################ Jeff normalization ###############################

library("tidyverse")
setwd("~/projects/mold/yellow_fever/")
devtools::load_all("pkg/", reset = T)
system("perl -pi -e 's/[pP](\\d*_\\d*)/P\\1/g' data/2019_05_09_Samples_for_Matrix_InVivo.csv")
samples <- read_csv("data/2019_05_09_Samples_for_Matrix_InVivo.csv")
genes <- read_csv("data/2019_05_09_Jeff_overlap_MB7_FAll.txt")
load(file = "results/counts.Rdata")

bad_F_cells <- paste0("P1292_", 1097:1192)
good_F_cells <- paste0("P1373_", 1097:1192)
samples <- samples[!(samples$Sample_ID %in% bad_F_cells), ]
samples <- rbind(samples, data.frame(Sample_ID = good_F_cells))


length(samples$Sample_ID)
table(scd$getfeature("id") %in% samples$Sample_ID)
samples$Sample_ID[!samples$Sample_ID %in% scd$getfeature("id")]
system("rm results/tmp/normalization_scnorm_jeff_selection_tmp.Rdata")
scd_norm <- normalize(
  scd = scd$select(b_cells = scd$getfeature("id") %in% samples$Sample_ID,
                   genes = genes$x),
  b_cells = TRUE,
  method = "SCnorm",
  cpus = 4,
  tmp_file = paste0("results/tmp/normalization_scnorm_jeff_selection_tmp.Rdata")
)

save(scd_norm,
     file = "results/QC/cells_counts_Jeff.Rdata"
)

load("results/QC/cells_counts_Jeff.Rdata")
load("results/cycling/cells_counts_QC_cycling.Rdata")

infos_M <- scd$getfeatures
rownames(infos_M) <- infos_M$id
counts_M <- scd_norm$getcounts
infos_M <- infos_M[rownames(infos_M) %in% rownames(counts_M), ]
counts_M <- counts_M[rownames(counts_M) %in% rownames(infos_M), ]
dim(infos_M)
dim(counts_M)
normalized_counts <- cbind(infos_M[match(rownames(infos_M), rownames(counts_M)), ],
                           counts_M
)
normalized_counts <- apply(normalized_counts, 2, as.vector)
normalized_counts <- t(normalized_counts)
write.csv(
  normalized_counts,
  file = paste0("results/QC/cells_counts_Jeff.csv")
)


load(file = "results/QC/cells_counts_Jeff.Rdata")
scRNAtools::pca_plot(
  scd_norm, color = "day", color_name = "clonality",
  shape = "sex",
  tmp_file = "results/tmp/pca_cells_counts_Jeff.Rdata")
ggsave(file = "results/QC/pca/pca_cells_counts_Jeff.pdf")

scRNAtools::mds_plot(
  scd_norm, color = "day", color_name = "clonality",
  shape = "sex",
  tmp_file = "results/tmp/mds_cells_counts_Jeff.Rdata")
ggsave(file = "results/QC/pca/mds_cells_counts_Jeff.pdf")

library("umap")
scd_umap <- umap(ascb(scd_norm$getcounts))
ggplot(data.frame(
        x = scd_umap$layout[, 1],
        y = scd_umap$layout[, 2],
        day = scd_norm$getfeature("day"),
        sex = scd_norm$getfeature("sex"),
        antigen = scd_norm$getfeature("antigen"),
        batch = scd_norm$getfeature("batch")
      )
  ) +
  geom_point(aes(x = x,
                 y = y,
                 color = as.factor(batch),
                 shape = day
                )) +
  theme_bw()
ggsave(file = "results/QC/pca/umap_cells_counts_Jeff.pdf")

system("rm results/tmp/normalization_jeff_all_tmp.Rdata")
scd_norm <- normalize(
  scd = scd$select(b_cells = scd$getfeature("id") %in% samples$Sample_ID),
  b_cells = TRUE,
  method = "SCnorm",
  cpus = 1,
  tmp_file = paste0("results/tmp/normalization_jeff_all_tmp.Rdata")
)

save(scd_norm,
     file = "results/QC/cells_counts_all_Jeff.Rdata"
)

load("results/cycling/cells_counts_QC_cycling.Rdata")
load("results/QC/cells_counts_all_Jeff.Rdata")

infos_M <- scd$getfeatures
rownames(infos_M) <- infos_M$id
counts_M <- scd_norm$getcounts
infos_M <- infos_M[rownames(infos_M) %in% rownames(counts_M), ]
counts_M <- counts_M[rownames(counts_M) %in% rownames(infos_M), ]

dim(infos_M)
dim(counts_M)
normalized_counts <- cbind(infos_M[match(rownames(infos_M), rownames(counts_M)), ],
                           counts_M
)
normalized_counts <- apply(normalized_counts, 2, as.vector)
normalized_counts <- t(normalized_counts)

write.csv(
  normalized_counts,
  file = paste0("results/QC/cells_counts_all_Jeff.csv")
)

load(file = "results/QC/cells_counts_all_Jeff.Rdata")
scRNAtools::pca_plot(
  scd_norm, color = "day", color_name = "clonality",
  shape = "sex",
  tmp_file = "results/tmp/pca_cells_counts_all_Jeff.Rdata")
ggsave(file = "results/QC/pca/pca_cells_counts_all_Jeff.pdf")

scRNAtools::mds_plot(
  scd_norm, color = "day", color_name = "clonality",
  shape = "sex",
  tmp_file = "results/tmp/mds_cells_counts_all_Jeff.Rdata")
ggsave(file = "results/QC/pca/mds_cells_counts_all_Jeff.pdf")

library("umap")
scd_umap <- umap(ascb(scd_norm$getcounts))
ggplot(data.frame(
        x = scd_umap$layout[, 1],
        y = scd_umap$layout[, 2],
        day = scd_norm$getfeature("day"),
        sex = scd_norm$getfeature("sex"),
        antigen = scd_norm$getfeature("antigen")
      )
  ) +
  geom_point(aes(x = x,
                 y = y,
                 color = day,
                 shape = sex
                )) +
  theme_bw()
ggsave(file = "results/QC/pca/umap_cells_counts_all_Jeff.pdf")


################################################################################
# QC of P1856 data

rm(list=ls())
setwd("~/projects/mold/yellow_fever/")
library(scRNAtools)
devtools::load_all("pkg/", reset = T)
load("results/QC/counts_QC_F.Rdata")
bad_F_cells <- paste0("P1292_", 1097:1192)
scd_old <- scd
load("results/count_P1856.Rdata")
scd$setfeature("QC_score", 0.5)
scd$setfeature("to_QC", T) 
scd$setfeature("batch", "P1856") 
scd$setfeature("cell_number", 1) 
colnames(scd_old$getfeatures)
for( feat in colnames(scd_old$getfeatures) ){
  if (!(feat %in% colnames(scd$getfeatures))){
    scd$setfeature(feat, NA)
  }
}
scd_old$setfeature("donor", NA)

scd <- scdata$new(
  infos = rbind(
    scd_old$getfeatures,
    scd$getfeatures
  ),
  counts = rbind(
    scd_old$getcounts,
    scd$getcounts[, colnames(scd$getcounts) %in% colnames(scd_old$getcounts)]
  )
)

table(scd$getfeature("to_QC"), scd$getfeature("batch"))
table(scd$getfeature("to_QC"), scd$getfeature("sex"))

scRNAtools::QC_classification(
  scd = scd,
  is_blank = scd$getfeature("QC_score") < 0.5
)
save(scd, file = "results/QC/counts_QC_P1856.Rdata")

b_cells <- scd$getfeature("batch") %in% "P1856"
summary(scd$select(b_cells = b_cells)$getfeature("QC_good"))


ggplot(data = data.frame(QC_good = (scd$select(b_cells = b_cells)$getfeature("QC_good")),
                         QC_score = (scd$select(b_cells = b_cells)$getfeature("QC_score")))) +
  geom_point(aes(x = QC_good, y = QC_score))

scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "QC_good", color_name = "antigen",
  tmp_file = "results/tmp/pca_counts_P1856_tmp.Rdata",
  main = "all day"
)

# cells effect normalization
load("results/QC/counts_QC_P1856.Rdata")
b_cells = scd$getfeature("QC_good") %in% T
scd <- normalize(
  scd = scd,
  b_cells = b_cells,
  method = "SCnorm",
  cpus = 4,
  tmp_file = paste0("results/tmp/normalization_P1856_tmp.Rdata")
)
save(scd, file = "results/QC/cells_counts_QC_P1856.Rdata")
load(file = "results/QC/cells_counts_QC_P1856.Rdata")

scRNAtools::pca_plot(
  scd$select(b_cells = b_cells), color = "clonality", color_name = "clonality",
  tmp_file = "results/tmp/pca_cells_QC_P1856_tmp.Rdata",
  main = "all day"
)
ggsave(file = "results/QC/pca/pca_cells_counts_P1856_QC_good.pdf")

# batch effect normalization
load(file = "results/QC/cells_counts_QC_P1856.Rdata")
devtools::load_all("../scRNAtools/", reset = T)
b_cells = scd$getfeature("QC_good") %in% T

scd <- normalize(
  scd = scd,
  b_cells = b_cells,
  method = "ComBat",
  cpus = 5,
  tmp_file = paste0("results/tmp/normalization_cells_ComBat_P1856_tmp.Rdata")
)
save(scd, file = "results/QC/CB_counts_QC_1856.Rdata")

data_counts <- scd$select(b_cells = b_cells & scd$getfeature("batch") %in% "P1856")$getcounts
write.csv(data_counts, "results/QC/CB_counts_QC_1856.csv")

data_infos <- scd$select(b_cells = b_cells & scd$getfeature("batch") %in% "P1856")$getfeatures
write.csv(data_infos, "results/QC/CB_counts_QC_1856_infos.csv")
