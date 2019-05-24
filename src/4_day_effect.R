rm(list = ls())
install.packages("testthat")
setwd("~/projects/mold/yellow_fever")
devtools::load_all("pkg/", reset = T)
load("results/cycling/cells_counts_QC_cycling.Rdata")

### Male donor
days <- c("D15", "D136", "D593")
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("day") %in% days &
  scd$getfeature("sex") %in% "M"

days <- ifelse(scd$getfeature("day") %in% "D15", "D15",
              "D100+")
scd$setfeature("days", days)


table(scd$select(b_cells = b_cells)$getfeature("days"),
      scd$select(b_cells = b_cells)$getfeature("DEA_cell_type"))

system(paste0("mkdir -p results/days/"))
mbatch_days_DEA <- DEA(
  scd = scd,
  formula_null = "y ~ (1|batch) + DEA_cell_type",
  formula_full = "y ~ (1|batch) + DEA_cell_type + days",
  b_cells = b_cells,
  cpus = 11,
  v = T,
  folder_name = paste0("results/days/mbatch_DEA_cell_type_days_DEA")
)
save(
  mbatch_days_DEA,
  file = paste0("results/days/mbatch_DEA_cell_type_days_DEA.Rdata")
)
system("~/scripts/sms.sh \"DEA done\"")
print(table(is.na(mbatch_days_DEA$padj)))
print(table(mbatch_days_DEA$padj < 0.05))
write.csv(
  mbatch_days_DEA,
  file = paste0("results/days/mbatch_DEA_cell_type_days_DEA.csv")
)

devtools::load_all("pkg/", reset = T)
ud_mbatch_days_DEA <- up_down(
  folder_name = paste0("results/days/mbatch_DEA_cell_type_days_DEA")
)
test_genes <- c("TSC22D3",
                "CXCR4",
                "NFKBIA",
                "MKI67"
                )
head(ud_mbatch_days_DEA)
ud_mbatch_days_DEA[test_genes, ]
write.csv(
  ud_mbatch_days_DEA,
  file = paste0("results/days/ud_mbatch_DEA_cell_type_days_DEA.csv")
)


### Male donor A2
days <- c("D15", "D136", "D593")
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("day") %in% days &
  scd$getfeature("sex") %in% "M" &
  scd$getfeature("antigen") %in% "A2"

days <- ifelse(scd$getfeature("day") %in% "D15", "D15",
              "D100+")
scd$setfeature("days", days)


table(scd$select(b_cells = b_cells)$getfeature("days"),
      scd$select(b_cells = b_cells)$getfeature("DEA_cell_type"))

table(scd$select(b_cells = b_cells)$getfeature("batch"),
      scd$select(b_cells = b_cells)$getfeature("DEA_cell_type"))

system(paste0("mkdir -p results/days/"))
system(paste0("rm -R results/days/mbatch_DEA_cell_type_days_A2_DEA.Rdata"))
mbatch_days_A2_DEA <- DEA(
  scd = scd,
  formula_null = "y ~ (1|batch) + DEA_cell_type",
  formula_full = "y ~ (1|batch) + DEA_cell_type + days",
  b_cells = b_cells,
  cpus = 11,
  v = T,
  folder_name = paste0("results/days/mbatch_DEA_cell_type_days_A2_DEA")
)
save(
  mbatch_days_A2_DEA,
  file = paste0("results/days/mbatch_DEA_cell_type_days_A2_DEA.Rdata")
)
system("~/scripts/sms.sh \"DEA done\"")
load( file = paste0("results/days/mbatch_DEA_cell_type_days_A2_DEA.Rdata"))
print(table(is.na(mbatch_days_A2_DEA$padj)))
print(table(mbatch_days_A2_DEA$padj < 0.05))
write.csv(
  mbatch_days_A2_DEA,
  file = paste0("results/days/mbatch_DEA_cell_type_days_A2_DEA.csv")
)

devtools::load_all("pkg/", reset = T)
ud_mbatch_days_A2_DEA <- up_down(
  folder_name = paste0("results/days/mbatch_DEA_cell_type_days_A2_DEA")
)

test_genes <- c("TSC22D3",
                "CXCR4",
                "NFKBIA",
                "MKI67"
                )
head(ud_mbatch_days_A2_DEA)
ud_mbatch_days_A2_DEA[test_genes, ]
write.csv(
  ud_mbatch_days_A2_DEA,
  file = paste0("results/days/ud_mbatch_DEA_cell_type_days_A2_DEA.csv")
)

### Male donor B7
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("day") %in% days &
  scd$getfeature("sex") %in% "M" &
  scd$getfeature("antigen") %in% "B7"

days <- ifelse(scd$getfeature("day") %in% "D15", "D15",
              "D100+")
scd$setfeature("days", days)


table(scd$select(b_cells = b_cells)$getfeature("days"),
      scd$select(b_cells = b_cells)$getfeature("DEA_cell_type"))

system(paste0("mkdir -p results/days/"))
system(paste0("rm -R results/days/batch_DEA_cell_type_days_B7_DEA"))
batch_days_B7_DEA <- DEA(
  scd = scd,
  formula_null = "y ~ batch + DEA_cell_type",
  formula_full = "y ~ batch + DEA_cell_type + days",
  b_cells = b_cells,
  cpus = 11,
  v = T,
  folder_name = paste0("results/days/batch_DEA_cell_type_days_B7_DEA")
)
save(
  batch_days_B7_DEA,
  file = paste0("results/days/batch_DEA_cell_type_days_B7_DEA.Rdata")
)
system("~/scripts/sms.sh \"DEA done\"")
load(file = paste0("results/days/batch_DEA_cell_type_days_B7_DEA.Rdata"))
print(table(is.na(batch_days_B7_DEA$padj)))
print(table(batch_days_B7_DEA$padj < 0.05))
write.csv(
  batch_days_B7_DEA,
  file = paste0("results/days/batch_DEA_cell_type_days_B7_DEA.csv")
)

ud_batch_days_B7_DEA <- up_down(
  folder_name = paste0("results/days/batch_DEA_cell_type_days_B7_DEA")
)
test_genes <- c("TSC22D3",
                "CXCR4",
                "NFKBIA",
                "MKI67"
                )
head(ud_batch_days_B7_DEA)
ud_batch_days_B7_DEA[test_genes, ]
write.csv(
  ud_batch_days_B7_DEA,
  file = paste0("results/days/ud_batch_DEA_cell_type_days_B7_DEA.csv")
)

# Male donor B7, D15 vs D136 for MEM
days <- c("D15", "D136")
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("day") %in% days &
  scd$getfeature("sex") %in% "M" &
  scd$getfeature("DEA_cell_type") %in% "MEM"

table(scd$select(b_cells = b_cells)$getfeature("day"),
      scd$select(b_cells = b_cells)$getfeature("batch"))

mbatch_days_B7_MEM_DEA <- DEA(
  scd = scd,
  formula_null = "y ~ (1|batch)",
  formula_full = "y ~ (1|batch) + day",
  b_cells = b_cells,
  cpus = 11,
  v = T,
  folder_name = paste0("results/days/mbatch_DEA_days_B7_MEM_DEA")
)
save(
  mbatch_days_B7_MEM_DEA,
  file = paste0("results/days/mbatch_DEA_days_B7_MEM_DEA.Rdata")
)
system("~/scripts/sms.sh \"DEA done\"")
load(file = paste0("results/days/mbatch_DEA_days_B7_MEM_DEA.Rdata"))
print(table(is.na(mbatch_days_B7_MEM_DEA$padj)))
print(table(mbatch_days_B7_MEM_DEA$padj < 0.05))
write.csv(
  mbatch_days_B7_MEM_DEA,
  file = paste0("results/days/mbatch_DEA_days_B7_MEM_DEA.csv")
)

ud_mbatch_day_B7_MEM_DEA <- up_down(
  folder_name = paste0("results/days/mbatch_DEA_days_B7_MEM_DEA")
)
test_genes <- c("TSC22D3",
                "CXCR4",
                "NFKBIA",
                "MKI67"
                )
head(ud_mbatch_day_B7_MEM_DEA)
ud_mbatch_day_B7_MEM_DEA[test_genes, ]
write.csv(
  ud_mbatch_day_B7_MEM_DEA,
  file = paste0("results/days/ud_mbatch_DEA_days_B7_MEM_DEA.csv")
)

# Male donor B7, D15 vs D136 for EFF
days <- c("D15", "D136")
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("day") %in% days &
  scd$getfeature("sex") %in% "M" &
  scd$getfeature("DEA_cell_type") %in% "EFF"

table(scd$select(b_cells = b_cells)$getfeature("day"),
      scd$select(b_cells = b_cells)$getfeature("batch"))

mbatch_days_B7_EFF_DEA <- DEA(
  scd = scd,
  formula_null = "y ~ (1|batch)",
  formula_full = "y ~ (1|batch) + day",
  b_cells = b_cells,
  cpus = 11,
  v = T,
  folder_name = paste0("results/days/mbatch_DEA_days_B7_EFF_DEA")
)
save(
  mbatch_days_B7_EFF_DEA,
  file = paste0("results/days/mbatch_DEA_days_B7_EFF_DEA.Rdata")
)
system("~/scripts/sms.sh \"DEA done\"")
load(file = paste0("results/days/mbatch_DEA_days_B7_EFF_DEA.Rdata"))
print(table(is.na(mbatch_days_B7_EFF_DEA$padj)))
print(table(mbatch_days_B7_EFF_DEA$padj < 0.05))
write.csv(
  mbatch_days_B7_EFF_DEA,
  file = paste0("results/days/mbatch_DEA_days_B7_EFF_DEA.csv")
)

ud_mbatch_day_B7_EFF_DEA <- up_down(
  folder_name = paste0("results/days/mbatch_DEA_days_B7_EFF_DEA")
)
test_genes <- c("TSC22D3",
                "CXCR4",
                "NFKBIA",
                "MKI67"
                )
head(ud_mbatch_day_B7_EFF_DEA)
ud_mbatch_day_B7_EFF_DEA[test_genes, ]
write.csv(
  ud_mbatch_day_B7_EFF_DEA,
  file = paste0("results/days/ud_mbatch_DEA_days_B7_EFF_DEA.csv")
)


### Female donor

days <- c("D15", "D90")
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("day") %in% days &
  scd$getfeature("sex") %in% "F"


table(scd$select(b_cells = b_cells)$getfeature("day"),
      scd$select(b_cells = b_cells)$getfeature("DEA_cell_type"))

table(scd$select(b_cells = b_cells)$getfeature("day"),
      scd$select(b_cells = b_cells)$getfeature("batch"))

system(paste0("mkdir -p results/days/"))
batch_day_DEA <- DEA(
  scd = scd,
  formula_null = "y ~ batch + DEA_cell_type",
  formula_full = "y ~ batch + DEA_cell_type + day",
  b_cells = b_cells,
  cpus = 11,
  v = T,
  folder_name = paste0("results/days/batch_DEA_cell_type_day_F_DEA")
)
save(
  batch_day_DEA,
  file = paste0("results/days/batch_DEA_cell_type_day_F_DEA.Rdata")
)
system("~/scripts/sms.sh \"DEA done\"")
load(file = paste0("results/days/batch_DEA_cell_type_day_F_DEA.Rdata"))
print(table(is.na(batch_day_DEA$padj)))
print(table(batch_day_DEA$padj < 0.05))
write.csv(
  batch_day_DEA,
  file = paste0("results/days/batch_DEA_cell_type_day_F_DEA.csv")
)

devtools::load_all("pkg/", reset = T)
ud_batch_day_F_DEA <- up_down(
  folder_name = paste0("results/days/batch_DEA_cell_type_day_F_DEA")
)
test_genes <- c("TSC22D3",
                "CXCR4",
                "NFKBIA",
                "MKI67"
                )
head(ud_batch_day_F_DEA)
ud_batch_day_F_DEA[test_genes, ]
ud_batch_day_F_DEA[batch_day_DEA$gene[batch_day_DEA$padj < 0.05], ]
write.csv(
  ud_batch_day_F_DEA,
  file = paste0("results/days/ud_batch_DEA_cell_type_day_F_DEA.csv")
)

### Female donor A2
days <- c("D15", "D90")
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("day") %in% days &
  scd$getfeature("sex") %in% "F" &
  scd$getfeature("antigen") %in% "A2"


table(scd$select(b_cells = b_cells)$getfeature("day"),
      scd$select(b_cells = b_cells)$getfeature("DEA_cell_type"))

table(scd$select(b_cells = b_cells)$getfeature("day"),
      scd$select(b_cells = b_cells)$getfeature("batch"))

system(paste0("mkdir -p results/days/"))
batch_day_A2_DEA <- DEA(
  scd = scd,
  formula_null = "y ~ batch + DEA_cell_type",
  formula_full = "y ~ batch + DEA_cell_type + day",
  b_cells = b_cells,
  cpus = 11,
  v = T,
  folder_name = paste0("results/days/batch_DEA_cell_type_day_A2_F_DEA")
)
save(
  batch_day_A2_DEA,
  file = paste0("results/days/batch_DEA_cell_type_day_A2_F_DEA.Rdata")
)
system("~/scripts/sms.sh \"DEA done\"")
load(file = paste0("results/days/batch_DEA_cell_type_day_A2_F_DEA.Rdata"))
print(table(is.na(batch_day_A2_DEA$padj)))
print(table(batch_day_A2_DEA$padj < 0.05))
write.csv(
  batch_day_DEA,
  file = paste0("results/days/batch_DEA_cell_type_day_A2_F_DEA.csv")
)

ud_mbatch_day_A2_F_DEA <- up_down(
  folder_name = paste0("results/days/mbatch_DEA_cell_type_day_A2_F_DEA")
)
test_genes <- c("TSC22D3",
                "CXCR4",
                "NFKBIA",
                "MKI67"
                )
head(ud_mbatch_day_A2_F_DEA)
ud_mbatch_day_A2_F_DEA[test_genes, ]
write.csv(
  ud_mbatch_day_A2_F_DEA,
  file = paste0("results/days/ud_mbatch_DEA_cell_type_day_A2_F_DEA.csv")
)

### Female donor B7
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("day") %in% days &
  scd$getfeature("sex") %in% "F" &
  scd$getfeature("antigen") %in% "B7"


table(scd$select(b_cells = b_cells)$getfeature("day"),
      scd$select(b_cells = b_cells)$getfeature("DEA_cell_type"))

table(scd$select(b_cells = b_cells)$getfeature("day"),
      scd$select(b_cells = b_cells)$getfeature("batch"))

system(paste0("mkdir -p results/days/"))
batch_day_B7_DEA <- DEA(
  scd = scd,
  formula_null = "y ~ batch + DEA_cell_type",
  formula_full = "y ~ batch + DEA_cell_type + day",
  b_cells = b_cells,
  cpus = 11,
  v = T,
  folder_name = paste0("results/days/batch_DEA_cell_type_day_B7_F_DEA")
)
save(
  batch_day_B7_DEA,
  file = paste0("results/days/batch_DEA_cell_type_day_B7_F_DEA.Rdata")
)
system("~/scripts/sms.sh \"DEA done\"")
load(file = paste0("results/days/batch_DEA_cell_type_day_B7_F_DEA.Rdata"))
print(table(is.na(batch_day_B7_DEA$padj)))
print(table(batch_day_B7_DEA$padj < 0.05))
write.csv(
  batch_day_DEA,
  file = paste0("results/days/batch_DEA_cell_type_day_B7_F_DEA.csv")
)

ud_mbatch_day_B7_F_DEA <- up_down(
  folder_name = paste0("results/days/batch_DEA_cell_type_day_B7_F_DEA")
)
test_genes <- c("TSC22D3",
                "CXCR4",
                "NFKBIA",
                "MKI67"
                )
head(ud_mbatch_day_B7_F_DEA)
ud_mbatch_day_B7_F_DEA[test_genes, ]
write.csv(
  ud_mbatch_day_B7_F_DEA,
  file = paste0("results/days/ud_mbatch_DEA_cell_type_day_B7_F_DEA.csv")
)

