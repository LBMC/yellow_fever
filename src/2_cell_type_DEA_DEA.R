setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)
load("results/cell_type/cells_counts_QC_DEA_cell_type.Rdata")

for (day in c("D15", "D136", "D593")) {
  system(
    paste0("mkdir -p results/cell_type/mbatch_", day, "_surface_cell_type_DEA")
  )
  b_cells <- scd$getfeature("QC_good") %in% T &
    !is.na(scd$getfeature("DEA_cell_type")) &
    scd$getfeature("day") %in% day
  mbatch_surface_cell_type_DEA <- DEA(
    scd = scd,
    formula_null = "y ~ (1|batch)",
    formula_full = "y ~ (1|batch) + DEA_cell_type",
    b_cells = b_cells,
    cpus = 16,
    v = F,
    folder_name = paste0("results/cell_type/mbatch_", day, "DEA_cell_type_DEA")
  )
  save(
    mbatch_surface_cell_type_DEA,
    file = paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA.Rdata")
  )
  system("~/scripts/sms.sh \"DEA done\"")
  print(day)
  print(table(is.na(mbatch_surface_cell_type_DEA$padj)))
  print(table(mbatch_surface_cell_type_DEA$padj < 0.05))
}

################################################################################
## DEA for the F donor

load("results/cell_type/cells_counts_QC_DEA_cell_type_F.Rdata")

for (day in c("D15", "D90")) {
  system(
    paste0("mkdir -p results/cell_type/mbatch_", day, "_surface_cell_type_DEA_F")
  )
  b_cells <- scd$getfeature("QC_good") %in% T &
    scd$getfeature("sex") %in% "F" &
    !is.na(scd$getfeature("DEA_cell_type")) &
    scd$getfeature("day") %in% day
  mbatch_surface_cell_type_DEA <- DEA(
    scd = scd,
    formula_null = "y ~ (1|batch)",
    formula_full = "y ~ (1|batch) + DEA_cell_type",
    b_cells = b_cells,
    cpus = 16,
    v = F,
    folder_name = paste0("results/cell_type/mbatch_", day, "DEA_cell_type_DEA_F")
  )
  save(
    mbatch_surface_cell_type_DEA,
    file = paste0("results/cell_type/mbatch_", day, "_DEA_cell_type_DEA_F.Rdata")
  )
  system("~/scripts/sms.sh \"DEA done\"")
  print(day)
  print(table(is.na(mbatch_surface_cell_type_DEA$padj)))
  print(table(mbatch_surface_cell_type_DEA$padj < 0.05))
}
