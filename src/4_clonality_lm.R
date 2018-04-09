setwd("~/projects/yellow_fever")
devtools::load_all("../scRNAtools/", reset = T)
load("results/cell_type/CB_counts_QC_DEA_cell_type.Rdata")

system("mkdir -p results/clonality/")

days <- c("D15", "D136", "D593")
b_cells <- scd$getfeature("QC_good") %in% T &
  !is.na(scd$getfeature("DEA_cell_type")) &
  scd$getfeature("day") %in% days &
  scd$getfeature("sex") %in% "M"


clone_names <- table(
  factorize(scd$select(b_cells = b_cells)$getfeature("clonality")),
  factorize(scd$select(b_cells = b_cells)$getfeature("day"))
)
clone_names <- clone_names >= 2
clone_names <- rownames(clone_names)[apply(clone_names, 1, function(x){
  x <- as.numeric(paste0(as.numeric(x), collapse=""))
  x %in% c(111,
           101,
           110,
           011)
  })]
data_tmp <- scd$select(b_cells = b_cells)$getfeatures[scd$select(b_cells = b_cells)$getfeature("clonality") %in% clone_names,]
data_tmp$clonality <- factorize(data_tmp$clonality)

data_tmp$D15_cell_type <- NA
for(clone in levels(data_tmp$clonality)){
  r_select <- data_tmp$clonality %in% clone
  last_day <- table(data_tmp$day[r_select])
  last_day <- ifelse(last_day[3] == 0, "D136", "D593")
  r_select <- data_tmp$clonality %in% clone & data_tmp$day %in% last_day
  D15_cell_type <- ifelse(mean(data_tmp$pDEA_cell_type[r_select]) >= 0.5,
                          "Memory", "Effector")
  r_select <- data_tmp$clonality %in% clone
  data_tmp$D15_cell_type[r_select] <- D15_cell_type
}

data_tmp$D15_cell_type[data_tmp$clonality==22] <- "Effector"
data_tmp$D15_cell_type <- factorize(data_tmp$D15_cell_type)
data_tmp$clonality <- factor(
  data_tmp$clonality,
  levels=levels(factorize(data_tmp$clonality))[
    order(by(data_tmp$pDEA_cell_type, data_tmp$clonality, mean))
  ]
)
data_tmp$day <- factor(data_tmp$day, levels = days)



g <- ggplot() +
geom_smooth(method="lm", se=F, data=data_tmp, aes(x=day, y=pDEA_cell_type, group=clonality, color=D15_cell_type)) +
geom_jitter(height=0, width=0.25, size=2, data=data_tmp, aes(x=day, y=pDEA_cell_type)) +
scale_color_manual(values=cell_type_palette(levels(factorize(data_tmp$D15_cell_type)))) +
theme_bw() +
labs(x = "day",
     y =  "probability of being an Memory cell",
     title = "Evolution of clone cell-type through time",
     colour="cell-type at D15")
print(g)
ggsave(file = "results/clonality/lm_pDEA_clone_M.pdf")

g <- ggplot() +
geom_smooth(method="lm", se=F, data=data_tmp, aes(x=day, y=pDEA_cell_type, group=clonality, color=D15_cell_type)) +
scale_color_manual(values=cell_type_palette(levels(factorize(data_tmp$D15_cell_type)))) +
geom_jitter(height=0, width=0.25, size=2, data=data_tmp, aes(x=day, y=pDEA_cell_type, fill=pDEA_cell_type), shape=21) +
scale_fill_gradient2(low="blue", mid="gray", high="red", midpoint=0.5) +
theme_bw() +
labs(x = "day",
     y =  "probability of being an Memory cell",
     title = "Evolution of clone cell-type through time",
     colour="cell-type at D15")
print(g)
ggsave(file = "results/clonality/lm_colored_pDEA_clone_M.pdf")

g <- ggplot() +
geom_smooth(method="lm", se=F, data=data_tmp, aes(x=day, y=pDEA_cell_type, group=clonality, color=D15_cell_type)) +
geom_jitter(height=0, width=0.25, size=2, data=data_tmp, aes(x=day, y=pDEA_cell_type)) +
scale_color_manual(values=cell_type_palette(levels(factorize(data_tmp$D15_cell_type)))) +
theme_bw() +
facet_wrap(~clonality) +
labs(x = "day",
     y =  "probability of being an Memory cell",
     title = "Evolution of clone cell-type through time",
     colour="cell-type at D15") +
theme(
    panel.margin = unit(0.5, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    strip.text.y = element_text(angle = 180),
    strip.background = element_blank()
  )
print(g)
ggsave(file = "results/clonality/lm_per_clone_pDEA_clone_M.pdf")
