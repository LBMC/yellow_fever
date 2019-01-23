install.packages("vegan")

rm(list=ls())
setwd("~/projects/mold/yellow_fever/")
devtools::load_all("pkg/", reset = T)
load("results/QC/counts_QC.Rdata")

# check TCR correspondance
tcr_data <- c()
for (csv in list.files("data/")[grepl("TCR_.*csv", list.files("data/"))]) {
  message(csv)
  x <- read.csv(paste0("data/", csv), header = TRUE, sep = ",")
  x <- x[!(x$Clone.Number %in% "Clone Number"),]
  x <- x[!is.na(x$Timepoint),]
  tcr_data <- rbind(tcr_data, x)
}

tcr_data$TCR.Alpha.Chain[tcr_data$TCR.Alpha.Chain == 'none'] <- NA
tcr_data$TCR.Beta.Chain[tcr_data$TCR.Beta.Chain == 'none'] <- NA
tcr_data$TCR.Alpha.Chain.n2[tcr_data$TCR.Alpha.Chain.n2 == 'none'] <- NA

tcr_data$TCR.Alpha.Chain <- as.factor(tcr_data$TCR.Alpha.Chain)
tcr_data$TCR.Beta.Chain <- as.factor(tcr_data$TCR.Beta.Chain)
tcr_data$TCR.Alpha.Chain.n2 <- as.factor(tcr_data$TCR.Alpha.Chain.n2)
tcr_data$donor <- as.factor(as.vector(tcr_data$Donor_ID))
tcr_data$clonality <- as.factor(tcr_data$Clone.Number)
tcr_data$clonality[tcr_data$clonality %in% c("", 0)] <- NA
tcr_data$antigen <- as.factor(as.vector(tcr_data$Epitope))
tcr_data$Timepoint <- as.factor(tcr_data$Timepoint)
tcr_data$day <- factor(tcr_data$Timepoint,
                             levels = c(10, 15, 30, 90, 136, 148, 605, 720))
tcr_data$TCRA1 <- as.numeric(tcr_data$TCR.Alpha.Chain)
tcr_data$TCRB1 <- as.numeric(tcr_data$TCR.Beta.Chain)
tcr_data$TCRA2 <- as.numeric(tcr_data$TCR.Alpha.Chain.n2)

# we assign a clone to cells without clones

r_select <- which(is.na(tcr_data$clonality))
tcr_data$clonality <- as.vector(tcr_data$clonality)
tcr_data$clonality[r_select] <- 10000+seq_len(length(r_select))
tcr_data$clonality <- as.factor(tcr_data$clonality)

infos <- tcr_data[,which(colnames(tcr_data) %in% c("donor", "day", "antigen", "clonality"))]
summary(infos)
table(infos$donor, infos$day, infos$antigen)

system("mkdir -p results/survival/")
save(infos, file="results/survival/clone_survival.Rdata")

load(file="results/survival/clone_survival.Rdata")

# we create clone that are present at later time-point
dig_clone <- function(infos, donor, antigen){
  infos$donor <- as.vector(infos$donor)
  infos$day <- as.vector(infos$day)
  infos$antigen <- as.vector(infos$antigen)
  infos$clonality <- as.vector(infos$clonality)
  r_select <- which(infos$antigen %in% antigen & infos$donor %in% donor)
  infos_tmp <- infos[r_select,]
  infos_tmp$day <- as.numeric(infos_tmp$day)
  infos_tmp <- infos_tmp[order(infos_tmp$day),]
  for(day in unique(infos_tmp$day)){
    r_select_day <- infos_tmp$day == day
    r_select_day_late <- infos_tmp$day > day
    clones <- levels(as.factor(
                     as.vector(infos_tmp$clonality[r_select_day])))
    clones_late <- levels(as.factor(
                          as.vector(infos_tmp$clonality[r_select_day_late])))
    for(clone in clones_late){
      if(!(clone %in% clones)){
        infos <- rbind(infos,
                       c(donor, clone, antigen, day))
      }
    }
  }
  infos$donor <- as.factor(infos$donor)
  infos$day <- as.factor(infos$day)
  infos$antigen <- as.factor(infos$antigen)
  infos$clonality <- as.numeric(infos$clonality)
  return(infos)
}

for(antigen in levels(infos$antigen)){
    for(donor in levels(infos$donor)){
      infos <- dig_clone(infos, donor, antigen)
  }
}
r_select <- infos$donor %in% "YFV2003" & infos$antigen %in% "A2"
infos <- infos[!r_select,]

tmp <- as.data.frame(infos)
tmp$day <- as.numeric(as.vector(tmp$day))
tmp <- tmp[order(tmp$day), ]
tmp <- tmp[order(tmp$donor, tmp$antigen, tmp$day), ]
tmp$sampling <- 0
tmp$sampling[tmp$day %in% 10] <- 0.001 * 1/50
tmp$sampling[tmp$day %in% c(15, 30)] <- 0.0025 * 1/50
tmp$sampling[tmp$day >= 90] <- 0.0005 * 1/50
tmp$time <- "early"
tmp$time[tmp$day %in% 10] <- "early"
tmp$time[tmp$day %in% c(15, 30)] <- "peak"
tmp$time[tmp$day %in% c(90, 136, 148)] <- "late"
tmp$time[tmp$day %in% c(605, 720)] <- "very_late"
tmp$time <- factor(tmp$time, c("early", "peak", "late", "very_late"))
tmp$type <- paste(tmp$donor, tmp$antigen, tmp$day, sep = "_")
type_levels <- unique(tmp$type)
tmp$type <- as.factor(tmp$type)
tmp$type <- factor(tmp$type, levels = type_levels)
tmp <- table(tmp$type, tmp$clonality)
rownames(tmp)

simpson_infos <- simpson_infos[-1, ]
simpson_infos
table(simpson_infos$donor, simpson_infos$day, simpson_infos$antigen)
simpson_infos$donor <- factorize(simpson_infos$donor)
simpson_infos$antigen <- factorize(simpson_infos$antigen)
simpson_infos$day <- as.numeric(simpson_infos$day)
simpson_infos$simpson <- as.numeric(simpson_infos$simpson)
summary(simpson_infos)
save(simpson_infos, simpson_boots,
     file="results/survival/simpson_sampling.Rdata")

data <- simpson_infos[!is.na(simpson_infos$simpson), ]
head(data)
ggplot(data=data, aes(x=day, y=simpson, color=donor, group=donor)) +
  geom_line() +
  theme_bw() +
  facet_wrap(~antigen, scale="free") +
  geom_point(size = 1) +
  labs(x = "time",
       y =  "Simpson index",
       title = "Clone diversity through time")


library("vegan")

clone_diversity <- data.frame(
  shannon = diversity(tmp, index = "shannon"),
  simpson = diversity(tmp, index = "simpson"),
  invsimpson = diversity(tmp, index = "invsimpson")
)
H <- diversity(tmp, index = "shannon")
clone_diversity$evenness <- H/log(specnumber(tmp))
clone_diversity$specnumber <- specnumber(tmp)
clone_diversity$fisher.alpha <- fisher.alpha(tmp)
clone_diversity$rarification <- rarefy(tmp, min(rowSums(tmp)))

preston_mode <- c()
preston_width <- c()
preston_S0 <- c()
for (i in 1:nrow(tmp)) {
  prest <- prestondistr(tmp[i, ])
  preston_mode <- c(preston_mode, prest$coefficients["mode"])
  preston_width <- c(preston_width, prest$coefficients["width"])
  preston_S0 <- c(preston_S0, prest$coefficients["S0"])
}
clone_diversity$preston_mode <- preston_mode
clone_diversity$preston_width <- preston_width
clone_diversity$preston_S0 <- preston_S0

Mandelbrot_c <- c()
Mandelbrot_gamma <- c()
Mandelbrot_beta <- c()
for (i in 1:nrow(tmp)) {
  rad <- radfit(tmp[i, ])
  Mandelbrot_c <- c(Mandelbrot_c,
                    rad$models$Mandelbrot$coefficients["c"])
  Mandelbrot_gamma <- c(Mandelbrot_gamma,
                        rad$models$Mandelbrot$coefficients["gamma"])
  Mandelbrot_beta <- c(Mandelbrot_beta,
                       rad$models$Mandelbrot$coefficients["beta"])
}
clone_diversity$Mandelbrot_c <- Mandelbrot_c
clone_diversity$Mandelbrot_gamma <- Mandelbrot_gamma
clone_diversity$Mandelbrot_beta <- Mandelbrot_beta

R_est <- t(data.frame(
                    S.obs = NA,
                    S.chao1 = NA,
                    se.chao1 = NA,
                    S.ACE = NA,
                    se.ACE = NA
                    ))
for (i in 1:nrow(tmp)) {
  R_est <- cbind(R_est, estimateR(tmp[i, ]))
}
R_est <- R_est[, -1]
colnames(R_est) <- row.names(tmp)
clone_diversity <- cbind(clone_diversity, t(R_est))

save(clone_diversity, file="results/survival/clone_diversity.Rdata")
write.csv(clone_diversity, file="results/survival/clone_diversity.csv")

pdf(file = "results/survival/clone_diversity_fisher_scaled.pdf", height = 10, width = 10)
par(mfrow=c(4,4))
for (i in 1:nrow(tmp)) {
  str(fisherfit(tmp[i, ]))
  plot(fisherfit(tmp[i, ]), main = rownames(tmp)[i])
}
dev.off()

pdf(file = "results/survival/clone_diversity_fisher.pdf", height = 10, width = 10)
par(mfrow=c(4,4))
x <- strsplit(rownames(tmp)[1], "_")[[1]]
prev_donor <- x[1]
prev_antigen <- x[2]
prev_specnumber <- specnumber(tmp[1, ])
for (i in 1:nrow(tmp)) {
  x <- strsplit(rownames(tmp)[i], "_")[[1]]
  donor <- x[1]
  antigen <- x[2]
  if (donor != prev_donor | antigen != prev_antigen) {
    prev_specnumber <- specnumber(tmp[i, ])
  }
  x <- fisherfit(tmp[i, ])
  freq <- as.numeric(names(x$fisher))
  plot(freq, x$fisher, main = rownames(tmp)[i], ylab = "Species",
       xlab = "Frequency",
       ylim = c(0, prev_specnumber), xlim = c(0.5, max(freq) + 0.5),
       type = "n")
  rect(freq - 0.5, 0, freq + 0.5, x$fisher, col = "skyblue")
  alpha <- x$estimate
  k <- x$nuisance
  curve(alpha * k^x/x, 1, max(freq), col = "red", lwd = 2,
      add = TRUE)
  prev_donor <- donor
  prev_antigen <- antigen
}
dev.off()

diversity_s <- function(infos, size = 50, samples = 10, ncore = 4,
                        sampling = TRUE, sample_size = 0){
  require("vegan")
  if(sampling){
    results <- as.list(1:samples)
    results <- lapply(results, FUN = function(x,
                                            infos,
                                            size){
        rarefy(round(sample(infos, size = size, replace = TRUE)),
               round(sample_size))
      },
      infos,
      size)
    return(mean(unlist(results)))
  }
  return(rarefy(infos, sample_size))
}

boot_simpson <- function(by, boot = 1000, ncore = 12,
                         sample_size = 1, size = 50, sampling = TRUE){
  require("vegan", "parallel")
  results <- as.list(1:boot)
  results <- mclapply(results,
                      FUN = function(x,
                                     by,
                                     sample_size,
                                     sampling,
                                     size){
      return(diversity_s(by,
                         sampling = sampling,
                         sample_size = sample_size,
                         size = size
                         )
      )
    },
    by,
    sample_size,
    sampling,
    size,
    mc.cores=ncore)
  return(unlist(results))
}

simpson_infos <- data.frame(id= NA, donor=NA, antigen=NA, day=NA,
                            rarefication=NA)
simpson_boots <- data.frame()
id <- 1
for(antigen in levels(infos$antigen)){
  for(donor in levels(infos$donor)){
    r_select <- which(infos$antigen %in% antigen & infos$donor %in% donor)
    sample_size <- as.vector(table(factorize(infos$day[r_select])))
    if (length(sample_size) > 1){
      print(antigen)
      print(donor)
      print(table(factorize(infos$day[r_select])))
      sample_size <- min(sample_size)
      for(day in levels(infos$day)){
        r_select <- which(infos$antigen %in% antigen &
                          infos$day %in% day & infos$donor %in% donor)
        if(length(r_select) > 0){
          infos_tmp <- (infos$clonality[r_select])
          x <- diversity_s(infos_tmp, size = length(r_select),
                           sample_size = round(sample_size))
          simpson_infos <- rbind(simpson_infos,
                                c(id, donor, antigen, day, x))
          simpson_boots <- rbind(simpson_boots,
                                boot_simpson(infos_tmp, boot = 1000,
                                             sample_size = sample_size,
                                             size = length(r_select)
                                             )
                                )
          id <- id + 1
        }
      }
    }
  }
}

simpson_infos <- simpson_infos[-1,]
simpson_infos$donor <- factorize(simpson_infos$donor)
simpson_infos$antigen <- factorize(simpson_infos$antigen)
simpson_infos$day <- as.numeric(simpson_infos$day)
simpson_infos$rarefication <- as.numeric(simpson_infos$rarefication)

simpson_infos$rarefication_min <- apply(simpson_boots, 1, FUN = function(x){quantile(as.numeric(x), 0.05)})
simpson_infos$rarefication_max <- apply(simpson_boots, 1, FUN = function(x){quantile(as.numeric(x), 0.95)})
simpson_infos
save(simpson_infos, simpson_boots, file="results/survival/rarefication.Rdata")
write.csv(simpson_infos, file = "results/survival/rarefication.csv")

infos$clonality_p <- paste0(infos$donor, infos$clonality)

for(antigen in levels(infos$antigen)){
  message(antigen)
  for(donor in levels(infos$donor)){
    message(donor)
    r_select <- infos$donor %in% donor & infos$antigen %in% antigen
    assign(paste0("infos_", donor, "_", antigen),
           factorize(infos[r_select,]),
           envir = .GlobalEnv)
  }
}

clone_size <- function(data, rm.unique = F){
  if (rm.unique) {
    data <- data[as.numeric(as.vector(data$clonality)) < 10000,]
  }
  data <- table(factorize(data$clonality_p), as.numeric(as.vector(data$day)))
  data <- as.matrix(data)
  for(i in seq_len(ncol(data)-1)){
    data[data[,i] == 0 & data[,i+1] != 0,i] <- 1
  }
  data <- cbind(data,
                data[,ncol(data)] != 0)
  colnames(data) <- c(paste0("D",colnames(data)[1:(ncol(data)-1)]), "LLC")
  data <- as.data.frame(data)
  data <- data[, -ncol(data)]
  return(data)
}
surival_prob <- function(coefficients){
  coefficients <- as.data.frame(coefficients)
  coefficients$OR <- exp(coefficients[, 1])
  coefficients$proba <- coefficients$OR / (1 + coefficients$OR)
  return(coefficients)
}

library("reshape2")
library("mice")
data <- clone_size(infos_YFV16_A2, T)
g <- summary(glm((D720!=0)~-1+D15+D90+D15:D90, data=data, family="binomial"))
g <- surival_prob(g$coefficients)
g
data <- clone_size(infos_YFV16_A2)
g <- summary(glm((D720!=0)~-1+D15+D90+D15:D90, data=data, family="binomial"))
g <- surival_prob(g$coefficients)
g
write.csv(g, file = "results/survival/survival_YFV16_A2.csv")
data <- clone_size(infos_YFV2001_A2, T)
g <- summary(glm((D605!=0)~-1+D15+D136, data=data, family="binomial"))
g <- surival_prob(g$coefficients)
g
data <- clone_size(infos_YFV2001_A2)
g <- summary(glm((D605!=0)~-1+D15+D136+D15:D136, data=data, family="binomial"))
g <- surival_prob(g$coefficients)
g
write.csv(g, file = "results/survival/survival_YFV2001_A2.csv")
data <- clone_size(infos_YFV5_A2, T)
g <- summary(glm((D90!=0)~-1+D15, data=data, family="binomial"))
g <- surival_prob(g$coefficients)
g
data <- clone_size(infos_YFV5_A2)
g <- summary(glm((D90!=0)~-1+D15, data=data, family="binomial"))
g <- surival_prob(g$coefficients)
g
write.csv(g, file = "results/survival/survival_YFV5_A2.csv")
data <- clone_size(infos_YFV2001_B7, T)
g <- summary(glm((D136!=0)~-1+D15, data=data, family="binomial"))
g <- surival_prob(g$coefficients)
g
data <- clone_size(infos_YFV2001_B7)
g <- summary(glm((D136!=0)~-1+D15, data=data, family="binomial"))
g <- surival_prob(g$coefficients)
g
write.csv(g, file = "results/survival/survival_YFV2001_B7.csv")
data <- clone_size(infos_YFV2003_B7, T)
g <- summary(glm((D148!=0)~-1+D15+D30+D90, data=data, family="binomial"))
g <- surival_prob(g$coefficients)
g
data <- clone_size(infos_YFV2003_B7)
g <- summary(glm((D148!=0)~-1+D15+D30+D90, data=data, family="binomial"))
g <- surival_prob(g$coefficients)
g
write.csv(g, file = "results/survival/survival_YFV2004_B7.csv")

# fisherplot of clone size
devtools::install_github("chrisamiller/fishplot")
library("fishplot")

fish_plot <- function(data, timepoints, title, min_size = 2) {
  frac.table <- as.matrix(data)
  frac.table <- frac.table[rowSums(frac.table) > min_size, ]
  frac.table <- apply(frac.table, 2, function(x){
    x <- x / sum(x) * 100
  })
  frac.table <- frac.table[order(frac.table[, 2], frac.table[, 3], frac.table[, 1]), ]
  parents <- rep(0, nrow(frac.table))
  fish <- createFishObject(frac.table, parents, timepoints = timepoints)
  fish <- layoutClones(fish)
  fish <- setCol(fish, rainbow(nrow(frac.table)))
  pdf(file = paste0("results/survival/fish_plot", title, ".pdf"),
      height = 10, width = 10)
  fishPlot(fish, shape = "spline", title.btm = title,
            vlines = timepoints,
            vlab = paste("day", timepoints))
  dev.off()
}

data <- clone_size(infos_YFV16_A2)
fish_plot(data, timepoints=c(15,90,720), "YF16_A2", min_size = 2)
data <- clone_size(infos_YFV2001_A2)
fish_plot(data, timepoints=c(15,136,605), "YF2001_A2", min_size = 2)

