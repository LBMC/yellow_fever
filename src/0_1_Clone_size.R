rm(list=ls())
setwd("~/projects/mold/yellow_fever")
library("tidyverse")
library("readxl")
library("vegan")
library("broom")
library("pbmcapply")
theme_set(theme_classic())
pal = "Set1"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}

# read frequency sampling data
freq_resp <- read_xlsx("data/2019_08_22_Frequency_Responses_Fig1_jmup.xlsx") %>%
  as_tibble() %>%
  gather(key = "donor", value = "percent", -"day") %>%
  separate(donor, into = c("donor", "antigen"), sep = "_") %>%
  drop_na() %>%
  mutate(day = strsplit(day, "D") %>% unlist() %>% as.numeric() %>% .[!is.na(.)])
freq_resp
# read clone data
clone <- read_xlsx("data/2019_08_22_Table_S1_Aug20.xlsx", sheet = 1) %>%
  as_tibble() %>%
  rename(donor="Donor_ID",
         day="Timepoint",
         antigen="Epitope",
         clone="clone_id"
  ) %>%
  mutate(day = strsplit(day, "D") %>% unlist() %>% as.numeric() %>% .[!is.na(.)]) %>%
  mutate(day = replace(day, day %in% "605", "593"),
         day = as.numeric(day)) %>%
  select(donor, day, antigen, clone)
clone
for (i in 2:6) {
  print(i)
  clone <- read_xlsx("data/2019_08_22_Table_S1_Aug20.xlsx", sheet = i) %>%
    as_tibble() %>%
    rename(donor = "Donor_ID",
          day = "Timepoint",
          antigen = "Epitope",
          clone = "clone_id"
    ) %>%
    mutate(day = replace(day, day %in% "605", "593"),
           day = as.numeric(day)) %>%
    select(donor, day, antigen, clone) %>%
    bind_rows(clone)
}
# generate size 1 clone
size_1_number <- clone %>% filter(clone %in% 0) %>% nrow()
highest_clone_number <- clone %>% drop_na() %>% pull(clone) %>% max()
size_1_name <- ( highest_clone_number + 1 ):(highest_clone_number + 1 + size_1_number)
clone <- clone %>%
  left_join(freq_resp) %>%
  drop_na() %>%
  mutate(donor = as.factor(donor),
         day = as.factor(day),
         antigen = as.factor(antigen),
         clone = replace(clone, clone == 0, size_1_name),
         clone = as.factor(clone)
  )
clone
donors <- clone %>% pull(donor) %>% levels()
days <-  clone %>% pull(day) %>% levels()
antigens <- clone %>% pull(antigen) %>% levels()
clone <- clone %>% filter(donor %in% "YFV5") %>%
  filter(!(clone == 1)) %>%
  mutate(donor = paste(donor, "nobig")) %>%
  bind_rows(clone %>% mutate(donor = as.vector(donor))) %>%
  mutate(donor = as.factor(donor))
clone
# reconstruct missing clone
# we create clone that are present at later time-point
dig_clone <- function(infos, donor, antigen){
  infos$donor <- as.vector(infos$donor)
  infos$day <- as.vector(infos$day)
  infos$antigen <- as.vector(infos$antigen)
  infos$clone <- as.vector(infos$clone)
  r_select <- which(infos$antigen %in% antigen & infos$donor %in% donor)
  infos_tmp <- infos[r_select,]
  infos_tmp$day <- as.numeric(infos_tmp$day)
  infos_tmp <- infos_tmp[order(infos_tmp$day),]
  for(day in unique(infos_tmp$day)){
    r_select_day <- infos_tmp$day == day
    r_select_day_late <- infos_tmp$day > day
    clones <- levels(as.factor(
                     as.vector(infos_tmp$clone[r_select_day])))
    clones_late <- levels(as.factor(
                          as.vector(infos_tmp$clone[r_select_day_late])))
    for(clone in clones_late){
      if(!(clone %in% clones)){
        infos <- rbind(infos,
                       c(donor, day, antigen, clone, infos_tmp$percent[r_select_day]))
      }
    }
  }
  infos$donor <- as.factor(infos$donor)
  infos$day <- as.factor(infos$day)
  infos$antigen <- as.factor(infos$antigen)
  infos$clone <- as.numeric(infos$clone)
  return(infos)
}
salt <- 0
for(antigen in levels(clone$antigen)) {
  for(donor in levels(clone$donor)) {
    clone <- dig_clone(clone, donor, antigen)
    salt <- salt + nrow(clone)
  }
}
clone <- clone %>%
  right_join(clone %>% group_by(donor, clone, day) %>% count()) %>%
  mutate(n = replace(n, clone == 0, 1))
alpha_f <- tibble(donor = NA, day = NA, antigen = NA,
         alpha = NA, alpha_prec = NA)
donors <- clone %>% pull(donor) %>% levels()
for (i in 1:length(donors)) {
  for (j in 1:length(days)) {
    for (k in 1:length(antigens)) {
      if (clone %>%
        filter(donor %in% donors[i],
               day %in% days[j],
               antigen %in% antigens[k]
              ) %>% nrow() > 0 ) {
        tmp <- clone %>%
          filter(donor %in% donors[i],
                day %in% days[j],
                antigen %in% antigens[k]
                ) %>%
          group_by(clone) %>%
          count() %>%
          pull(n)
        alpha_f_fit <- tmp %>% fisherfit()
        alpha_f <- tibble(donor = donors[i], day = days[j], antigen = antigens[k],
                alpha = alpha_f_fit$estimate,
                alpha_prec = alpha_f_fit$estim.prec,
                shannon = diversity(tmp, index = "shannon"),
                simpson = diversity(tmp, index = "simpson"),
                invsimpson = diversity(tmp, index = "invsimpson"),
                evenness = diversity(tmp, index = "shannon")/log(specnumber(tmp)),
                specnumber = specnumber(tmp),
                fisher.alpha = fisher.alpha(tmp),
                rarification = rarefy(tmp, min(sum(tmp)))
                ) %>%
          bind_rows(alpha_f)
        }
    }
  }
}

alpha_f %>% select(-alpha_prec) %>%
  drop_na() %>%
  write.csv(file="results/survival/clone_diversity.csv")

alpha_f %>% select(-alpha_prec) %>%
  drop_na() %>%
  mutate(day = as.numeric(day),
         response = paste(donor, antigen)) %>%
  lm(alpha ~ -1 + day + donor, data = .) %>%
  summary()

# test for alpha decreasse through time accounting for response effect
library(lme4)
library(pbkrtest)
alpha_regression <- alpha_f %>% select(-alpha_prec) %>%
  drop_na() %>%
  mutate(day = as.numeric(day),
         response = paste(donor, antigen)) %>%
  lmer(log(alpha) ~ day + (1|response), data = .)
alpha_regression %>% summary() %>% coef() %>% as_tibble() %>%
  rename("estimate" = "Estimate",
         "st.error" = "Std. Error",
         "t.value" = "t value") %>%
  mutate(df = get_ddf_Lb(alpha_regression, fixef(alpha_regression)),
         pvalue = (1 - pt(t.value, df, lower.tail = F)))

clone <- alpha_f %>%
  drop_na() %>%
  mutate(day = as.numeric(as.vector(day))) %>%
  right_join(clone %>%
    mutate(donor = as.vector(donor),
           day = as.numeric(as.vector(day)),
           antigen = as.vector(antigen)
    )
  ) %>%
  mutate(donor = as.factor(donor),
         day = as.vector(day),
         antigen = as.factor(antigen),
         clone = as.factor(clone)
  ) %>%
  distinct()

clone %>% filter(donor %in% "YFV5" & day == 15) %>%
  group_by(clone) %>%
  count() %>%
  arrange(desc(nn))
clone %>% filter(donor %in% "YFV5" & day == 15) %>%
  arrange(desc(n))

clone %>% mutate(Donor = paste(donor, antigen, day, sep = "_")) %>%
  select(Donor, alpha) %>%
  distinct()

clone %>%
  mutate(percent = as.numeric(as.vector(percent))) %>%
  ggplot() +
  geom_histogram(aes(x = n, group = donor, fill = antigen)) +
  facet_wrap(~donor + antigen + day, scale = "free") +
  scale_y_log10() +
  labs(y = "clone size")
ggsave("results/survival/clone_size_histogram_scaled.pdf", dpi = 400)

par(mfrow=c(1,1))
clone %>%
  mutate(percent = as.numeric(as.vector(percent))) %>%
  gather(percent, alpha, key="measure", value = "value") %>%
  mutate(value = as.numeric(as.vector(value)),
         measure = replace(measure, measure %in% "alpha", "Fisher's Alpha"),
         measure = replace(measure, measure %in% "percent", "percentage of CD8+ T cells")) %>%
  ggplot() +
  geom_point(aes(x = day, y = value, group = donor, color = antigen, shape = donor), size = 4) +
  geom_line(aes(x = day, y = value, group = paste0(donor, antigen), color = antigen)) +
  facet_wrap(~measure, scale = "free_y") +
  labs(y = "")
ggsave("results/survival/fisher_vs_sampling_vs_time.pdf")

clone %>%
  mutate(percent = as.numeric(as.vector(percent))) %>%
  ggplot() +
  geom_point(aes(x = alpha, y = percent, group = donor, color = antigen, shape = donor), size = 4) +
  geom_line(aes(x = alpha, y = percent, group = paste0(donor, antigen), color = antigen)) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Fisher's Alpha",
       y = "percentage of CD8+ T cells"
  )
ggsave("results/survival/fisher_vs_sampling_log10.pdf")

clone %>%
  mutate(percent = as.numeric(as.vector(percent))) %>%
  ggplot() +
  geom_point(aes(x = alpha, y = percent, group = donor, color = antigen, shape = donor), size = 4) +
  geom_line(aes(x = alpha, y = percent, group = paste0(donor, antigen), color = antigen)) +
  labs(x = "Fisher's Alpha",
       y = "percentage of CD8+ T cells"
  )
ggsave("results/survival/fisher_vs_sampling.pdf")

clone %>%
  mutate(percent = as.numeric(as.vector(percent))) %>%
  filter(day %in% 15) %>%
  ggplot() +
  geom_point(aes(x = alpha, y = percent, group = donor, color = antigen, shape = donor), size = 4) +
  geom_line(aes(x = alpha, y = percent, group = paste0(donor, antigen), color = antigen)) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Fisher's Alpha",
       y = "percentage of CD8+ T cells"
  )
ggsave("results/survival/fisher_vs_sampling_log10_D15.pdf")


clone %>%
  select(alpha, day, percent, antigen) %>%
  distinct() %>%
  mutate(percent = as.numeric(as.vector(percent))) %>%
  lm(data = ., log10(alpha) ~ log10(percent) * antigen) %>%
  summary()

clone %>%
  select(alpha, day, percent, antigen) %>%
  distinct() %>%
  mutate(percent = as.numeric(as.vector(percent))) %>%
  lm(data = ., log10(percent) ~ log10(alpha) * antigen) %>%
  summary()

clone %>%
  select(alpha, day, percent, antigen) %>%
  filter(day %in% 15) %>%
  distinct() %>%
  mutate(percent = as.numeric(as.vector(percent))) %>%
  lm(data = ., alpha ~ percent )%>%
  summary()

clone %>%
  select(alpha, day, percent, antigen) %>%
  filter(day %in% 15) %>%
  distinct() %>%
  mutate(percent = as.numeric(as.vector(percent))) %>%
  t.test(data = ., alpha ~ antigen, alt = "greater")

clone %>%
  select(alpha, day, percent, antigen) %>%
  distinct() %>%
  mutate(percent = as.numeric(as.vector(percent))) %>%
  glm(data = ., antigen ~ log10(percent) * log10(alpha), family = "binomial") %>%
  summary()

clone %>%
  select(alpha, day, percent, antigen) %>%
  filter(day %in% 15) %>%
  distinct() %>%
  mutate(percent = as.numeric(as.vector(percent))) %>%
  glm(data = ., antigen ~ log10(percent) + log10(alpha), family = "binomial") %>%
  summary()

par(mfrow=c(2,2))
clone %>%
  filter(day %in% c("D15")) %>%
  glm(data = ., antigen ~ n * percent, family = "binomial") %>%
  plot()

clone %>%
  filter(donor=="Donor A" & antigen=="A2") %>%
  select(donor, day, n) %>%
  glm(data = ., n~day, family = "binomial") %>%
  summary()

clone_size <- function(clone, select_donor, select_antigen){
  data <- clone %>%
    filter(donor==select_donor & antigen==select_antigen) %>%
    select(day, clone, n)
  result <- tibble(clone = data %>% pull(clone) %>% unique())
  for(nday in data %>% arrange(day) %>% pull(day) %>% unique()){
    cday <- paste0("D", nday)
    result <- result %>% add_column(!!(cday) := 0)
    for (clone_number in result %>% pull(clone)) {
      n <- data %>%
        filter(clone==clone_number & day==as.numeric(nday)) %>%
        pull(n)
      if (length(n) != 0) {
        result[result$clone == clone_number, cday] <- n
        prev_col <- ( result %>% ncol() )-1
        if(result[result$clone == clone_number, prev_col ] == 0) {
          result[result$clone == clone_number, prev_col ] <- 1
        }
      }
    }
  }
  return(result)
}
clone_size(clone=clone, select_donor="Donor A", select_antigen="A2") %>%
  glm(data = ., (D593!=0)~-1+D15*D136, family = "binomial") %>%
  tidy() %>%
  rename(`z value`="statistic",
         `Pr(|>z|)`="p.value"
  ) %>%
  mutate(OR = exp(estimate),
         proba = OR / (1 + OR)
  ) %>%
  write.csv(file="results/survival/DonorA_A2_survival_analysis.csv")


library("fishplot")
fish_plot <- function(data, timepoints, title, min_size = function(x){any(x > 1)}) {
  colnames(data) <- as.numeric(substring(colnames(data), 2))
  frac.table <- as.matrix(data)
  frac.table <- frac.table[apply(frac.table, 1, min_size ), ]
  print(head(frac.table))
  frac.table <- apply(frac.table, 2, function(x){
    x <- x / sum(x) * 100
  })
  print(head(frac.table))
  print(colSums(frac.table))
  table_order <- hclust(dist(frac.table))$order
  # frac.table <- frac.table[order(frac.table[, 3], frac.table[, 2], frac.table[, 1]), ]
  frac.table <-
    frac.table[order(frac.table[, 3] > 0, frac.table[, 3] - frac.table[, 1]),]
  parents <- rep(0, nrow(frac.table))
  fish <-
    createFishObject(frac.table, parents, timepoints = timepoints)
  fish <- layoutClones(fish)
  fish <- setCol(fish, rainbow(nrow(frac.table)))
  pdf(file = paste0("results/survival/fish_plot", title, "_3-2-1.pdf"),
      height = 10, width = 10)
  fishPlot(
    fish,
    shape = "spline",
    title.btm = title,
    vlines = timepoints,
    vlab = paste("day", timepoints)
  )
  dev.off()
}

data <- clone_size(clone=clone, select_donor="Donor A", select_antigen="A2") %>%
  select(D15, D136, D593)
fish_plot(data, timepoints=c(15,136,593), "Donor A A2", min_size = function(x){any(x > 3)})

data <- clone_size(clone=clone, select_donor="Donor D", select_antigen="A2") %>%
  select(D15, D90, D720)
fish_plot(data, timepoints=c(15,90,720), "Donor D A2", min_size = function(x){any(x > 3)})


library(lme4)

# fisher alpha subsampling experiment
data <- clone %>% 
  mutate(day = fct_reorder(day, as.numeric(as.vector(day)))) %>% 
  group_by(donor, day, antigen) %>% 
  select(-percent) %>% 
  nest() %>% 
  mutate(alpha = lapply(data, function(x){
      n_sample <- 100
      tibble(
        sampling = seq(from = 0.1, to = 1, by = 0.1) %>% rep(n_sample),
        sample = rep(1:n_sample, each = 10)
        ) %>%
        mutate(
          alpha = lapply(sampling, function(y, x){
            x %>%
            pull(n) %>%
            sample(round(length(.) * y)) %>%
            fisherfit(.) %>%
            .$estimate
          }, x = x) %>% unlist()
        )
    })
  ) %>% 
  unnest(alpha) %>% 
  mutate(
    sample = as.factor(sample),
    percent = sampling * 100
  ) %>% 
  filter(alpha <= 2.0e9) %>% 
  ggplot() +
  geom_line(aes(
    x = alpha,
    y = percent,
    color = day,
    group = str_c(sample, day)
    ),
    alpha = 0.1
  ) +
  geom_smooth(aes(
    x = alpha,
    y = percent,
    color = day,
    group = day
    ),
    se = F
  ) +
  facet_wrap(~donor + antigen, scales = "free")
  

# abs number of cells
clone %>%
  mutate(day = fct_reorder(day, as.numeric(as.vector(day)))) %>% 
  select(-percent) %>% 
  group_by(donor, antigen, day) %>% 
  nest() %>% 
  mutate(alpha = pbmcapply::pbmclapply(data, function(data){
      n_sample <- 10
      tibble(
        sampling = seq(
          from = 20,
          to = 2000,
          length.out = 20) %>% rep(n_sample),
        sample = rep(1:n_sample, each = 20)
        ) %>%
        mutate(
          alpha = lapply(sampling, function(sampling, data){
            data %>%
            pull(n) %>%
            sample(round(sampling), replace = T) %>%
            fisherfit(.) %>%
            .$estimate
          }, data = data) %>% unlist()
        )
    },
    mc.cores = 10,
    ignore.interactive = T )) %>% 
  unnest(alpha) %>% 
  mutate(
    sample = as.factor(sample),
    n_cell = sampling
  ) %>% 
  filter(alpha <= 2.0e9) %>% 
  ggplot() +
  geom_line(aes(
    x = alpha,
    y = n_cell,
    color = day,
    group = str_c(sample, day)
    ),
    alpha = 0.1
  ) +
  geom_smooth(aes(
    x = alpha,
    y = n_cell,
    color = day,
    group = day
    ),
    se = F
  ) +
  facet_wrap(~donor + antigen, scales = "free")
  

# clone number subsampling experiment
data <- clone %>% 
  mutate(day = fct_reorder(day, as.numeric(as.vector(day)))) %>% 
  group_by(donor, day, antigen) %>% 
  select(-percent) %>% 
  nest() %>% 
  mutate(alpha = lapply(data, function(x){
      n_sample <- 100
      tibble(
        sampling = seq(from = 0.1, to = 1, by = 0.1) %>% rep(n_sample),
        sample = rep(1:n_sample, each = 10)
        ) %>%
        mutate(
          alpha = lapply(sampling, function(y, x){
            x %>%
            pull(n) %>%
            sample(round(length(.) * y)) %>%
            fisherfit(.) %>%
            .$estimate
          }, x = x) %>% unlist()
        )
    })
  ) %>% 
  unnest(alpha) %>% 
  mutate(
    sample = as.factor(sample),
    percent = sampling * 100
  )
data %>% 
  ggplot() +
  geom_line(aes(
    x = alpha,
    y = percent,
    color = day,
    group = str_c(sample, day)
    ),
    alpha = 0.1
  ) +
  geom_smooth(aes(
    x = alpha,
    y = percent,
    color = day,
    group = day
    ),
    se = F
  ) +
  facet_wrap(~donor + antigen, scales = "free")


# abs number of cells
clone
data <- clone %>% 
  mutate(day = fct_reorder(day, as.numeric(as.vector(day)))) %>% 
  group_by(donor, day, antigen) %>% 
  select(-percent) %>% 
  nest() %>% 
  mutate(detected_clone = lapply(data, function(data){
      n_sample <- 1000
      tibble(
        n_cell = seq(
          from = 100,
          to = max(500, nrow(data) + 10),
          step = 1) %>%
          rep(n_sample),
        sample = rep(
          1:n_sample,
          each = (
            seq(
              from = 100,
              to = max(500, nrow(data) + 10),
              step = 1) %>%
                length()
          )),
        day_size = nrow(data)
      ) %>%
        mutate(
          detected_clone = pbmcapply::pbmclapply(n_cell, function(n_cell, data){
            data %>%  
            select(clone) %>% 
            .[sample(1:nrow(.), round(n_cell), replace = T), ] %>%
            distinct() %>% 
            nrow()
          }, data = data,
        mc.cores = 10,
        ignore.interactive = T) %>% unlist(),
          day_clone = data %>%  
            select(clone) %>% 
            distinct() %>% 
            nrow()
        )
    }
    )) %>% 
  unnest(detected_clone) %>% 
  mutate(
    sample = as.factor(sample),
  ) %>% 
  group_by(donor, antigen, n_cell) %>% 
  nest() %>% 
  mutate(pval = lapply(data, function(data){
    data %>% 
    group_by(day) %>% 
    mutate(
      ecdf = ecdf(detected_clone)(detected_clone)
      ) %>% 
    filter(!duplicated(detected_clone)) %>% 
    group_by(detected_clone) %>% 
    mutate(
      ecdf = ecdf / length(levels(day)),
      s_ecdf = sum(ecdf)) %>% 
    group_by(day) %>%
    mutate(s_ecdf = s_ecdf - ecdf) %>%
    group_by(detected_clone) %>% 
    mutate(pval = max(sum(s_ecdf))) %>%
    pull(pval) %>% 
    max()
  })) %>% 
  unnest(data, pval) %>% 
  group_by(donor, antigen) %>% 
  mutate(pval_signif = max(n_cell[pval > 0.05])) %>% 
  select(-data)

save(data, file = "results/2020_10_30_clone_diversity_bootstrap.Rdata")


p <- ggplot(data %>%
         filter(n_cell < max(pval_signif, day_size))) +
  geom_vline(
    aes(
      xintercept = pval_signif
    ),
    color = "gray50",
    linetype = 1,
    size = 1.5
  ) +
  geom_line(data = data %>% 
              filter(n_cell < day_size),
    aes(
      x = n_cell,
      y = detected_clone,
      color = day,
      group = str_c(sample, day)
    ),
    alpha = 0.1,
  ) +
  scale_fill_viridis_d() +
  geom_smooth(data = data %>%
         filter(n_cell < max(pval_signif, day_size) + 10),
    aes(
      x = n_cell,
      y = detected_clone,
      group = day
    ),
    method = "loess",
    formula = 'y ~ x',
    color = "black",
    se = F
  ) +
  labs(x = "number of cells",
       y = "number of clone detected") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  facet_wrap(~ antigen + donor, scales = "free", ncol = 4)

ggsave(plot = p, filename = "results/2020_10_30_clone_diversity_bootstrap.png", width = 30, height = 15, units = "cm")
ggsave(plot = p, filename = "results/2020_10_30_clone_diversity_bootstrap.pdf", width = 30, height = 15, units = "cm")

load(file = "results/2020_10_29_clone_diversity_bootstrap.Rdata")
