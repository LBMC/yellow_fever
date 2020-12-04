rm(list=ls())
setwd("~/projects/mold/yellow_fever")

########################### init ##############################################
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

######################## optional #############################################
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

####################################### end of init ###########################

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


####################### test for time differences accross donor / antigen #####
# fisher alpha slopes analysis
data <- clone %>%
  mutate(day = fct_reorder(day, as.numeric(as.vector(day)))) %>% 
  group_by(donor, day, antigen, clone) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  mutate(
    time = ifelse(as.numeric(as.vector(day)) > 16, "late", "early"),
    time = as.factor(time)
    ) %>% 
  mutate(day = fct_reorder(day, as.numeric(as.vector(day)))) %>% 
  filter(
    !(donor %in% "Donor C" & antigen %in% "A2" & clone == 1)
  ) %>% 
  group_by(donor, day, antigen) %>% 
  select(-percent) %>% 
  mutate(day_size = n()) %>%
  group_by(donor, antigen) %>%
  mutate(days_size = max(day_size)) %>%
  group_by(antigen, donor, time) %>% 
  mutate(time_size = n(),
         min_time_size = min(time_size)) %>% 
  nest() %>% 
  mutate(alpha = lapply(data, function(data){
      n_sample <- 10000
      tibble(
        sample = 1:n_sample,
        min_time_size = data %>% pull(min_time_size) %>% .[1]
      ) %>%
        mutate(
          alpha = pbmcapply::pbmclapply(
            sample, function(sample, data){
            data %>%
            pull(n) %>%
            sample(round(min_time_size), replace = T) %>%
            fisher.alpha()
          },
          data = data,
          mc.cores = 10,
          ignore.interactive = T
        ) %>% unlist()
        ) %>% 
        select(-min_time_size)
    }
    )) %>% 
  ungroup() %>% 
  unnest(data) %>% 
  unnest(alpha) %>% 
  filter(alpha < 2000)

# alpha table without boostrap
clone %>%
  group_by(donor, day, antigen, clone) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  mutate(
    time = ifelse(as.numeric(as.vector(day)) > 16, "late", "early"),
    time = as.factor(time)
    ) %>% 
  mutate(day = fct_reorder(day, as.numeric(as.vector(day)))) %>% 
  filter(
    !(donor %in% "Donor C" & antigen %in% "A2" & clone == 1)
  ) %>% 
  group_by(donor, day, antigen) %>% 
  select(-percent) %>% 
  mutate(day_size = n()) %>%
  group_by(donor, antigen) %>%
  mutate(days_size = max(day_size)) %>%
  group_by(antigen, donor, time) %>% 
  mutate(time_size = n(),
         min_time_size = min(time_size)) %>%
  nest() %>% 
  mutate(alpha = lapply(data, function(data){
      n_sample <- 1
      tibble(
        sample = 1:n_sample,
        min_time_size = data %>% pull(min_time_size) %>% .[1]
      ) %>%
        mutate(
          alpha = pbmcapply::pbmclapply(
            sample, function(sample, data){
            data %>%
              pull(n) %>%
              sample(round(min_time_size), replace = T) %>%
              fisher.alpha()
          },
          data = data,
          mc.cores = 10,
          ignore.interactive = T
        ) %>% unlist()
        ) %>% 
        select(-min_time_size)
    }
    )) %>% 
  ungroup() %>% 
  unnest(data) %>% 
  unnest(alpha) %>% 
  filter(alpha < 2000) %>% 
  write_csv(
    path = "results/2020_11_18_alpha_diversity.csv"
  )
data %>%
  write_csv(
    path = "results/2020_11_18_alpha_diversity_bootstrap.csv"
  )

qqplot_resid <- function(x){
  x %>% 
  broom.mixed::augment() %>% 
  ggplot() +
  geom_qq(aes(sample = .resid)) +
  geom_qq_line(aes(sample = .resid))
}

fitted_vs_value <- function(x){
  x %>% 
  broom.mixed::augment() %>% 
  janitor::clean_names() %>% 
  ggplot() +
  geom_point(
    aes(
      x = alpha,
      y = fitted
    )
  )
}

# 3 way anova
data %>% 
  mutate(alpha = sqrt(alpha + 3/8)) %>% 
  lm(alpha ~ time + antigen + donor, data = .) %>% 
  broom::tidy() %>% 
  mutate(
    estimate = estimate^2 - 3/8,
    std.error = std.error^2 - 3/8,
  ) %>% 
  write_csv(
    path = "results/2020_11_18_alpha_diversity_model_intercept.csv"
  )

# histogram
data %>% 
  ggplot() +
  geom_histogram(
    aes(
      x = alpha,
      fill = time
    ),
    bins = 100
  ) +
  facet_wrap(~antigen + donor, scales = "free_x") +
  scale_x_log10() +
  theme_classic()
ggsave(filename = "results/2020_11_18_alpha_diversity_bootstrap_hist.pdf", width = 30, height = 15, units = "cm")


# boxplot
data %>% 
  ggplot() +
  geom_boxplot(
    aes(
      x = time,
      y = alpha,
      fill = NA,
      color = time
    ),
  ) +
  facet_wrap(~antigen + donor, scales = "free_x") +
  scale_y_log10() +
  theme_classic()
ggsave(filename = "results/2020_11_18_alpha_diversity_bootstrap_boxplot.pdf", width = 30, height = 15, units = "cm")


###################### fisher alpha subsampling experiment ###################
# boostrap
data <- clone %>%
  mutate(day = fct_reorder(day, as.numeric(as.vector(day)))) %>% 
  filter(
    !(donor %in% "Donor C" & antigen %in% "A2" & clone == 1)
  ) %>% 
  group_by(donor, antigen, clone, day) %>% 
  mutate(n = n()) %>%
  group_by(donor, antigen, day) %>% 
  select(-percent) %>% 
  mutate(day_size = n()) %>%
  group_by(donor, antigen) %>%
  mutate(days_size = max(day_size)) %>%
  group_by(donor, antigen, day) %>% 
  nest() %>% 
  mutate(alpha = lapply(data, function(data){
      n_sample <- 1000
      tibble(
        n_cell = seq(
          from = 20,
          to = max(600, (data %>% pull(days_size) %>% max()) * 1.1),
          by = 1) %>%
          rep(n_sample),
        sample = rep(
          1:n_sample,
          each = (
            seq(
              from = 20,
              to = max(600, (data %>% pull(days_size) %>% max()) * 1.1),
              by = 1) %>%
                length()
          )),
        day_size = (data %>% pull(day_size) %>%  max()),
        days_size = (data %>% pull(days_size) %>%  max())
      ) %>%
        mutate(
          alpha = pbmcapply::pbmclapply(n_cell, function(n_cell, data){
            data %>%
            pull(n) %>%
            sample(round(n_cell), replace = T) %>%
            fisher.alpha()
          },
          data = data,
          mc.cores = 10,
          ignore.interactive = T
        ) %>% unlist(),
          day_clone = data %>%
            select(clone) %>%
            distinct() %>%
            nrow()
        )
    }
    )) %>% 
  unnest(alpha) %>% 
  mutate(
    sample = as.factor(sample),
  )

data <- data %>%
  select(-c(pval)) %>%
  group_by(donor, antigen, n_cell) %>% 
  nest() %>% 
  mutate(pval = pbmcapply::pbmclapply(data, function(data){
    data %>% 
    group_by(day) %>% 
    mutate(
      ecdf = ecdf(alpha)(alpha)
      ) %>% 
    filter(!duplicated(alpha)) %>% 
    group_by(alpha) %>% 
    mutate(
      ecdf = ecdf / length(levels(day)),
      s_ecdf = sum(ecdf)) %>% 
    group_by(day) %>%
    mutate(s_ecdf = s_ecdf - ecdf) %>%
    group_by(alpha) %>% 
    mutate(pval = max(sum(s_ecdf))) %>%
    pull(pval) %>% 
    max()
  },
  mc.cores = 10,
  ignore.interactive = T)) %>% 
  unnest(data, pval) %>% 
  group_by(donor, antigen) %>% 
  mutate(pval_signif = max(n_cell[pval > 0.05])) %>% 
  select(-data) %>% 
  group_by(donor, antigen, day, n_cell) %>% 
  mutate(
    alpha_min = quantile(alpha, 0.05),
    alpha_max = quantile(alpha, 0.95)
  )

# p-value computation between early and late
#data <- 

data <- data %>%
  select(-c(pval_early_late)) %>%
  group_by(donor, antigen) %>% 
  mutate(
    time = as.numeric(as.vector(day)),
    time = ifelse(time == min(time),
                  "early",
                  ifelse(time == max(time),
                         "late",
                         NA)),
    time = as.factor(time)
    ) %>% 
  group_by(donor, antigen, n_cell) %>% 
  nest() %>% 
  mutate(pval_early_late = pbmcapply::pbmclapply(
    data, function(data){
    data %>% 
      filter(!is.na(time)) %>% 
      group_by(time) %>% 
      mutate(
        ecdf = ecdf(alpha)(alpha)
        ) %>% 
      filter(!duplicated(alpha)) %>% 
      group_by(alpha) %>% 
      mutate(
        ecdf = ecdf / length(levels(time)),
        s_ecdf = sum(ecdf)) %>% 
      group_by(time) %>%
      mutate(s_ecdf = s_ecdf - ecdf) %>%
      group_by(alpha) %>% 
      mutate(pval = max(sum(s_ecdf))) %>%
      pull(pval) %>% 
      max()
  },
  mc.cores = 10,
  ignore.interactive = T)) %>% 
  unnest(c(data, pval_early_late)) %>% 
  group_by(donor, antigen) %>% 
  mutate(pval_early_late_signif = max(n_cell[pval_early_late > 0.05]))

data %>%
  group_by(donor, antigen, day, n_cell) %>% 
  filter(!(donor == "Donor B" & antigen == "A2" & alpha == max(alpha))) %>%
  filter(donor == "Donor B", antigen == "A2") %>%
  ggplot(
    aes(
    x = n_cell,
    y = alpha,
    )
  ) +
  geom_point() +
  facet_wrap(~day)

data <- data %>%
  filter(alpha < 1e8) %>%
  group_by(donor, antigen, day, n_cell) %>% 
  filter(!(donor == "Donor B" & antigen == "A2" & alpha == max(alpha))) %>%
  group_by(donor, antigen, day, n_cell) %>% 
  mutate(
    alpha_min = quantile(alpha, 0.05),
    alpha_max = quantile(alpha, 0.95)
  )


save(data, file = "results/2021_11_20_fisher_diversity_bootstrap.Rdata")
load(file = "results/2021_11_20_fisher_diversity_bootstrap.Rdata")

# plot

p <- ggplot()+
  geom_vline(data = data %>%
         ungroup() %>%
         filter(pval_signif < days_size * 1.1) %>%
          slice_sample(n = 10000),
    aes(
      xintercept = pval_signif
    ),
    color = "gray50",
    linetype = 1,
    size = 2
  ) +
  geom_label(
    data = data %>%
      ungroup() %>%
      filter(pval_signif < days_size * 1.1) %>%
      slice_sample(n = 10000) %>%
      group_by(antigen, donor) %>%
      select(antigen, donor, pval_signif) %>%
      distinct(),
    aes(
      x = pval_signif,
      y = 30,
      label = pval_signif,
    ),
    fill = "gray50",
    check_overlap = T
  ) +
  geom_vline(data = data %>%
         ungroup() %>%
         filter(pval_early_late_signif < days_size * 1.1) %>%
         slice_sample(n = 10000),
    aes(
      xintercept = pval_early_late_signif
    ),
    color = "black",
    linetype = 1,
    size = 1
  ) +
  geom_label(
    data = data %>%
      ungroup() %>%
      filter(pval_early_late_signif < days_size * 1.1) %>%
      slice_sample(n = 10000) %>%
      group_by(antigen, donor) %>%
      select(antigen, donor, pval_early_late_signif) %>%
      distinct(),
    aes(
      x = pval_early_late_signif,
      y = 0,
      label = pval_early_late_signif
    ),
    check_overlap = T
  ) +
  geom_smooth(data = data %>%
         ungroup() %>%
         filter(n_cell < days_size * 1.1) %>%
         slice_sample(n = 10000),
    aes(
      x = n_cell,
      y = alpha,
      group = day
    ),
    method = "loess",
    formula = 'y ~ x',
    color = "black",
    se = F,
    size = 1.5
  ) +
  geom_ribbon(data = data %>% 
         ungroup() %>%
         filter(n_cell < day_size) %>%
         slice_sample(n = 10000),
    aes(
      x = n_cell,
      ymin = alpha_min,
      ymax = alpha_max,
      fill = day,
      group = day
    ),
    alpha = 0.6
  ) +
  geom_smooth(data = data %>%
         ungroup() %>%
         filter(n_cell < days_size * 1.1) %>%
         slice_sample(n = 10000),
    aes(
      x = n_cell,
      y = alpha,
      color = day,
      group = day
    ),
    method = "loess",
    formula = 'y ~ x',
    se = F,
    size = 0.5
  ) +
  labs(x = "number of cells",
       y = "Fisher's Alpha") +
  guides(
    colour = guide_legend(override.aes = list(alpha = 1)),
    fill = guide_legend(override.aes = list(alpha = 1))
    ) +
  facet_wrap(~ antigen + donor, ncol = 4, scales = "free") +
  theme_classic()
print(p)

ggsave(plot = p, filename = "results/2020_11_20_alpha_diversity_bootstrap.png", width = 30, height = 15, units = "cm")
ggsave(plot = p, filename = "results/2020_11_20_alpha_diversity_bootstrap.pdf", width = 30, height = 15, units = "cm")
  

###################### clone number subsampling experiment ###################

# boostrap
data <- clone %>%
  mutate(day = fct_reorder(day, as.numeric(as.vector(day)))) %>% 
  filter(
    !(donor %in% "Donor C" & antigen %in% "A2" & clone == 1)
  ) %>% 
  group_by(donor, day, antigen) %>% 
  select(-percent) %>% 
  mutate(day_size = n()) %>%
  group_by(donor, antigen) %>%
  mutate(days_size = max(day_size)) %>%
  group_by(donor, day, antigen) %>% 
  nest() %>% 
  mutate(detected_clone = lapply(data, function(data){
      n_sample <- 1000
      tibble(
        n_cell = seq(
          from = 20,
          to = max(1000, (data %>% pull(days_size) %>% max()) + 10),
          by = 1) %>%
          rep(n_sample),
        sample = rep(
          1:n_sample,
          each = (
            seq(
              from = 20,
              to = max(1000, (data %>% pull(days_size) %>% max()) + 10),
              by = 1) %>%
                length()
          )),
        day_size = (data %>% pull(day_size) %>%  max()),
        days_size = (data %>% pull(days_size) %>%  max())
      ) %>%
        mutate(
          detected_clone = pbmcapply::pbmclapply(n_cell, function(n_cell, data){
            data %>%
            select(clone) %>%
            .[sample(1:nrow(.), round(n_cell), replace = T), ] %>%
            distinct() %>%
            nrow()
          },
          data = data,
          mc.cores = 10,
          ignore.interactive = T
        ) %>% unlist(),
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
  mutate(pval = pbmcapply::pbmclapply(data, function(data){
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
  },
  mc.cores = 10,
  ignore.interactive = T)) %>% 
  unnest(data, pval) %>% 
  group_by(donor, antigen) %>% 
  mutate(pval_signif = max(n_cell[pval > 0.05])) %>% 
  select(-data) %>% 
  group_by(donor, antigen, day, n_cell) %>% 
  mutate(
    detected_clone_min = quantile(detected_clone, 0.05),
    detected_clone_max = quantile(detected_clone, 0.95)
  )

# p-value computation between early and late
data <- data %>%
  select(-c(pval_early_late)) %>%
  group_by(donor, antigen) %>% 
  mutate(
    time = as.numeric(as.vector(day)),
    time = ifelse(time == min(time),
                  "early",
                  ifelse(time == max(time),
                         "late",
                         NA)),
    time = as.factor(time)
    ) %>% 
  group_by(donor, antigen, n_cell) %>% 
  nest() %>% 
  mutate(pval_early_late = pbmcapply::pbmclapply(
    data, function(data){
    data %>% 
      filter(!is.na(time)) %>% 
      group_by(time) %>% 
      mutate(
        ecdf = ecdf(detected_clone)(detected_clone)
        ) %>% 
      filter(!duplicated(detected_clone)) %>% 
      group_by(detected_clone) %>% 
      mutate(
        ecdf = ecdf / length(levels(time)),
        s_ecdf = sum(ecdf)) %>% 
      group_by(time) %>%
      mutate(s_ecdf = s_ecdf - ecdf) %>%
      group_by(detected_clone) %>% 
      mutate(pval = max(sum(s_ecdf))) %>%
      pull(pval) %>% 
      max()
  },
  mc.cores = 10,
  ignore.interactive = T)) %>% 
  unnest(c(data, pval_early_late)) %>% 
  group_by(donor, antigen) %>% 
  mutate(pval_early_late_signif = max(n_cell[pval_early_late > 0.05]))

save(data, file = "results/2020_11_20_clone_diversity_bootstrap.Rdata")
load(file = "results/2020_11_20_clone_diversity_bootstrap.Rdata")


# plot
p <- ggplot() +
  geom_vline(data = data %>%
         filter(pval_signif < days_size * 1.1) %>%
         slice_sample(n = 100000),
    aes(
      xintercept = pval_signif
    ),
    color = "gray50",
    linetype = 1,
    size = 2
  ) +
  geom_label(
    data = data %>%
      ungroup() %>%
      filter(pval_signif < days_size * 1.1) %>%
      slice_sample(n = 10000) %>%
      group_by(antigen, donor) %>%
      select(antigen, donor, pval_signif) %>%
      distinct(),
    aes(
      x = pval_signif,
      y = 30,
      label = pval_signif,
    ),
    fill = "gray50",
    check_overlap = T
  ) +
  geom_vline(data = data %>%
         ungroup() %>%
         filter(pval_early_late_signif < days_size * 1.1) %>%
         slice_sample(n = 100000),
    aes(
      xintercept = pval_early_late_signif
    ),
    color = "black",
    linetype = 1,
    size = 1
  ) +
  geom_label(
    data = data %>%
      ungroup() %>%
      filter(pval_early_late_signif < days_size * 1.1) %>%
      slice_sample(n = 10000) %>%
      group_by(antigen, donor) %>%
      select(antigen, donor, pval_early_late_signif) %>%
      distinct(),
    aes(
      x = pval_early_late_signif,
      y = 0,
      label = pval_early_late_signif
    ),
    check_overlap = T
  ) +
  geom_smooth(data = data %>%
         ungroup() %>%
         filter(n_cell < days_size * 1.1) %>%
         slice_sample(n = 100000),
    aes(
      x = n_cell,
      y = detected_clone,
      group = day
    ),
    method = "loess",
    formula = 'y ~ x',
    color = "black",
    se = F,
    size = 1.5
  ) +
  geom_ribbon(data = data %>% 
              ungroup() %>%
              filter(n_cell < day_size) %>%
              slice_sample(n = 1000000),
    aes(
      x = n_cell,
      ymin = detected_clone_min,
      ymax = detected_clone_max,
      fill = day,
      group = day
    ),
    alpha = 0.6
  ) +
  geom_smooth(data = data %>%
         ungroup() %>%
         filter(n_cell < days_size * 1.1) %>%
         slice_sample(n = 100000),
    aes( 
      x = n_cell,
      y = detected_clone,
      color = day,
      group = day
    ),
    method = "loess",
    formula = 'y ~ x',
    se = F,
    size = 0.5
  ) +
  labs(x = "number of cells",
       y = "number of clone detected") +
  guides(colour = guide_legend(override.aes = list(alpha = 1)),
         fill = guide_legend(override.aes = list(alpha = 1))) +
  facet_wrap(~ antigen + donor, scales = "free", ncol = 4) +
  theme_classic()
print(p)

ggsave(plot = p, filename = "results/2020_11_20_clone_diversity_bootstrap.pdf", width = 30, height = 15, units = "cm")
ggsave(plot = p, filename = "results/2020_11_20_clone_diversity_bootstrap.png", width = 30, height = 15, units = "cm")

############################### FIg 1 boostraped ##############################

data <- clone %>%
  mutate(day = fct_reorder(day, as.numeric(as.vector(day)))) %>% 
  group_by(donor, day, antigen, clone) %>% 
  mutate(n = n()) %>%
  group_by(donor, day, antigen) %>% 
  group_by(donor, day, antigen, percent) %>% 
  nest() %>% 
  mutate(alpha = pbmcapply::pbmclapply(data, function(data){
      n_sample <- 1000
      tibble(
        sampling = nrow(data) %>% rep(n_sample),
        sample = seq(1:n_sample)
        ) %>%
        mutate(
          alpha = data %>%
            pull(n) %>%
            fisher.alpha(),
          alpha_boot = lapply(sampling, function(sampling, data){
            data %>%
            pull(n) %>%
            sample(round(sampling), replace = T) %>%
            fisher.alpha()
          }, data = data) %>% unlist(),
          alpha_min = quantile(alpha_boot, 0.05),
          alpha_max = quantile(alpha_boot, 0.95)
        )
    },
    mc.cores = 10,
    ignore.interactive = T)
  ) %>% 
  unnest(alpha) %>% 
  unnest(data)
  mutate(percent = as.numeric(as.vector(percent)))

data %>%
  write.csv(file="results/survival/2020_12_04_fisher_vs_sampling_vs_time.csv")

data  %>%
  ggplot() +
  geom_point(aes(x = day, y = alpha, group = donor, color = antigen, shape = donor), size = 4) +
  geom_linerange(
    aes(x = day,
        ymin = alpha_min,
        ymax = alpha_max, group = donor),
    color = "gray50"
  ) +
  geom_line(aes(x = day, y = alpha, group = paste0(donor, antigen), color = antigen)) +
  labs(y = "Fisher's Alpha")
ggsave("results/survival/2020_11_05_fisher_vs_sampling_vs_time.pdf")

data <- clone %>%
  mutate(day = fct_reorder(day, as.numeric(as.vector(day)))) %>% 
  group_by(donor, day, antigen, clone) %>% 
  mutate(n = n()) %>%
  group_by(donor, day, antigen) %>% 
  group_by(donor, day, antigen, percent) %>% 
  nest() %>% 
  mutate(alpha = pbmcapply::pbmclapply(data, function(data){
      n_sample <- 1000
      tibble(
        sampling = nrow(data) %>% rep(n_sample),
        sample = seq(1:n_sample)
        ) %>%
        mutate(
          alpha = data %>%
            pull(n) %>%
            fisher.alpha(),
          alpha_boot = lapply(sampling, function(sampling, data){
            data %>%
            pull(n) %>%
            sample(round(sampling), replace = T) %>%
            fisher.alpha()
          }, data = data) %>% unlist(),
          alpha_min = quantile(alpha_boot, 0.05),
          alpha_max = quantile(alpha_boot, 0.95),
        )
    },
    mc.cores = 10,
    ignore.interactive = T )) %>% 
  unnest(alpha) %>% 
  unnest(data) %>%
  mutate(percent = as.numeric(as.vector(percent)))

data %>%
  write.csv(file="results/survival/2020_12_04_fisher_vs_sampling_log10.csv")

data %>%
  ggplot() +
  geom_point(aes(y = alpha, x = percent, group = donor, color = antigen, shape = donor), size = 4) +
  geom_linerange(
    aes(x = percent,
        ymin = alpha_min,
        ymax = alpha_max, group = donor),
    color = "gray50"
  ) +
  geom_line(aes(y = alpha, x = percent, group = paste0(donor, antigen), color = antigen)) +
  scale_x_log10() +
  scale_y_log10() +
  labs(y = "Fisher's Alpha",
       x = "percentage of CD8+ T cells"
  )
ggsave("results/survival/2020_11_05_fisher_vs_sampling_log10.pdf")

data <- clone %>%
  mutate(day = fct_reorder(day, as.numeric(as.vector(day)))) %>% 
  group_by(donor, day, antigen, clone) %>% 
  mutate(n = n()) %>%
  group_by(donor, day, antigen) %>% 
  filter(day == 15) %>% 
  nest() %>% 
  mutate(alpha = pbmcapply::pbmclapply(data, function(data){
      n_sample <- 1000
      tibble(
        sampling = nrow(data) %>% rep(n_sample),
        sample = seq(1:n_sample)
        ) %>%
        mutate(
          alpha = data %>%
            pull(n) %>%
            fisher.alpha(),
          alpha_boot = lapply(sampling, function(sampling, data){
            data %>%
            pull(n) %>%
            sample(round(sampling), replace = T) %>%
            fisher.alpha()
          }, data = data) %>% unlist(),
          alpha_min = quantile(alpha_boot, 0.05),
          alpha_max = quantile(alpha_boot, 0.95),
        )
    },
    mc.cores = 10,
    ignore.interactive = T )) %>% 
  unnest(alpha) %>% 
  unnest(data) %>%
  mutate(percent = as.numeric(as.vector(percent)))

data %>%
  write.csv(file="results/survival/2020_12_04_fisher_vs_sampling_log10_D15.csv")

data %>%
  ggplot() +
  geom_point(aes(y = alpha, x = percent, group = donor, color = antigen, shape = donor), size = 4) +
  geom_linerange(
    aes(x = percent,
        ymin = alpha_min,
        ymax = alpha_max, group = donor),
    color = "gray50"
  ) +
  geom_line(aes(y = alpha, x = percent, group = paste0(donor, antigen), color = antigen)) +
  labs(y = "Fisher's Alpha",
       x = "percentage of CD8+ T cells"
  )
ggsave("results/survival/2020_11_05_fisher_vs_sampling_log10_D15.pdf")
