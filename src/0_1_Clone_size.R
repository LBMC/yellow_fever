library("tidyverse") library("readxl") library("vegan")
theme_set(theme_classic())
pal = "Set1"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}

# read frequency sampling data
freq_resp <- read_xlsx("data/2019_01_25_Frequency_Responses_Fig1_jmup.xlsx") %>%
  as_tibble() %>%
  rename(day = "% of Total CD8") %>%
  gather(key = "donor", value = "percent", -"day") %>%
  separate(donor, into = c("donor", "antigen"), sep = " ") %>%
  drop_na() %>%
  mutate(day = strsplit(day, "D") %>% unlist() %>% as.numeric() %>% .[!is.na(.)])
# read clone data
clone <- read_xlsx("data/2019_01_22_Table_S1 copy.xlsx", sheet = 1) %>%
  as_tibble() %>%
  rename(donor = "Donor_ID",
         day = Timepoint,
         antigen = Epitope,
         clone = "Clone Number"
  ) %>%
  mutate(day = replace(day, day %in% "605", "593"),
         day = as.numeric(day)) %>%
  select(donor, day, antigen, clone)
for (i in 2:6) {
  clone <- read_xlsx("data/2019_01_22_Table_S1 copy.xlsx", sheet = i) %>%
    as_tibble() %>%
    rename(donor = "Donor_ID",
          day = Timepoint,
          antigen = Epitope,
          clone = "Clone Number"
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
         clone = replace(clone, clone %in% 0, size_1_name),
         clone = as.factor(clone)
  )
donors <- clone %>% pull(donor) %>% levels()
days <-  clone %>% pull(day) %>% levels()
antigens <- clone %>% pull(antigen) %>% levels()
# reconstruct missing clone
# we create clone that are present at later time-point
dig_clone <- function(infos, donor, antigen, salt = 0){
  infos$donor <- as.vector(infos$donor)
  infos$day <- as.vector(infos$day)
  infos$antigen <- as.vector(infos$antigen)
  infos$clone <- as.vector(infos$clone)
  r_select <- which(infos$antigen %in% antigen & infos$donor %in% donor)
  infos$clone[r_select] <- as.numeric( infos$clone[r_select] ) + salt
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
    clone <- dig_clone(clone, donor, antigen, salt)
    salt <- salt + nrow(clone)
  }
}
clone <- clone %>%
  right_join(clone %>% group_by(clone) %>% count()) %>%
  mutate(n = replace(n, n == max(n), 1))
alpha_f <- tibble(donor = NA, day = NA, antigen = NA,
         alpha = NA, alpha_prec = NA)
for (i in 1:length(donors)) {
  for (j in 1:length(days)) {
    for (k in 1:length(antigens)) {
      if (clone %>%
        filter(donor %in% donors[i],
               day %in% days[j],
               antigen %in% antigens[k]
              ) %>% nrow() > 0 ) {
        alpha_f_fit <- clone %>%
          filter(donor %in% donors[i],
                day %in% days[j],
                antigen %in% antigens[k]
                ) %>%
          group_by(clone) %>%
          count() %>%
          pull(nn) %>%
          fisherfit()
        alpha_f <- tibble(donor = donors[i], day = days[j], antigen = antigens[k],
                alpha = alpha_f_fit$estimate, alpha_prec = alpha_f_fit$estim.prec) %>%
          bind_rows(alpha_f)
        }
    }
  }
}
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

clone %>%
  mutate(percent = as.numeric(as.vector(percent))) %>%
  ggplot() +
  geom_histogram(aes(x = n, group = donor, fill = antigen)) +
  facet_wrap(~donor + antigen + day)
ggsave("results/survival/clone_size_histogram.pdf")

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
  lm(data = ., log10(alpha) ~ antigen) %>%
  anova()

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
