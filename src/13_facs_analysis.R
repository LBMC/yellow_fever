library(tidyverse)
library(readxl)

path <- "data/2020_02_27_clonal progeny_all.xlsx"
data <- read_excel(path, sheet = 1) %>% 
  rename(exp = ...1) %>% 
  rename_all(tolower) %>%   
  pivot_longer(
    cols = -c(exp, day, donor),
    names_to = "cell_type",
    values_to = "size",
    values_drop_na = T
  ) %>% 
  mutate(exp = str_replace(exp, "\\S*\\s(\\S*)\\s\\S*", "\\1"),
         day = factor(day, c(90, 136, 593, 1400)),
         dday = fct_recode(day, early = "90", early = "136", late = "593", late = "1400"),
         exp = factor(exp),
         donor = factor(donor),
         cell_type = factor(cell_type, c("scm", "cd127- emra", "cd127+ emra")),
         size = size - 10,
         size = ifelse(size < 0, 0, size),
         lsize = log1p(size))
data %>% summary()

data %>% 
  ggplot() +
  geom_histogram(aes(x = lsize, fill = donor)) +
  facet_wrap(~day + cell_type, ncol = 3) +
  theme_bw()

data %>%
  filter(lsize != 0) %>% 
  lm(data = ., lsize ~ dday + cell_type + donor) %>% 
  anova()

library(MASS)
data %>%
  mutate(size = round(size)) %>% 
  filter(size != 0) %>% 
  glm.nb(size ~ dday * cell_type * donor, data = .) %>% 
  summary()
data %>%
  mutate(size = round(size)) %>% 
  filter(size != 0) %>% 
  glm.nb(size ~ dday * cell_type * donor, data = .) %>% 
  anova()
data %>%
  mutate(size = round(size)) %>% 
  filter(size != 0) %>% 
  glm.nb(size ~ dday * cell_type * donor, data = .) %>% 
  coef() %>% 
  exp()

data %>%
  mutate(lsize = lsize != 0 ) %>% 
  glm(data = ., lsize ~ dday + cell_type + donor, family = binomial(link='logit')) %>% 
  anova()

library(randomForest)
train <- sample(1:nrow(data), floor(nrow(data) * 0.9))
data_rf <- randomForest(lsize ~ dday + cell_type + donor ,
                        data = data , subset = train)
data_rf
plot(data_rf)

library(censReg)
data_creg <- censReg(lsize ~ dday + cell_type + donor, data = data)
coef(data_creg)
vcov(data_creg)

