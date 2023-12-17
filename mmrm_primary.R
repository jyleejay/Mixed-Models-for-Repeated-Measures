library(tidyverse)
library(sas7bdat)
library(officer)
library(reshape2)
library(haven)
library(gtsummary)
library(mmrm)
library(emmeans)
library(ggplot2)
library(kableExtra)
library(ggeasy)
library(flextable)
theme_gtsummary_compact(set_theme = TRUE, font_size = NULL)

# set working directory
setwd("C:/Users/hi/Desktop/A002-A4")

# Save CSV Files
list.files('sasdata')

# CSV file exist
subDir  <- "CSV"
sas.files <- list.files('sasdata')[grepl(".sas7bdat", list.files('sasdata'))]
for(files_ in sas.files){
  name_ <- toupper(strsplit(files_, ".", fixed = TRUE)[[1]][1])
  outfile <- write.csv(as.data.frame(read_sas(paste0('sasdata/',files_)), NULL),
                       file = paste0(file.path(subDir), "/", name_, ".csv"))
}

# list.files(getwd())

# Load CSV Files
data.path <- paste0(getwd(), "/A002-A4dataset")
data.files <- list.files(data.path)
# data.files

# 1. PANSS total score (change from baseline by visit)
panss <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<adpanss\\>", data.files)]), NULL))
# move placebo into first
panss$ARM <- fct_relevel(panss$ARM, "Placebo")

# whole study period
whole <- panss %>% 
  filter(FASFL == "Y") %>%
  filter(XPNSID == 33)  

# stack the table
entire_panss <- whole %>%
  select(ARM, AVISIT, XPNSSCRE, XPNSCHG) %>% 
  filter(!is.na(AVISIT)) %>% 
  mutate(AVISIT = factor(AVISIT, levels = c(-1, 10, 20, 30, 40, 50, 60), labels = c("Baseline", "Week 1", "Week 2", "Week 3", "Week 4", "Week 5", "Week 6"))) %>% 
  tbl_strata(strata = c(AVISIT),
             .tbl_fun =
               ~.x %>%
               tbl_summary(by = ARM,
                           type = list(where(is.numeric) ~ "continuous2"),
                           missing = "no",
                           statistic = list(all_continuous() ~ c("{mean} ± ({sd})",
                                                                 "{median}",
                                                                 "({min}, {max})")),
                           digits = everything() ~ 1),
             .combine_with = "tbl_stack")
entire_panss <- as_flex_table(entire_panss)

# -----------------------------------------------------
# MMRM
panss_mmrm <- panss %>%
  filter(FASFL == "Y") %>%
  group_by(SUBJID) %>% 
  filter(XPNSID == 33) %>% 
  filter(AVISIT != -1) %>% 
  filter(!(is.na(XPNSSCRE))) %>%
  select(SUBJID, ARM, REGION, AVISIT, XPNSSCRBL, XPNSCHG) %>%
  mutate_at(vars(REGION, AVISIT, ARM), as.factor)

# fit the model
panss_fit <- mmrm(
  formula = XPNSCHG ~ ARM + REGION + AVISIT + ARM:AVISIT + XPNSSCRBL + XPNSSCRBL:AVISIT + us(AVISIT|SUBJID),
  data = panss_mmrm,
  method = "Kenward-Roger")

# lsmean by ARM
panss_lsm_value <- emmeans(panss_fit, ~ ARM|AVISIT)
lsmean_panss <- summary(panss_lsm_value)
lsmean_panss <- flextable(lsmean_panss)

# lsmean difference
panss_fit.rg <- ref_grid(panss_fit)
panss_lsm <- pairs(emmeans(panss_fit.rg, ~ ARM|AVISIT), adjust = "none", reverse = T)

panss_lsm # for p-value
confint(panss_lsm) # display 95% CI

# make a df for lsmean differ table
ci_lsm <- confint(panss_lsm) %>%
  arrange(AVISIT, contrast) %>% 
  select(contrast, AVISIT, lower.CL, upper.CL) %>% 
  filter(str_detect(contrast, "Placebo"))

df.lsm <- summary(panss_lsm) %>%
  arrange(AVISIT, contrast) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  left_join(ci_lsm, by = c("contrast", "AVISIT")) %>% 
  mutate(AVISIT = factor(AVISIT, levels = c(10, 20, 30, 40, 50, 60), labels = c("Week 1", "Week 2", "Week 3", "Week 4", "Week 5", "Week 6")))

# generate tables by visit
lsm_table <- flextable(df.lsm)

########################## make plots by country
# Korea
# stack the table
korea_panss <- whole %>%
  filter(REGION == 2) %>%
  select(ARM, AVISIT, XPNSSCRE, XPNSCHG) %>% 
  mutate(AVISIT = factor(AVISIT, levels = c(-1, 10, 20, 30, 40, 50, 60), labels = c("Baseline", "Week 1", "Week 2", "Week 3", "Week 4", "Week 5", "Week 6"))) %>% 
  tbl_strata(strata = c(AVISIT),
             .tbl_fun =
               ~.x %>%
               tbl_summary(by = ARM,
                           type = list(where(is.numeric) ~ "continuous2"),
                           missing = "no",
                           statistic = list(all_continuous() ~ c("{mean} ± ({sd})",
                                                                 "{median}",
                                                                 "({min}, {max})")),
                           digits = everything() ~ 1),
             .combine_with = "tbl_stack")
korea_panss <- as_flex_table(korea_panss)

# MMRM
panss_k <- panss %>%
  filter(FASFL == "Y") %>%
  group_by(SUBJID) %>% 
  filter(XPNSID == 33) %>% 
  filter(AVISIT != -1) %>%
  filter(REGION == 2) %>%
  filter(!(is.na(XPNSSCRE))) %>%
  select(SUBJID, ARM, AVISIT, XPNSSCRBL, XPNSCHG) %>%
  mutate_at(vars(AVISIT, ARM), as.factor)

# fit the model
panss_fit_k <- mmrm(
  formula = XPNSCHG ~ ARM + AVISIT + ARM:AVISIT + XPNSSCRBL + XPNSSCRBL:AVISIT + us(AVISIT|SUBJID),
  data = panss_k,
  method = "Kenward-Roger")

# lsmean by ARM
panss_lsm_value_k <- emmeans(panss_fit_k, ~ ARM|AVISIT)
lsmean_panss_k <- summary(panss_lsm_value_k)
lsmean_panss_k <- flextable(lsmean_panss_k)

# lsmean difference
panss_fit.rg_k <- ref_grid(panss_fit_k)
panss_lsm_k <- pairs(emmeans(panss_fit.rg_k, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_k <- confint(panss_lsm_k) %>%
  arrange(AVISIT, contrast) %>% 
  select(contrast, AVISIT, lower.CL, upper.CL) %>% 
  filter(str_detect(contrast, "Placebo"))

df.lsm_k <- summary(panss_lsm_k) %>%
  arrange(AVISIT, contrast) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  left_join(ci_lsm_k, by = c("contrast", "AVISIT")) %>% 
  mutate(AVISIT = factor(AVISIT, levels = c(10, 20, 30, 40, 50, 60), labels = c("Week 1", "Week 2", "Week 3", "Week 4", "Week 5", "Week 6")))

# generate tables by visit
lsm_table_k <- flextable(df.lsm_k)

# dataset for plot
df_plot_k <- summary(panss_lsm_value_k) %>% 
  select(-c(SE, df)) %>% 
  mutate(AVISIT = factor(AVISIT, levels = c(10, 20, 30, 40, 50, 60), labels = c("Week 1", "Week 2", "Week 3", "Week 4", "Week 5", "Week 6")))

# plotting
plot_k <- ggplot(df_plot_k, aes(x = AVISIT, y = emmean, group = ARM, color = ARM)) +
  theme_minimal() +
  geom_line() +
  geom_point(aes(y = emmean), size = 3, shape = 21) +
  easy_remove_legend_title() +
  labs(x = "Visit number", y = "lsmean", title = "lsmean result over time - Korea") +
  ggeasy::easy_center_title() +
  easy_plot_title_size(20)

# Japan
# stack the table
japan_panss <- whole %>%
  filter(REGION == 1) %>%
  select(ARM, AVISIT, XPNSSCRE, XPNSCHG) %>% 
  mutate(AVISIT = factor(AVISIT, levels = c(-1, 10, 20, 30, 40, 50, 60), labels = c("Baseline", "Week 1", "Week 2", "Week 3", "Week 4", "Week 5", "Week 6"))) %>% 
  tbl_strata(strata = c(AVISIT),
             .tbl_fun =
               ~.x %>%
               tbl_summary(by = ARM,
                           type = list(where(is.numeric) ~ "continuous2"),
                           missing = "no",
                           statistic = list(all_continuous() ~ c("{mean} ± ({sd})",
                                                                 "{median}",
                                                                 "({min}, {max})")),
                           digits = everything() ~ 1),
             .combine_with = "tbl_stack")
japan_panss <- as_flex_table(japan_panss)

# MMRM
panss_j <- panss %>%
  filter(FASFL == "Y") %>%
  group_by(SUBJID) %>% 
  filter(XPNSID == 33) %>% 
  filter(AVISIT != -1) %>%
  filter(REGION == 1) %>%
  filter(!(is.na(XPNSSCRE))) %>%
  select(SUBJID, ARM, AVISIT, XPNSSCRBL, XPNSCHG) %>%
  mutate_at(vars(AVISIT, ARM), as.factor)

# fit the model
panss_fit_j <- mmrm(
  formula = XPNSCHG ~ ARM + AVISIT + ARM:AVISIT + XPNSSCRBL + XPNSSCRBL:AVISIT + us(AVISIT|SUBJID),
  data = panss_j,
  method = "Kenward-Roger")

# lsmean by ARM
panss_lsm_value_j <- emmeans(panss_fit_j, ~ ARM|AVISIT)
lsmean_panss_j <- summary(panss_lsm_value_j)
lsmean_panss_j <- flextable(lsmean_panss_j)

# lsmean difference
panss_fit.rg_j <- ref_grid(panss_fit_j)
panss_lsm_j <- pairs(emmeans(panss_fit.rg_j, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_j <- confint(panss_lsm_j) %>%
  arrange(AVISIT, contrast) %>% 
  select(contrast, AVISIT, lower.CL, upper.CL) %>% 
  filter(str_detect(contrast, "Placebo"))

df.lsm_j <- summary(panss_lsm_j) %>%
  arrange(AVISIT, contrast) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  left_join(ci_lsm_j, by = c("contrast", "AVISIT")) %>% 
  mutate(AVISIT = factor(AVISIT, levels = c(10, 20, 30, 40, 50, 60), labels = c("Week 1", "Week 2", "Week 3", "Week 4", "Week 5", "Week 6")))

# generate tables by visit
lsm_table_j <- flextable(df.lsm_j)

# dataset for plot
df_plot_j <- summary(panss_lsm_value_j) %>% 
  select(-c(SE, df)) %>% 
  mutate(AVISIT = factor(AVISIT, levels = c(10, 20, 30, 40, 50, 60), labels = c("Week 1", "Week 2", "Week 3", "Week 4", "Week 5", "Week 6")))

# plotting
plot_j <- ggplot(df_plot_j, aes(x = AVISIT, y = emmean, group = ARM, color = ARM)) +
  theme_minimal() +
  geom_line() +
  geom_point(aes(y = emmean), size = 3, shape = 21) +
  easy_remove_legend_title() +
  labs(x = "Visit number", y = "lsmean", title = "lsmean result over time - Japan") +
  ggeasy::easy_center_title() +
  easy_plot_title_size(20)


# 2. CGI-I total score (change from baseline at w6)
cgi <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<adcgi\\>", data.files)]), NULL))
# move placebo into first
cgi$ARM <- fct_relevel(cgi$ARM, "Placebo")

cgi_w6 <- cgi %>%
  filter(FASFL == "Y") %>%
  filter(XCGICD == 2) %>% 
  filter(AVISIT == 60) %>%
  mutate_at(vars(XCGISCR), as.numeric)

tb_cgi_w6 <- cgi_w6 %>%
  select(ARM, XCGISCR) %>% 
  tbl_strata(strata = c(ARM),
             .tbl_fun =
               ~.x %>%
               tbl_summary(type = list(where(is.numeric) ~ "continuous2"),
                           missing = "no",
                           statistic = list(all_continuous() ~ c("{mean} ± ({sd})",
                                                                 "{median}",
                                                                 "({min}, {max})")),
                           digits = everything() ~ 1) %>%
               add_n(statistic = "{N_nonmiss}",
                     col_label = "**N responses**")%>%
               modify_caption("**CGI-I Score at Week 6(MMRM)**"))

tb_cgi_w6 <- as_flex_table(tb_cgi_w6)

# MMRM
cgi_mmrm <- cgi %>%
  filter(FASFL == "Y") %>%
  group_by(SUBJID) %>% 
  filter(XCGICD == 2) %>% 
  filter(AVISIT != -1) %>% 
  filter(!(is.na(XCGISCR))) %>%
  select(SUBJID, ARM, REGION, AVISIT, XCGISCRBL, XCGISCR) %>%
  mutate_at(vars(REGION, AVISIT, ARM), as.factor)

# fit the model
cgi_fit <- mmrm(
  formula = XCGISCR ~ ARM + REGION + AVISIT + ARM:AVISIT + XCGISCRBL + XCGISCRBL:AVISIT + us(AVISIT|SUBJID),
  data = cgi_mmrm,
  method = "Kenward-Roger")

# lsmean by ARM
cgi_lsm_value <- emmeans(cgi_fit, ~ ARM|AVISIT)
lsmean_cgi <- summary(cgi_lsm_value) %>% filter(AVISIT == 60)
lsmean_cgi <- flextable(lsmean_cgi)

# lsmean difference
cgi_fit.rg <- ref_grid(cgi_fit)
cgi_lsm <- pairs(emmeans(cgi_fit.rg, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_cgi <- confint(cgi_lsm) %>%
  filter(AVISIT == 60) %>% 
  select(contrast, lower.CL, upper.CL)

df.lsm_cgi <- summary(cgi_lsm) %>%
  filter(AVISIT == 60) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  left_join(ci_lsm_cgi, by = c("contrast"))

# generate a table
lsm_table_cgi <- flextable(df.lsm_cgi)

# Korea
tb_cgi_w6_k <- cgi_w6 %>%
  filter(REGION == 2) %>%
  select(ARM, XCGISCR) %>% 
  tbl_strata(strata = c(ARM),
             .tbl_fun =
               ~.x %>%
               tbl_summary(type = list(where(is.numeric) ~ "continuous2"),
                           missing = "no",
                           statistic = list(all_continuous() ~ c("{mean} ± ({sd})",
                                                                 "{median}",
                                                                 "({min}, {max})")),
                           digits = everything() ~ 1) %>%
               add_n(statistic = "{N_nonmiss}",
                     col_label = "**N responses**")%>%
               modify_caption("**CGI-I Score at Week 6(MMRM)**"))

tb_cgi_w6_k <- as_flex_table(tb_cgi_w6_k)

# MMRM
cgi_mmrm_k <- cgi %>%
  filter(FASFL == "Y") %>%
  group_by(SUBJID) %>% 
  filter(XCGICD == 2) %>% 
  filter(AVISIT != -1) %>%
  filter(REGION == 2) %>%
  filter(!(is.na(XCGISCR))) %>%
  select(SUBJID, ARM, AVISIT, XCGISCRBL, XCGISCR) %>%
  mutate_at(vars(AVISIT, ARM), as.factor)

# fit the model
cgi_fit_k <- mmrm(
  formula = XCGISCR ~ ARM + AVISIT + ARM:AVISIT + XCGISCRBL + XCGISCRBL:AVISIT + us(AVISIT|SUBJID),
  data = cgi_mmrm_k,
  method = "Kenward-Roger")

# lsmean by ARM
cgi_lsm_value_k <- emmeans(cgi_fit_k, ~ ARM|AVISIT)
lsmean_cgi_k <- summary(cgi_lsm_value_k) %>% filter(AVISIT == 60)
lsmean_cgi_k <- flextable(lsmean_cgi_k)

# lsmean difference
cgi_fit.rg_k <- ref_grid(cgi_fit_k)
cgi_lsm_k <- pairs(emmeans(cgi_fit.rg_k, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_cgi_k <- confint(cgi_lsm_k) %>%
  filter(AVISIT == 60) %>% 
  select(contrast, lower.CL, upper.CL)

df.lsm_cgi_k <- summary(cgi_lsm_k) %>%
  filter(AVISIT == 60) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  left_join(ci_lsm_cgi_k, by = c("contrast"))

# generate a table
lsm_table_cgi_k <- flextable(df.lsm_cgi_k)

# Japan
tb_cgi_w6_j <- cgi_w6 %>%
  filter(REGION == 1) %>%
  select(ARM, XCGISCR) %>% 
  tbl_strata(strata = c(ARM),
             .tbl_fun =
               ~.x %>%
               tbl_summary(type = list(where(is.numeric) ~ "continuous2"),
                           missing = "no",
                           statistic = list(all_continuous() ~ c("{mean} ± ({sd})",
                                                                 "{median}",
                                                                 "({min}, {max})")),
                           digits = everything() ~ 1) %>%
               add_n(statistic = "{N_nonmiss}",
                     col_label = "**N responses**")%>%
               modify_caption("**CGI-I Score at Week 6(MMRM)**"))

tb_cgi_w6_j <- as_flex_table(tb_cgi_w6_j)

# MMRM
cgi_mmrm_j <- cgi %>%
  filter(FASFL == "Y") %>%
  group_by(SUBJID) %>% 
  filter(XCGICD == 2) %>% 
  filter(AVISIT != -1) %>%
  filter(REGION == 1) %>%
  filter(!(is.na(XCGISCR))) %>%
  select(SUBJID, ARM, AVISIT, XCGISCRBL, XCGISCR) %>%
  mutate_at(vars(AVISIT, ARM), as.factor)

# fit the model
cgi_fit_j <- mmrm(
  formula = XCGISCR ~ ARM + AVISIT + ARM:AVISIT + XCGISCRBL + XCGISCRBL:AVISIT + us(AVISIT|SUBJID),
  data = cgi_mmrm_j,
  method = "Kenward-Roger")

# lsmean by ARM
cgi_lsm_value_j <- emmeans(cgi_fit_j, ~ ARM|AVISIT)
lsmean_cgi_j <- summary(cgi_lsm_value_j) %>% filter(AVISIT == 60)
lsmean_cgi_j <- flextable(lsmean_cgi_j)

# lsmean difference
cgi_fit.rg_j <- ref_grid(cgi_fit_j)
cgi_lsm_j <- pairs(emmeans(cgi_fit.rg_j, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_cgi_j <- confint(cgi_lsm_j) %>%
  filter(AVISIT == 60) %>% 
  select(contrast, lower.CL, upper.CL)

df.lsm_cgi_j <- summary(cgi_lsm_j) %>%
  filter(AVISIT == 60) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  left_join(ci_lsm_cgi_j, by = c("contrast"))

# generate a table
lsm_table_cgi_j <- flextable(df.lsm_cgi_j)

####################################################
# 3.NSA-16 total score (change from baseline at w6)
nsa <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<adnsa\\>", data.files)]), NULL))
# move placebo into first
nsa$ARM <- fct_relevel(nsa$ARM, "Placebo")

nsa_w6 <- nsa %>%
  filter(FASFL == "Y") %>%
  filter(XNSACD == 18) %>% 
  filter(AVISIT == 60)

tb_nsa_w6 <- nsa_w6 %>%
  select(ARM, XNSASCR, XNSACHG) %>% 
  tbl_strata(strata = c(ARM),
             .tbl_fun =
               ~.x %>%
               tbl_summary(type = list(where(is.numeric) ~ "continuous2"),
                           missing = "no",
                           statistic = list(all_continuous() ~ c("{mean} ± ({sd})",
                                                                 "{median}",
                                                                 "({min}, {max})")),
                           digits = everything() ~ 1) %>%
               add_n(statistic = "{N_nonmiss}",
                     col_label = "**N responses**")%>%
               modify_caption("**NSA-16 Total Score at Week 6(MMRM)**"))

tb_nsa_w6 <- as_flex_table(tb_nsa_w6)

# MMRM
nsa_mmrm <- nsa %>%
  filter(FASFL == "Y") %>%
  group_by(SUBJID) %>% 
  filter(XNSACD == 18) %>% 
  filter(AVISIT != -1) %>% 
  filter(!(is.na(XNSASCR))) %>%
  select(SUBJID, ARM, REGION, AVISIT, XNSASCRBL, XNSACHG) %>%
  mutate_at(vars(REGION, AVISIT, ARM), as.factor)

# fit the model
nsa_fit <- mmrm(
  formula = XNSACHG ~ ARM + REGION + AVISIT + ARM:AVISIT + XNSASCRBL + XNSASCRBL:AVISIT + us(AVISIT|SUBJID),
  data = nsa_mmrm,
  method = "Kenward-Roger")

# lsmean by ARM
nsa_lsm_value <- emmeans(nsa_fit, ~ ARM|AVISIT)
lsm_nsa <- summary(nsa_lsm_value) %>% filter(AVISIT == 60)
lsm_nsa <- flextable(lsm_nsa)

# lsmean difference
nsa_fit.rg <- ref_grid(nsa_fit)
lsmean_nsa <- pairs(emmeans(nsa_fit.rg, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_nsa <- confint(lsmean_nsa) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(contrast, lower.CL, upper.CL) %>% 
  filter(str_detect(contrast, "Placebo"))

df.lsm_nsa <- summary(lsmean_nsa) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  left_join(ci_lsm_nsa, by = c("contrast"))

# generate tables by visit
lsm_table_nsa <- flextable(df.lsm_nsa)

# KOREA
tb_nsa_w6_k <- nsa_w6 %>%
  filter(REGION == 2) %>%
  select(ARM, XNSASCR, XNSACHG) %>% 
  tbl_strata(strata = c(ARM),
             .tbl_fun =
               ~.x %>%
               tbl_summary(type = list(where(is.numeric) ~ "continuous2"),
                           missing = "no",
                           statistic = list(all_continuous() ~ c("{mean} ± ({sd})",
                                                                 "{median}",
                                                                 "({min}, {max})")),
                           digits = everything() ~ 1) %>%
               add_n(statistic = "{N_nonmiss}",
                     col_label = "**N responses**")%>%
               modify_caption("**NSA-16 Total Score at Week 6(MMRM)**"))

tb_nsa_w6_k <- as_flex_table(tb_nsa_w6_k)

# MMRM
nsa_mmrm_k <- nsa %>%
  filter(FASFL == "Y") %>%
  filter(REGION == 2) %>%
  group_by(SUBJID) %>% 
  filter(XNSACD == 18) %>% 
  filter(AVISIT != -1) %>% 
  filter(!(is.na(XNSASCR))) %>%
  select(SUBJID, ARM, AVISIT, XNSASCRBL, XNSACHG) %>%
  mutate_at(vars(AVISIT, ARM), as.factor)

# fit the model
nsa_fit_k <- mmrm(
  formula = XNSACHG ~ ARM + AVISIT + ARM:AVISIT + XNSASCRBL + XNSASCRBL:AVISIT + us(AVISIT|SUBJID),
  data = nsa_mmrm_k,
  method = "Kenward-Roger")

# lsmean by ARM
nsa_lsm_value_k <- emmeans(nsa_fit_k, ~ ARM|AVISIT)
lsm_nsa_k <- summary(nsa_lsm_value_k) %>% filter(AVISIT == 60)
lsm_nsa_k <- flextable(lsm_nsa_k)

# lsmean difference
nsa_fit.rg_k <- ref_grid(nsa_fit_k)
lsmean_nsa_k <- pairs(emmeans(nsa_fit.rg_k, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_nsa_k <- confint(lsmean_nsa_k) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(contrast, lower.CL, upper.CL) %>% 
  filter(str_detect(contrast, "Placebo"))

df.lsm_nsa_k <- summary(lsmean_nsa_k) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  left_join(ci_lsm_nsa_k, by = c("contrast"))

# generate tables by visit
lsm_table_nsa_k <- flextable(df.lsm_nsa_k)


# Japan
tb_nsa_w6_j <- nsa_w6 %>%
  filter(REGION == 1) %>%
  select(ARM, XNSASCR, XNSACHG) %>% 
  tbl_strata(strata = c(ARM),
             .tbl_fun =
               ~.x %>%
               tbl_summary(type = list(where(is.numeric) ~ "continuous2"),
                           missing = "no",
                           statistic = list(all_continuous() ~ c("{mean} ± ({sd})",
                                                                 "{median}",
                                                                 "({min}, {max})")),
                           digits = everything() ~ 1) %>%
               add_n(statistic = "{N_nonmiss}",
                     col_label = "**N responses**")%>%
               modify_caption("**NSA-16 Total Score at Week 6(MMRM)**"))

tb_nsa_w6_j <- as_flex_table(tb_nsa_w6_j)

# MMRM
nsa_mmrm_j <- nsa %>%
  filter(FASFL == "Y") %>%
  filter(REGION == 1) %>%
  group_by(SUBJID) %>% 
  filter(XNSACD == 18) %>% 
  filter(AVISIT != -1) %>% 
  filter(!(is.na(XNSASCR))) %>%
  select(SUBJID, ARM, AVISIT, XNSASCRBL, XNSACHG) %>%
  mutate_at(vars(AVISIT, ARM), as.factor)

# fit the model
nsa_fit_j <- mmrm(
  formula = XNSACHG ~ ARM + AVISIT + ARM:AVISIT + XNSASCRBL + XNSASCRBL:AVISIT + us(AVISIT|SUBJID),
  data = nsa_mmrm_j,
  method = "Kenward-Roger")

# lsmean by ARM
nsa_lsm_value_j <- emmeans(nsa_fit_j, ~ ARM|AVISIT)
lsm_nsa_j <- summary(nsa_lsm_value_j) %>% filter(AVISIT == 60)
lsm_nsa_j <- flextable(lsm_nsa_j)

# lsmean difference
nsa_fit.rg_j <- ref_grid(nsa_fit_j)
lsmean_nsa_j <- pairs(emmeans(nsa_fit.rg_j, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_nsa_j <- confint(lsmean_nsa_j) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(contrast, lower.CL, upper.CL) %>% 
  filter(str_detect(contrast, "Placebo"))

df.lsm_nsa_j <- summary(lsmean_nsa_j) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  left_join(ci_lsm_nsa_j, by = c("contrast"))

# generate tables by visit
lsm_table_nsa_j <- flextable(df.lsm_nsa_j)

##############
# 4. NSA-16 global negative symptoms rating (change from baseline at w6)
nsa_neg_w6 <- nsa %>%
  filter(FASFL == "Y") %>%
  filter(XNSACD == 17) %>% 
  filter(AVISIT == 60)

tb_nsa_neg_w6 <- nsa_neg_w6 %>%
  select(ARM, XNSASCR, XNSACHG) %>% 
  tbl_strata(strata = c(ARM),
             .tbl_fun =
               ~.x %>%
               tbl_summary(type = list(where(is.numeric) ~ "continuous2"),
                           missing = "no",
                           statistic = list(all_continuous() ~ c("{mean} ± ({sd})",
                                                                 "{median}",
                                                                 "({min}, {max})")),
                           digits = everything() ~ 1) %>%
               add_n(statistic = "{N_nonmiss}",
                     col_label = "**N responses**")%>%
               modify_caption("**NSA-16 Global negative symptoms rating at Week 6(MMRM)**"))

tb_nsa_neg_w6 <- as_flex_table(tb_nsa_neg_w6)

# MMRM
nsa_mmrm_neg <- nsa %>%
  filter(FASFL == "Y") %>%
  group_by(SUBJID) %>% 
  filter(XNSACD == 17) %>% 
  filter(AVISIT != -1) %>% 
  filter(!(is.na(XNSASCR))) %>%
  select(SUBJID, ARM, REGION, AVISIT, XNSASCRBL, XNSACHG) %>%
  mutate_at(vars(REGION, AVISIT, ARM), as.factor)

# fit the model
nsa_fit_neg <- mmrm(
  formula = XNSACHG ~ ARM + REGION + AVISIT + ARM:AVISIT + XNSASCRBL + XNSASCRBL:AVISIT + us(AVISIT|SUBJID),
  data = nsa_mmrm_neg,
  method = "Kenward-Roger")

# lsmean by ARM
nsa_lsm_value_neg <- emmeans(nsa_fit_neg, ~ ARM|AVISIT)
lsm_nsa_neg <- summary(nsa_lsm_value_neg) %>% filter(AVISIT == 60)
lsm_nsa_neg <- flextable(lsm_nsa_neg)

# lsmean difference
nsa_fit.rg_neg <- ref_grid(nsa_fit_neg)
lsmean_nsa_neg <- pairs(emmeans(nsa_fit.rg_neg, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_nsa_neg <- confint(lsmean_nsa_neg) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(contrast, lower.CL, upper.CL) %>% 
  filter(str_detect(contrast, "Placebo"))

df.lsm_nsa_neg <- summary(lsmean_nsa_neg) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  left_join(ci_lsm_nsa_neg, by = c("contrast"))

# generate tables by visit
lsm_table_nsa_neg <- flextable(df.lsm_nsa_neg)

# KOREA
tb_nsa_w6_neg_k <- nsa_neg_w6 %>%
  filter(REGION == 2) %>%
  select(ARM, XNSASCR, XNSACHG) %>% 
  tbl_strata(strata = c(ARM),
             .tbl_fun =
               ~.x %>%
               tbl_summary(type = list(where(is.numeric) ~ "continuous2"),
                           missing = "no",
                           statistic = list(all_continuous() ~ c("{mean} ± ({sd})",
                                                                 "{median}",
                                                                 "({min}, {max})")),
                           digits = everything() ~ 1) %>%
               add_n(statistic = "{N_nonmiss}",
                     col_label = "**N responses**")%>%
               modify_caption("**NSA-16 Global negative symptoms rating at Week 6(MMRM)**"))

tb_nsa_w6_neg_k <- as_flex_table(tb_nsa_w6_neg_k)

# MMRM
nsa_mmrm_k_neg <- nsa %>%
  filter(FASFL == "Y") %>%
  filter(REGION == 2) %>%
  group_by(SUBJID) %>% 
  filter(XNSACD == 17) %>% 
  filter(AVISIT != -1) %>% 
  filter(!(is.na(XNSASCR))) %>%
  select(SUBJID, ARM, AVISIT, XNSASCRBL, XNSACHG) %>%
  mutate_at(vars(AVISIT, ARM), as.factor)

# fit the model
nsa_fit_k_neg <- mmrm(
  formula = XNSACHG ~ ARM + AVISIT + ARM:AVISIT + XNSASCRBL + XNSASCRBL:AVISIT + us(AVISIT|SUBJID),
  data = nsa_mmrm_k_neg,
  method = "Kenward-Roger")

# lsmean by ARM
nsa_lsm_value_k_neg <- emmeans(nsa_fit_k_neg, ~ ARM|AVISIT)
lsm_nsa_k_neg <- summary(nsa_lsm_value_k_neg) %>% filter(AVISIT == 60)
lsm_nsa_k_neg <- flextable(lsm_nsa_k_neg)

# lsmean difference
nsa_fit.rg_k_neg <- ref_grid(nsa_fit_k_neg)
lsmean_nsa_k_neg <- pairs(emmeans(nsa_fit.rg_k_neg, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_nsa_k_neg <- confint(lsmean_nsa_k_neg) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(contrast, lower.CL, upper.CL) %>% 
  filter(str_detect(contrast, "Placebo"))

df.lsm_nsa_k_neg <- summary(lsmean_nsa_k_neg) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  left_join(ci_lsm_nsa_k_neg, by = c("contrast"))

# generate tables by visit
lsm_table_nsa_k_neg <- flextable(df.lsm_nsa_k_neg)

# Japan
tb_nsa_w6_j_neg <- nsa_neg_w6 %>%
  filter(REGION == 1) %>%
  select(ARM, XNSASCR, XNSACHG) %>% 
  tbl_strata(strata = c(ARM),
             .tbl_fun =
               ~.x %>%
               tbl_summary(type = list(where(is.numeric) ~ "continuous2"),
                           missing = "no",
                           statistic = list(all_continuous() ~ c("{mean} ± ({sd})",
                                                                 "{median}",
                                                                 "({min}, {max})")),
                           digits = everything() ~ 1) %>%
               add_n(statistic = "{N_nonmiss}",
                     col_label = "**N responses**")%>%
               modify_caption("**NSA-16 Global negative symptoms rating at Week 6(MMRM)**"))

tb_nsa_w6_j_neg <- as_flex_table(tb_nsa_w6_j_neg)

# MMRM
nsa_mmrm_j_neg <- nsa %>%
  filter(FASFL == "Y") %>%
  filter(REGION == 1) %>%
  group_by(SUBJID) %>% 
  filter(XNSACD == 17) %>% 
  filter(AVISIT != -1) %>% 
  filter(!(is.na(XNSASCR))) %>%
  select(SUBJID, ARM, AVISIT, XNSASCRBL, XNSACHG) %>%
  mutate_at(vars(AVISIT, ARM), as.factor)

# fit the model
nsa_fit_j_neg <- mmrm(
  formula = XNSACHG ~ ARM + AVISIT + ARM:AVISIT + XNSASCRBL + XNSASCRBL:AVISIT + us(AVISIT|SUBJID),
  data = nsa_mmrm_j_neg,
  method = "Kenward-Roger")

# lsmean by ARM
nsa_lsm_value_j_neg <- emmeans(nsa_fit_j_neg, ~ ARM|AVISIT)
lsm_nsa_j_neg <- summary(nsa_lsm_value_j_neg) %>% filter(AVISIT == 60)
lsm_nsa_j_neg <- flextable(lsm_nsa_j_neg)

# lsmean difference
nsa_fit.rg_j_neg <- ref_grid(nsa_fit_j_neg)
lsmean_nsa_j_neg <- pairs(emmeans(nsa_fit.rg_j_neg, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_nsa_j_neg <- confint(lsmean_nsa_j_neg) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(contrast, lower.CL, upper.CL) %>% 
  filter(str_detect(contrast, "Placebo"))

df.lsm_nsa_j_neg <- summary(lsmean_nsa_j_neg) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  left_join(ci_lsm_nsa_j_neg, by = c("contrast"))

# generate tables by visit
lsm_table_nsa_j_neg <- flextable(df.lsm_nsa_j_neg)

# 5. PANSS positive score (change from baseline at w6)
# summary of positive score at week 6
whole_pos <- panss %>% 
  filter(FASFL == "Y") %>%
  filter(XPNSID == 31) 

panss_pos <- whole_pos %>%
  filter(AVISIT == 60) %>% 
  select(ARM, XPNSSCRE, XPNSCHG) %>% 
  tbl_summary(by = ARM,
              type = list(where(is.numeric) ~ "continuous2"),
              missing = "no",
              statistic = list(all_continuous() ~ c("{mean} ± ({sd})",
                                                    "{median}",
                                                    "({min}, {max})")),
              digits = everything() ~ 1)
panss_pos <- as_flex_table(panss_pos)

# MMRM
panss_mmrm_pos <- panss %>%
  filter(FASFL == "Y") %>%
  group_by(SUBJID) %>% 
  filter(XPNSID == 31) %>% 
  filter(AVISIT != -1) %>% 
  filter(!(is.na(XPNSSCRE))) %>%
  select(SUBJID, ARM, REGION, AVISIT, XPNSSCRBL, XPNSCHG) %>%
  mutate_at(vars(REGION, AVISIT, ARM), as.factor)

# fit the model
panss_fit_pos <- mmrm(
  formula = XPNSCHG ~ ARM + REGION + AVISIT + ARM:AVISIT + XPNSSCRBL + XPNSSCRBL:AVISIT + us(AVISIT|SUBJID),
  data = panss_mmrm_pos,
  method = "Kenward-Roger")

# lsmean by ARM
panss_lsm_value_pos <- emmeans(panss_fit_pos, ~ ARM|AVISIT)
lsm_panss_pos <- summary(panss_lsm_value_pos) %>% filter(AVISIT == 60)
lsm_panss_pos <- flextable(lsm_panss_pos)

# lsmean difference
panss_fit.rg_pos <- ref_grid(panss_fit_pos)
lsmean_panss_pos <- pairs(emmeans(panss_fit.rg_pos, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_panss_pos <- confint(lsmean_panss_pos) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(contrast, lower.CL, upper.CL) %>% 
  filter(str_detect(contrast, "Placebo"))

df.lsm_panss_pos <- summary(lsmean_panss_pos) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  left_join(ci_lsm_panss_pos, by = c("contrast"))

# generate a table
lsm_table_panss_pos <- flextable(df.lsm_panss_pos)

#### Korea
# summary table
panss_pos_k <- whole_pos %>%
  filter(REGION == 2) %>%
  filter(AVISIT == 60) %>% 
  select(ARM, XPNSSCRE, XPNSCHG) %>% 
  tbl_summary(by = ARM,
              type = list(where(is.numeric) ~ "continuous2"),
              missing = "no",
              statistic = list(all_continuous() ~ c("{mean} ± ({sd})",
                                                    "{median}",
                                                    "({min}, {max})")),
              digits = everything() ~ 1)
panss_pos_k <- as_flex_table(panss_pos_k)

# MMRM
panss_mmrm_pos_k <- panss %>%
  filter(FASFL == "Y") %>%
  filter(REGION == 2) %>%
  group_by(SUBJID) %>% 
  filter(XPNSID == 31) %>% 
  filter(AVISIT != -1) %>% 
  filter(!(is.na(XPNSSCRE))) %>%
  select(SUBJID, ARM, AVISIT, XPNSSCRBL, XPNSCHG) %>%
  mutate_at(vars(AVISIT, ARM), as.factor)

# fit the model
panss_fit_pos_k <- mmrm(
  formula = XPNSCHG ~ ARM + AVISIT + ARM:AVISIT + XPNSSCRBL + XPNSSCRBL:AVISIT + us(AVISIT|SUBJID),
  data = panss_mmrm_pos_k,
  method = "Kenward-Roger")

# lsmean by ARM
panss_lsm_value_pos_k <- emmeans(panss_fit_pos_k, ~ ARM|AVISIT)
lsm_panss_pos_k <- summary(panss_lsm_value_pos_k) %>% filter(AVISIT == 60)
lsm_panss_pos_k <- flextable(lsm_panss_pos_k)

# lsmean difference
panss_fit.rg_pos_k <- ref_grid(panss_fit_pos_k)
lsmean_panss_pos_k <- pairs(emmeans(panss_fit.rg_pos_k, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_panss_pos_k <- confint(lsmean_panss_pos_k) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(contrast, lower.CL, upper.CL) %>% 
  filter(str_detect(contrast, "Placebo"))

df.lsm_panss_pos_k <- summary(lsmean_panss_pos_k) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  left_join(ci_lsm_panss_pos_k, by = c("contrast"))

# generate a table
lsm_table_panss_pos_k <- flextable(df.lsm_panss_pos_k)

#### Japan
# summary table
panss_pos_j <- whole_pos %>%
  filter(REGION == 1) %>%
  filter(AVISIT == 60) %>% 
  select(ARM, XPNSSCRE, XPNSCHG) %>% 
  tbl_summary(by = ARM,
              type = list(where(is.numeric) ~ "continuous2"),
              missing = "no",
              statistic = list(all_continuous() ~ c("{mean} ± ({sd})",
                                                    "{median}",
                                                    "({min}, {max})")),
              digits = everything() ~ 1)
panss_pos_j <- as_flex_table(panss_pos_j)

# MMRM
panss_mmrm_pos_j <- panss %>%
  filter(FASFL == "Y") %>%
  filter(REGION == 1) %>%
  group_by(SUBJID) %>% 
  filter(XPNSID == 31) %>% 
  filter(AVISIT != -1) %>% 
  filter(!(is.na(XPNSSCRE))) %>%
  select(SUBJID, ARM, AVISIT, XPNSSCRBL, XPNSCHG) %>%
  mutate_at(vars(AVISIT, ARM), as.factor)

# fit the model
panss_fit_pos_j <- mmrm(
  formula = XPNSCHG ~ ARM + AVISIT + ARM:AVISIT + XPNSSCRBL + XPNSSCRBL:AVISIT + us(AVISIT|SUBJID),
  data = panss_mmrm_pos_j,
  method = "Kenward-Roger")

# lsmean by ARM
panss_lsm_value_pos_j <- emmeans(panss_fit_pos_j, ~ ARM|AVISIT)
lsm_panss_pos_j <- summary(panss_lsm_value_pos_j) %>% filter(AVISIT == 60)
lsm_panss_pos_j <- flextable(lsm_panss_pos_j)

# lsmean difference
panss_fit.rg_pos_j <- ref_grid(panss_fit_pos_j)
lsmean_panss_pos_j <- pairs(emmeans(panss_fit.rg_pos_j, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_panss_pos_j <- confint(lsmean_panss_pos_j) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(contrast, lower.CL, upper.CL) %>% 
  filter(str_detect(contrast, "Placebo"))

df.lsm_panss_pos_j <- summary(lsmean_panss_pos_j) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  left_join(ci_lsm_panss_pos_j, by = c("contrast"))

# generate a table
lsm_table_panss_pos_j <- flextable(df.lsm_panss_pos_j)

# 6. PANSS negative score (change from baseline at w6) 
whole_neg <- panss %>% 
  filter(FASFL == "Y") %>%
  filter(XPNSID == 32) 

panss_neg <- whole_neg %>%
  filter(AVISIT == 60) %>% 
  select(ARM, XPNSSCRE, XPNSCHG) %>% 
  tbl_summary(by = ARM,
              type = list(where(is.numeric) ~ "continuous2"),
              missing = "no",
              statistic = list(all_continuous() ~ c("{mean} ± ({sd})",
                                                    "{median}",
                                                    "({min}, {max})")),
              digits = everything() ~ 1)
panss_neg <- as_flex_table(panss_neg)

# MMRM
panss_mmrm_neg <- panss %>%
  filter(FASFL == "Y") %>%
  group_by(SUBJID) %>% 
  filter(XPNSID == 32) %>% 
  filter(AVISIT != -1) %>% 
  filter(!(is.na(XPNSSCRE))) %>%
  select(SUBJID, ARM, REGION, AVISIT, XPNSSCRBL, XPNSCHG) %>%
  mutate_at(vars(REGION, AVISIT, ARM), as.factor)

# fit the model
panss_fit_neg <- mmrm(
  formula = XPNSCHG ~ ARM + REGION + AVISIT + ARM:AVISIT + XPNSSCRBL + XPNSSCRBL:AVISIT + us(AVISIT|SUBJID),
  data = panss_mmrm_neg,
  method = "Kenward-Roger")

# lsmean by ARM
panss_lsm_value_neg <- emmeans(panss_fit_neg, ~ ARM|AVISIT)
lsm_panss_neg <- summary(panss_lsm_value_neg) %>% filter(AVISIT == 60)
lsm_panss_neg <- flextable(lsm_panss_neg)

# lsmean difference
panss_fit.rg_neg <- ref_grid(panss_fit_neg)
lsmean_panss_neg <- pairs(emmeans(panss_fit.rg_neg, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_panss_neg <- confint(lsmean_panss_neg) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(contrast, lower.CL, upper.CL) %>% 
  filter(str_detect(contrast, "Placebo"))

df.lsm_panss_neg <- summary(lsmean_panss_neg) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  left_join(ci_lsm_panss_neg, by = c("contrast"))

# generate a table
lsm_table_panss_neg <- flextable(df.lsm_panss_neg)

#### Korea
# summary table
panss_neg_k <- whole_neg %>%
  filter(REGION == 2) %>%
  filter(AVISIT == 60) %>% 
  select(ARM, XPNSSCRE, XPNSCHG) %>% 
  tbl_summary(by = ARM,
              type = list(where(is.numeric) ~ "continuous2"),
              missing = "no",
              statistic = list(all_continuous() ~ c("{mean} ± ({sd})",
                                                    "{median}",
                                                    "({min}, {max})")),
              digits = everything() ~ 1)
panss_neg_k <- as_flex_table(panss_neg_k)

# MMRM
panss_mmrm_neg_k <- panss %>%
  filter(FASFL == "Y") %>%
  filter(REGION == 2) %>%
  group_by(SUBJID) %>% 
  filter(XPNSID == 32) %>% 
  filter(AVISIT != -1) %>% 
  filter(!(is.na(XPNSSCRE))) %>%
  select(SUBJID, ARM, AVISIT, XPNSSCRBL, XPNSCHG) %>%
  mutate_at(vars(AVISIT, ARM), as.factor)

# fit the model
panss_fit_neg_k <- mmrm(
  formula = XPNSCHG ~ ARM + AVISIT + ARM:AVISIT + XPNSSCRBL + XPNSSCRBL:AVISIT + us(AVISIT|SUBJID),
  data = panss_mmrm_neg_k,
  method = "Kenward-Roger")

# lsmean by ARM
panss_lsm_value_neg_k <- emmeans(panss_fit_neg_k, ~ ARM|AVISIT)
lsm_panss_neg_k <- summary(panss_lsm_value_neg_k) %>% filter(AVISIT == 60)
lsm_panss_neg_k <- flextable(lsm_panss_neg_k)

# lsmean difference
panss_fit.rg_neg_k <- ref_grid(panss_fit_neg_k)
lsmean_panss_neg_k <- pairs(emmeans(panss_fit.rg_neg_k, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_panss_neg_k <- confint(lsmean_panss_neg_k) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(contrast, lower.CL, upper.CL) %>% 
  filter(str_detect(contrast, "Placebo"))

df.lsm_panss_neg_k <- summary(lsmean_panss_neg_k) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  left_join(ci_lsm_panss_neg_k, by = c("contrast"))

# generate a table
lsm_table_panss_neg_k <- flextable(df.lsm_panss_neg_k)

#### Japan
# summary table
panss_neg_j <- whole_neg %>%
  filter(REGION == 1) %>%
  filter(AVISIT == 60) %>% 
  select(ARM, XPNSSCRE, XPNSCHG) %>% 
  tbl_summary(by = ARM,
              type = list(where(is.numeric) ~ "continuous2"),
              missing = "no",
              statistic = list(all_continuous() ~ c("{mean} ± ({sd})",
                                                    "{median}",
                                                    "({min}, {max})")),
              digits = everything() ~ 1)
panss_neg_j <- as_flex_table(panss_neg_j)

# MMRM
panss_mmrm_neg_j <- panss %>%
  filter(FASFL == "Y") %>%
  filter(REGION == 1) %>%
  group_by(SUBJID) %>% 
  filter(XPNSID == 32) %>% 
  filter(AVISIT != -1) %>% 
  filter(!(is.na(XPNSSCRE))) %>%
  select(SUBJID, ARM, AVISIT, XPNSSCRBL, XPNSCHG) %>%
  mutate_at(vars(AVISIT, ARM), as.factor)

# fit the model
panss_fit_neg_j <- mmrm(
  formula = XPNSCHG ~ ARM + AVISIT + ARM:AVISIT + XPNSSCRBL + XPNSSCRBL:AVISIT + us(AVISIT|SUBJID),
  data = panss_mmrm_neg_j,
  method = "Kenward-Roger")

# lsmean by ARM
panss_lsm_value_neg_j <- emmeans(panss_fit_neg_j, ~ ARM|AVISIT)
lsm_panss_neg_j <- summary(panss_lsm_value_neg_j) %>% filter(AVISIT == 60)
lsm_panss_neg_j <- flextable(lsm_panss_neg_j)

# lsmean difference
panss_fit.rg_neg_j <- ref_grid(panss_fit_neg_j)
lsmean_panss_neg_j <- pairs(emmeans(panss_fit.rg_neg_j, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_panss_neg_j <- confint(lsmean_panss_neg_j) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(contrast, lower.CL, upper.CL) %>% 
  filter(str_detect(contrast, "Placebo"))

df.lsm_panss_neg_j <- summary(lsmean_panss_neg_j) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  left_join(ci_lsm_panss_neg_j, by = c("contrast"))

# generate a table
lsm_table_panss_neg_j <- flextable(df.lsm_panss_neg_j)

#------------------------ for importing as doc
word_export <- read_docx()
word_export <- body_add_par(word_export, "PANSS total score by visit time", style = "heading 3")
word_export <- body_add_flextable(word_export, entire_panss)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsmean_panss)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table)
word_export <- body_add_break(word_export)
word_export <- body_add_par(word_export, "PANSS total score change at week 6 - Korea", style = "heading 2")
word_export <- body_add_flextable(word_export, korea_panss)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsmean_panss_k)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_k)
word_export <- body_add_break(word_export)
word_export <- body_add_par(word_export, "PANSS total score change at week 6 - Japan", style = "heading 2")
word_export <- body_add_flextable(word_export, japan_panss)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsmean_panss_j)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_j)
word_export <- body_add_break(word_export)

word_export <- body_add_par(word_export, "CGI-I score at week 6", style = "heading 3")
word_export <- body_add_flextable(word_export, tb_cgi_w6)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsmean_cgi)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_cgi)
word_export <- body_add_break(word_export)
word_export <- body_add_par(word_export, "CGI-I score at week 6 - Korea", style = "heading 2")
word_export <- body_add_flextable(word_export, tb_cgi_w6_k)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsmean_cgi_k)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_cgi_k)
word_export <- body_add_break(word_export)
word_export <- body_add_par(word_export, "CGI-I score at week 6 - Japan", style = "heading 2")
word_export <- body_add_flextable(word_export, tb_cgi_w6_j)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsmean_cgi_j)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_cgi_j)
word_export <- body_add_break(word_export)

word_export <- body_add_par(word_export, "NSA-16 total score change at week 6", style = "heading 3")
word_export <- body_add_flextable(word_export, tb_nsa_w6)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_nsa)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_nsa)
word_export <- body_add_break(word_export)
word_export <- body_add_par(word_export, "NSA-16 total score change at week 6 - Korea", style = "heading 2")
word_export <- body_add_flextable(word_export, tb_nsa_w6_k)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_nsa_k)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_nsa_k)
word_export <- body_add_break(word_export)
word_export <- body_add_par(word_export, "NSA-16 total score change at week 6 - Japan", style = "heading 2")
word_export <- body_add_flextable(word_export, tb_nsa_w6_j)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_nsa_j)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_nsa_j)
word_export <- body_add_break(word_export)

word_export <- body_add_par(word_export, "NSA-16 negative change at week 6", style = "heading 3")
word_export <- body_add_flextable(word_export, tb_nsa_neg_w6)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_nsa_neg)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_nsa_neg)
word_export <- body_add_break(word_export)
word_export <- body_add_par(word_export, "NSA-16 negative change at week 6 - Korea", style = "heading 2")
word_export <- body_add_flextable(word_export, tb_nsa_w6_neg_k)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_nsa_k_neg)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_nsa_k_neg)
word_export <- body_add_break(word_export)
word_export <- body_add_par(word_export, "NSA-16 negative change at week 6 - Japan", style = "heading 2")
word_export <- body_add_flextable(word_export, tb_nsa_w6_j_neg)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_nsa_j_neg)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_nsa_j_neg)
word_export <- body_add_break(word_export)

word_export <- body_add_par(word_export, "PANSS positive change at week 6", style = "heading 3")
word_export <- body_add_flextable(word_export, panss_pos)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_panss_pos)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_panss_pos)
word_export <- body_add_break(word_export)
word_export <- body_add_par(word_export, "PANSS positive change at week 6 - Korea", style = "heading 2")
word_export <- body_add_flextable(word_export, panss_pos_k)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_panss_pos_k)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_panss_pos_k)
word_export <- body_add_break(word_export)
word_export <- body_add_par(word_export, "PANSS positive change at week 6 - Japan", style = "heading 2")
word_export <- body_add_flextable(word_export, panss_pos_j)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_panss_pos_j)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_panss_pos_j)
word_export <- body_add_break(word_export)

word_export <- body_add_par(word_export, "PANSS negative change at week 6", style = "heading 3")
word_export <- body_add_flextable(word_export, panss_neg)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_panss_neg)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_panss_neg)
word_export <- body_add_break(word_export)
word_export <- body_add_par(word_export, "PANSS negative change at week 6 - Korea", style = "heading 2")
word_export <- body_add_flextable(word_export, panss_neg_k)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_panss_neg_k)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_panss_neg_k)
word_export <- body_add_break(word_export)
word_export <- body_add_par(word_export, "PANSS negative change at week 6 - Japan", style = "heading 2")
word_export <- body_add_flextable(word_export, panss_neg_j)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_panss_neg_j)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_panss_neg_j)
word_export <- body_add_break(word_export)

print(word_export, 'A002-A4_tables.docx')