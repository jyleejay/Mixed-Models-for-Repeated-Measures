library(tidyverse)
library(sas7bdat)
library(officer)
library(reshape2)
library(haven)
library(gtsummary)
library(mmrm)
library(emmeans)
library(kableExtra)
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

# Load CSV Files
data.path <- paste0(getwd(), "/A002-A4dataset")
data.files <- list.files(data.path)

# 1. PANSS total score (change from baseline by visit)
panss <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<adpanss\\>", data.files)]), NULL))
# move placebo into first
panss$ARM <- fct_relevel(panss$ARM, "Placebo")

# whether investigator training   - from additional SAP p.4
# Site ID Site Name Investigator Meeting Attendance
# 2001 Inha University Hospital Y
# 2002 Gachon University Gil Hospital Y
# 2003 Chonbuk National University Hospital Y
# 2004 Eulji General Hospital Y
# 2005 Inje University Haeundae Paik Hospital Y
# 2006 The Catholic University of Korea Yeouido St.Mary's Hospital Y
# 2007 Severance Hospital, Yonsei University Health System Y
# 2008 Chungnam National University Hospital Y
# 2009 Korea University Ansan Hospital Y
# 2010 Konkuk University Medical Center Y
# 2011 Hallym University Hangang Sacred Heart Hospital Y
# 2012 Kyungpook National University Hospital Y
# 2013 Jeju national university hospital Y
# 2014 Kyung Hee University Hospital Y
# 2015 DongGuk University Gyongju Hospital Y
# 2016 National Health Insurance Service Ilsan Hospital N
# 2017 Seoul National University Bundang Hospital N
# 2018 The Catholic University of Korea St. Vincent's Hospital N
# 2019 Hanyang University Guri Hospital N
# 2020 Inje University Busan Paik Hospital N
# 2021 Chonnam National University Hospital N
# 2022 Yong-In Mental Hospital N
# 2023 Naju National Hospital N
# 2024 Pusan National University Yangsan Hospital N
# 2025 Hallym Univerisity Sacred Heart Hospital N

# whether training
adsl <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<adsl\\>", data.files)]), NULL))

# move placebo into first
adsl$ARM <- fct_relevel(adsl$ARM, "Placebo")

# site numbers
tr_yes <- c(2001:2015)
tr_no <- c(2016:2025)

tr_var <- adsl %>% 
  filter(REGION == 2) %>% 
  mutate(inv_tr = ifelse(SITEID %in% tr_no, "N",
                         ifelse(SITEID %in% tr_yes, "Y", 0))) %>% 
  select(SUBJID, inv_tr)

# Korea
training <- panss %>%
  filter(FASFL == "Y") %>%
  filter(REGION == 2) %>%
  select(SUBJID, ARM, REGION, AVISIT, XPNSID, XPNSSCRBL, XPNSSCRE, XPNSCHG) %>% 
  left_join(tr_var, by = "SUBJID")

count <- training %>% 
  filter(XPNSID == 33) %>% 
  filter(AVISIT == 60) %>% 
  group_by(SUBJID)%>%
  select(SUBJID, ARM, inv_tr) %>%
  unique() %>% 
  mutate(inv_tr = factor(inv_tr, levels = c("Y", "N"), labels = c("Yes", "No"))) %>% 
  tbl_cross(row = ARM, col = inv_tr,
            label = list(inv_tr ~ "**PANSS Investigator Training in Korea**"))
count <- as_flex_table(count)

# japan, Taiwan
other <- panss %>%
  filter(FASFL == "Y") %>%
  filter(REGION != 2) %>% 
  select(SUBJID, ARM, REGION, AVISIT, XPNSID, XPNSSCRBL, XPNSSCRE, XPNSCHG) %>% 
  mutate(inv_tr = "Y")

# make a entire dataset
tr_df <- training %>%
  bind_rows(other) %>% 
  mutate_at(vars(REGION, AVISIT, ARM, inv_tr), as.factor)

tb_inv_tr <- tr_df %>%
  filter(XPNSID == 33) %>%
  filter(AVISIT == 60) %>% 
  select(inv_tr, REGION, ARM) %>% 
  mutate(inv_tr = factor(inv_tr, levels = c("Y", "N"), labels = c("Yes", "No")),
         REGION = factor(REGION, levels = c(1, 2, 3), labels = c("Japan", "Korea", "Taiwan"))) %>%
  tbl_strata(strata = c(inv_tr),
             .tbl_fun =
               ~.x %>%
               tbl_summary(by = REGION,
                           digits = everything() ~ 1) %>%
               modify_caption("**PANSS Investigator Training**") %>% 
               add_overall())
tb_inv_tr <- as_flex_table(tb_inv_tr)

# get rid of not trained ones
tr_y_df <- tr_df %>% 
  filter(inv_tr == "Y") %>%
  filter(XPNSID == 33) %>%
  filter(AVISIT == 60) %>% 
  mutate(REGION = as.factor(ifelse(REGION !=2, 1,
                         ifelse(REGION == 2, 2, 0))))

tb_tr_y <- tr_y_df %>%
  select(REGION, ARM) %>% 
  mutate(REGION = factor(REGION, levels = c(1, 2), labels = c("Japan or Taiwan", "Korea"))) %>%
  tbl_summary(by = REGION,
              digits = everything() ~ 1) %>%
  modify_caption("**Investigator trained only**") %>% 
  add_overall()

tb_tr_y <- as_flex_table(tb_tr_y)

#---------------------------------------
# panss total change at W6
# Korea trained only
panss_mmrm <- tr_df %>%
  filter(inv_tr == "Y") %>% 
  filter(AVISIT != -1) %>%
  filter(REGION == 2) %>% 
  filter(XPNSID == 33) %>%
  select(SUBJID, ARM, AVISIT, XPNSSCRBL, XPNSSCRE, XPNSCHG) %>%
  mutate_at(vars(AVISIT, ARM), as.factor)

tb_train <- panss_mmrm %>% 
  filter(AVISIT == 60) %>% 
  select(ARM, XPNSSCRE, XPNSCHG) %>% 
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
               modify_caption("**PANSS Score at Week 6(MMRM) - Korea trained**"))
tb_train <- as_flex_table(tb_train)

# fit the model
tr_fit <- mmrm(
  formula = XPNSCHG ~ ARM + AVISIT + ARM:AVISIT + XPNSSCRBL + XPNSSCRBL:AVISIT + us(AVISIT|SUBJID),
  data = panss_mmrm,
  method = "Kenward-Roger")

# lsmean by ARM
panss_lsm_value <- emmeans(tr_fit, ~ ARM|AVISIT)
lsm_panss <- summary(panss_lsm_value) %>% filter(AVISIT == 60)
lsm_panss <- flextable(lsm_panss)

# lsmean difference
panss_fit.rg <- ref_grid(tr_fit)
panss_lsm <- pairs(emmeans(panss_fit.rg, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm <- confint(panss_lsm) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(contrast, lower.CL, upper.CL) %>% 
  filter(str_detect(contrast, "Placebo"))

df.lsm <- summary(panss_lsm) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  left_join(ci_lsm, by = c("contrast"))

# generate a table
lsm_table <- flextable(df.lsm)


# Korea no-trained
panss_mmrm_no <- tr_df %>%
  filter(inv_tr == "N") %>% 
  filter(AVISIT != -1) %>%
  filter(REGION == 2) %>% 
  filter(XPNSID == 33) %>%
  select(SUBJID, ARM, AVISIT, REGION, XPNSSCRBL, XPNSSCRE, XPNSCHG) %>%
  mutate_at(vars(AVISIT, ARM), as.factor)

tb_notrain <- panss_mmrm_no %>% 
  filter(AVISIT == 60) %>% 
  select(ARM, XPNSSCRE, XPNSCHG) %>% 
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
               modify_caption("**PANSS Score at Week 6(MMRM) - Korea no trained**"))

tb_notrain <- as_flex_table(tb_notrain)

# fit the model
tr_fit_no <- mmrm(
  formula = XPNSCHG ~ ARM + AVISIT + ARM:AVISIT + XPNSSCRBL + XPNSSCRBL:AVISIT + us(AVISIT|SUBJID),
  data = panss_mmrm_no,
  method = "Kenward-Roger")

# lsmean by ARM
panss_lsm_no <- emmeans(tr_fit_no, ~ ARM|AVISIT)
lsm_panss_no <- summary(panss_lsm_no) %>%
  filter(AVISIT == 60) %>% 
  select(-c(df))

lsm_panss_no <- flextable(lsm_panss_no)
lsm_panss_no <- set_caption(lsm_panss_no, caption = "PANSS total score change at Week 6", 
                           style = "Table Caption")

# lsmean difference
panss_fit_no <- ref_grid(tr_fit_no)
panss_lsm_no <- pairs(emmeans(panss_fit_no, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_no <- confint(panss_lsm_no) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(contrast, lower.CL, upper.CL) %>% 
  filter(str_detect(contrast, "Placebo"))

df.lsm_no <- summary(panss_lsm_no) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  left_join(ci_lsm_no, by = c("contrast"))

# generate a table
lsm_table_no <- flextable(df.lsm_no)
lsm_table_no

######## whole study training only
tr_df_whole <- tr_df %>% 
  filter(inv_tr == "Y") %>% 
  filter(AVISIT != -1) %>% 
  filter(XPNSID == 33) %>% 
  select(SUBJID, ARM, REGION, AVISIT, XPNSSCRBL, XPNSCHG, inv_tr) %>%
  mutate_at(vars(REGION, AVISIT, ARM, inv_tr), as.factor)

# fit the model
tr_fit_w <- mmrm(
  formula = XPNSCHG ~ ARM + REGION + AVISIT + ARM:AVISIT + XPNSSCRBL + XPNSSCRBL:AVISIT + us(AVISIT|SUBJID),
  data = tr_df_whole,
  method = "Kenward-Roger")

# lsmean by ARM
panss_lsm_whole <- emmeans(tr_fit_w, ~ ARM|AVISIT)
lsm_panss_w <- summary(panss_lsm_whole) %>% filter(AVISIT == 60)
lsm_panss_w <- flextable(lsm_panss_w)

# lsmean difference
panss_fit_w <- ref_grid(tr_fit_w)
panss_lsm_w <- pairs(emmeans(panss_fit_w, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_w <- confint(panss_lsm_w) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(contrast, lower.CL, upper.CL) %>% 
  filter(str_detect(contrast, "Placebo"))

df.lsm_w <- summary(panss_lsm_w) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  left_join(ci_lsm_w, by = c("contrast"))

# generate a table
lsm_table_w <- flextable(df.lsm_w)

#----------------------------------
# panss positive score
# Korea
training_pos <- training %>%
  filter(XPNSID == 31) %>% 
  select(SUBJID, ARM, REGION, AVISIT, XPNSSCRBL, XPNSCHG) %>% 
  left_join(tr_var, by = "SUBJID")

# japan, Taiwan
other_pos <- panss %>%
  filter(FASFL == "Y") %>%
  filter(REGION != 2) %>% 
  filter(XPNSID == 31) %>%
  select(SUBJID, ARM, REGION, AVISIT, XPNSSCRBL, XPNSCHG) %>% 
  mutate(inv_tr = "Y")

# make a entire dataset
tr_df_pos <- training_pos %>%
  bind_rows(other_pos) %>% 
  mutate_at(vars(REGION, AVISIT, ARM, inv_tr), as.factor)

# Korea trained only
panss_mmrm_pos <- training_pos %>%
  filter(AVISIT != -1) %>% 
  select(SUBJID, ARM, AVISIT, XPNSSCRBL, XPNSCHG) %>%
  mutate_at(vars(AVISIT, ARM), as.factor)

# fit the model
tr_fit_pos <- mmrm(
  formula = XPNSCHG ~ ARM + AVISIT + ARM:AVISIT + XPNSSCRBL + XPNSSCRBL:AVISIT + us(AVISIT|SUBJID),
  data = panss_mmrm_pos,
  method = "Kenward-Roger")

# lsmean by ARM
panss_lsm_pos <- emmeans(tr_fit_pos, ~ ARM|AVISIT)
lsm_panss_pos <- summary(panss_lsm_pos) %>% filter(AVISIT == 60)
lsm_panss_pos <- flextable(lsm_panss_pos)

# lsmean difference
panss_fit.pos <- ref_grid(tr_fit_pos)
panss_lsm_pos <- pairs(emmeans(panss_fit.pos, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_pos <- confint(panss_lsm_pos) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(contrast, lower.CL, upper.CL) %>% 
  filter(str_detect(contrast, "Placebo"))

df.lsm_pos <- summary(panss_lsm_pos) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  left_join(ci_lsm_pos, by = c("contrast"))

# generate a table
lsm_table_pos <- flextable(df.lsm_pos)

# Korea no-trained
panss_pos_no <- training %>%
  filter(inv_tr == "N") %>% 
  filter(XPNSID == 31) %>% 
  filter(AVISIT != -1) %>% 
  select(SUBJID, ARM, AVISIT, XPNSSCRBL, XPNSCHG) %>%
  mutate_at(vars(AVISIT, ARM), as.factor)

# fit the model
tr_fit_no_pos <- mmrm(
  formula = XPNSCHG ~ ARM + AVISIT + ARM:AVISIT + XPNSSCRBL + XPNSSCRBL:AVISIT + us(AVISIT|SUBJID),
  data = panss_pos_no,
  method = "Kenward-Roger")

# lsmean by ARM
panss_lsm_no_pos <- emmeans(tr_fit_no_pos, ~ ARM|AVISIT)
lsm_panss_no_pos <- summary(panss_lsm_no_pos) %>% filter(AVISIT == 60)
lsm_panss_no_pos <- flextable(lsm_panss_no_pos)

# lsmean difference
panss_fit_no_pos <- ref_grid(tr_fit_no_pos)
panss_lsm_no_pos <- pairs(emmeans(panss_fit_no_pos, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_no_pos <- confint(panss_lsm_no_pos) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(contrast, lower.CL, upper.CL) %>% 
  filter(str_detect(contrast, "Placebo"))

df.lsm_no_pos <- summary(panss_lsm_no_pos) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  left_join(ci_lsm_no_pos, by = c("contrast"))

# generate a table
lsm_table_no_pos <- flextable(df.lsm_no_pos)

######## whole study training only
tr_df_whole_pos <- tr_df_pos %>% 
  filter(inv_tr == "Y") %>% 
  filter(AVISIT != -1) %>% 
  select(SUBJID, ARM, REGION, AVISIT, XPNSSCRBL, XPNSCHG) %>%
  mutate_at(vars(REGION, AVISIT, ARM), as.factor)

# fit the model
tr_fit_w_pos <- mmrm(
  formula = XPNSCHG ~ ARM + REGION + AVISIT + ARM:AVISIT + XPNSSCRBL + XPNSSCRBL:AVISIT + us(AVISIT|SUBJID),
  data = tr_df_whole_pos,
  method = "Kenward-Roger")

# lsmean by ARM
panss_lsm_whole_pos <- emmeans(tr_fit_w_pos, ~ ARM|AVISIT)
lsm_panss_w_pos <- summary(panss_lsm_whole_pos) %>% filter(AVISIT == 60)
lsm_panss_w_pos <- flextable(lsm_panss_w_pos)

# lsmean difference
panss_fit_w_pos <- ref_grid(tr_fit_w_pos)
panss_lsm_w_pos <- pairs(emmeans(panss_fit_w_pos, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_w_pos <- confint(panss_lsm_w_pos) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(contrast, lower.CL, upper.CL) %>% 
  filter(str_detect(contrast, "Placebo"))

df.lsm_w_pos <- summary(panss_lsm_w_pos) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  left_join(ci_lsm_w_pos, by = c("contrast"))

# generate a table
lsm_table_w_pos <- flextable(df.lsm_w_pos)

# panss negative score
# Korea
training_neg <- training %>%
  filter(XPNSID == 32) %>% 
  select(SUBJID, ARM, REGION, AVISIT, XPNSSCRBL, XPNSCHG) %>% 
  left_join(tr_var, by = "SUBJID")

# japan, Taiwan
other_neg <- panss %>%
  filter(FASFL == "Y") %>%
  filter(REGION != 2) %>% 
  filter(XPNSID == 32) %>%
  select(SUBJID, ARM, REGION, AVISIT, XPNSSCRBL, XPNSCHG) %>% 
  mutate(inv_tr = "Y")

# make a entire dataset
tr_df_neg <- training_neg %>%
  bind_rows(other_neg) %>% 
  mutate_at(vars(REGION, AVISIT, ARM, inv_tr), as.factor)

# Korea trained only
panss_mmrm_neg <- training_neg %>%
  filter(AVISIT != -1) %>% 
  select(SUBJID, ARM, AVISIT, XPNSSCRBL, XPNSCHG) %>%
  mutate_at(vars(AVISIT, ARM), as.factor)

# fit the model
tr_fit_neg <- mmrm(
  formula = XPNSCHG ~ ARM + AVISIT + ARM:AVISIT + XPNSSCRBL + XPNSSCRBL:AVISIT + us(AVISIT|SUBJID),
  data = panss_mmrm_neg,
  method = "Kenward-Roger")

# lsmean by ARM
panss_lsm_neg <- emmeans(tr_fit_neg, ~ ARM|AVISIT)
lsm_panss_neg <- summary(panss_lsm_neg) %>% filter(AVISIT == 60)
lsm_panss_neg <- flextable(lsm_panss_neg)

# lsmean difference
panss_fit.pos <- ref_grid(tr_fit_neg)
panss_lsm_neg <- pairs(emmeans(panss_fit.pos, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_neg <- confint(panss_lsm_neg) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(contrast, lower.CL, upper.CL) %>% 
  filter(str_detect(contrast, "Placebo"))

df.lsm_neg <- summary(panss_lsm_neg) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  left_join(ci_lsm_neg, by = c("contrast"))

# generate a table
lsm_table_neg <- flextable(df.lsm_neg)

# Korea no-trained
panss_neg_no <- training %>%
  filter(inv_tr == "N") %>% 
  filter(XPNSID == 32) %>% 
  filter(AVISIT != -1) %>% 
  select(SUBJID, ARM, AVISIT, XPNSSCRBL, XPNSCHG) %>%
  mutate_at(vars(AVISIT, ARM), as.factor)

# fit the model
tr_fit_no_neg <- mmrm(
  formula = XPNSCHG ~ ARM + AVISIT + ARM:AVISIT + XPNSSCRBL + XPNSSCRBL:AVISIT + us(AVISIT|SUBJID),
  data = panss_neg_no,
  method = "Kenward-Roger")

# lsmean by ARM
panss_lsm_no_neg <- emmeans(tr_fit_no_neg, ~ ARM|AVISIT)
lsm_panss_no_neg <- summary(panss_lsm_no_neg) %>% filter(AVISIT == 60)
lsm_panss_no_neg <- flextable(lsm_panss_no_neg)

# lsmean difference
panss_fit_no_neg <- ref_grid(tr_fit_no_neg)
panss_lsm_no_neg <- pairs(emmeans(panss_fit_no_neg, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_no_neg <- confint(panss_lsm_no_neg) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(contrast, lower.CL, upper.CL) %>% 
  filter(str_detect(contrast, "Placebo"))

df.lsm_no_neg <- summary(panss_lsm_no_neg) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  left_join(ci_lsm_no_neg, by = c("contrast"))

# generate a table
lsm_table_no_neg <- flextable(df.lsm_no_neg)

######## whole study training only
# Korea
tr_df_whole_neg <- tr_df_neg %>% 
  filter(inv_tr == "Y") %>% 
  filter(AVISIT != -1) %>% 
  select(SUBJID, ARM, REGION, AVISIT, XPNSSCRBL, XPNSCHG) %>%
  mutate_at(vars(REGION, AVISIT, ARM), as.factor)

# fit the model
tr_fit_w_neg <- mmrm(
  formula = XPNSCHG ~ ARM + REGION + AVISIT + ARM:AVISIT + XPNSSCRBL + XPNSSCRBL:AVISIT + us(AVISIT|SUBJID),
  data = tr_df_whole_neg,
  method = "Kenward-Roger")

# lsmean by ARM
panss_lsm_whole_neg <- emmeans(tr_fit_w_neg, ~ ARM|AVISIT)
lsm_panss_w_neg <- summary(panss_lsm_whole_neg) %>% filter(AVISIT == 60)
lsm_panss_w_neg <- flextable(lsm_panss_w_neg)

# lsmean difference
panss_fit_w_neg <- ref_grid(tr_fit_w_neg)
panss_lsm_w_neg <- pairs(emmeans(panss_fit_w_neg, ~ ARM|AVISIT), adjust = "none", reverse = T)

# make a df for lsmean differ table
ci_lsm_w_neg <- confint(panss_lsm_w_neg) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(contrast, lower.CL, upper.CL) %>% 
  filter(str_detect(contrast, "Placebo"))

df.lsm_w_neg <- summary(panss_lsm_w_neg) %>%
  arrange(AVISIT, contrast) %>% 
  filter(AVISIT == 60) %>% 
  select(-c(t.ratio, SE, df)) %>% 
  filter(str_detect(contrast, "Placebo")) %>% 
  left_join(ci_lsm_w_neg, by = c("contrast"))

# generate a table
lsm_table_w_neg <- flextable(df.lsm_w_neg)

#------------------------ for importing as doc
word_export <- read_docx()
word_export <- body_add_par(word_export, "Whether PANSS investigator training", style = "heading 3")
word_export <- body_add_flextable(word_export, tb_inv_tr)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, tb_tr_y)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, count)
word_export <- body_add_break(word_export)
word_export <- body_add_par(word_export, "PANSS investigator - Trained - Korea", style = "heading 2")
word_export <- body_add_flextable(word_export, tb_train)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_panss)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table)
word_export <- body_add_break(word_export)
word_export <- body_add_par(word_export, "PANSS investigator - Not trained - Korea", style = "heading 2")
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, tb_notrain)
word_export <- body_add_flextable(word_export, lsm_panss_no)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_no)
word_export <- body_add_break(word_export)
word_export <- body_add_par(word_export, "PANSS investigator trained - Entire group", style = "heading 2")
word_export <- body_add_flextable(word_export, lsm_panss_w)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_w)
word_export <- body_add_break(word_export)
word_export <- body_add_par(word_export, "PANSS positive score - Trained - Korea", style = "heading 2")
word_export <- body_add_flextable(word_export, lsm_panss_pos)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_pos)
word_export <- body_add_break(word_export)
word_export <- body_add_par(word_export, "PANSS positive score - Not trained - Korea", style = "heading 2")
word_export <- body_add_flextable(word_export, lsm_panss_no_pos)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_no_pos)
word_export <- body_add_break(word_export)
word_export <- body_add_par(word_export, "PANSS positive score - trained - Entire group", style = "heading 2")
word_export <- body_add_flextable(word_export, lsm_panss_w_pos)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_w_pos)
word_export <- body_add_break(word_export)
word_export <- body_add_par(word_export, "PANSS negtive score - Trained - Korea", style = "heading 2")
word_export <- body_add_flextable(word_export, lsm_panss_neg)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_neg)
word_export <- body_add_break(word_export)
word_export <- body_add_par(word_export, "PANSS negative score - Not trained - Korea", style = "heading 2")
word_export <- body_add_flextable(word_export, lsm_panss_no_neg)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_no_neg)
word_export <- body_add_break(word_export)
word_export <- body_add_par(word_export, "PANSS negative score - trained - Entire group", style = "heading 2")
word_export <- body_add_flextable(word_export, lsm_panss_w_neg)
word_export <- body_add_break(word_export)
word_export <- body_add_flextable(word_export, lsm_table_w_neg)
print(word_export, 'A002-A4_tr_1103.docx')