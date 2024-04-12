### Replication of BENDITA Results
library(readr)
library(tidyverse)
library(boot)

source("chagas_adam.R")
source("Functions.R")

dm_chagas <- read_csv("~/Chagas/Data 2024.03.08/DM 2024-03-08.csv", guess_max = Inf)
in_chagas <- read_csv("~/Chagas/Data 2024.03.08/IN 2024-03-08.csv", guess_max = Inf)
mb_chagas <- read_csv("~/Chagas/Data 2024.03.08/MB 2024-03-08.csv", guess_max = Inf)
pc_chagas <- read_csv("~/Chagas/Data 2024.03.08/PC 2024-03-08.csv", guess_max = Inf)
ti_chagas <- read_csv("~/Chagas/Data 2024.03.08/TI.csv", guess_max = Inf)
ts_chagas <- read_csv("~/Chagas/Data 2024.03.08/TS.csv", guess_max = Inf)
tv_chagas <- read_csv("~/Chagas/Data 2024.03.08/TV.csv", guess_max = Inf)

pcr_chagas <- chagas_adam(
  dm_chagas,
  in_chagas,
  mb_chagas,
  ts_chagas
)

{
pcr_chagas$pcr_band <- cut(pcr_chagas$CT, c(19, 34, 36, 38, 40, 42),
                           right = FALSE)

pcr_chagas$score = NA
for (i in 1:nrow(pcr_chagas)) {
  if(pcr_chagas$pcr_band[i] == "[19,34)"){
    pcr_chagas$score[i] = 4
  } else if (pcr_chagas$pcr_band[i] == "[34,36)"){
    pcr_chagas$score[i] = 3
  } else if (pcr_chagas$pcr_band[i] == "[36,38)"){
    pcr_chagas$score[i] = 2
  } else if (pcr_chagas$pcr_band[i] == "[38,40)"){
    pcr_chagas$score[i] = 1
  } else {
    pcr_chagas$score[i] = 0
  }
}

pcr_score <- pcr_chagas %>% 
#  filter(ID != 'CGTNWOV_347') %>%  #person who only took 3 days
  group_by(USUBJID, ARMCD, VISITNUM) %>% 
  summarise(total_score = sum(score),
            avg_score = mean(score))

pcr_visit_score <- pcr_score %>% 
  arrange(VISITNUM) %>% 
  pivot_wider(id_cols = c(USUBJID, ARMCD),
              names_from = VISITNUM,
              values_from = total_score,
              names_glue = "VISITNUM_{.name}_SCORE") %>% 
  mutate(BASELINE_SCORE = VISITNUM_3_SCORE,
         W2_AFTER_TRT = if_else(ARMCD == "PLACEBO", VISITNUM_6_SCORE,
                                if_else(ARMCD == "150BZN4W" | ARMCD == "150BZN4WFOS" | ARMCD == "300BZN4W", 
                                        VISITNUM_10_SCORE,
                                        if_else(ARMCD == "300BZN2W", VISITNUM_8_SCORE, VISITNUM_13_SCORE))),
         DIFF = W2_AFTER_TRT - BASELINE_SCORE) 
}
# 202/210 completed 
sapply(pcr_visit_score, function(x) table(is.na(x), useNA = "ifany"))

table1::table1(~ AGE | ARMCD, data = dm_chagas)

pcr_visit_score$cum_score_EOT_M6 <- NA
pcr_visit_score$cum_score_EOT_M12 <- NA
for (i in 1:nrow(pcr_visit_score)) {
  col_6M = "VISITNUM_16_SCORE"
  col_12M = "VISITNUM_17_SCORE"
  if (pcr_visit_score$ARMCD[i] == "PLACEBO"){
    col_EOT = "VISITNUM_3_SCORE"
    pcr_visit_score$cum_score_EOT_M6[i] <- 
      sum(as.numeric(pcr_visit_score[i, which(names(pcr_visit_score) == col_EOT):which(names(pcr_visit_score) == col_6M)]), 
          na.rm = TRUE)
    
    pcr_visit_score$cum_score_EOT_M12[i] <- 
      sum(as.numeric(pcr_visit_score[i, which(names(pcr_visit_score) == col_EOT):which(names(pcr_visit_score) == col_12M)]), 
          na.rm = TRUE)
  } else if (pcr_visit_score$ARMCD[i] == "150BZN4W" | pcr_visit_score$ARMCD[i] == "150BZN4WFOS" | pcr_visit_score$ARMCD[i] == "300BZN4W") {
    col_EOT = "VISITNUM_10_SCORE"
    pcr_visit_score$cum_score_EOT_M6[i] <- 
      sum(as.numeric(pcr_visit_score[i, which(names(pcr_visit_score) == col_EOT):which(names(pcr_visit_score) == col_6M)]), 
          na.rm = TRUE)
    
    pcr_visit_score$cum_score_EOT_M12[i] <- 
      sum(as.numeric(pcr_visit_score[i, which(names(pcr_visit_score) == col_EOT):which(names(pcr_visit_score) == col_12M)]), 
          na.rm = TRUE)
  } else if (pcr_visit_score$ARMCD[i] == "300BZN2W") {
    col_EOT = "VISITNUM_7_SCORE"
    pcr_visit_score$cum_score_EOT_M6[i] <- 
      sum(as.numeric(pcr_visit_score[i, which(names(pcr_visit_score) == col_EOT):which(names(pcr_visit_score) == col_6M)]), 
          na.rm = TRUE)
    
    pcr_visit_score$cum_score_EOT_M12[i] <- 
      sum(as.numeric(pcr_visit_score[i, which(names(pcr_visit_score) == col_EOT):which(names(pcr_visit_score) == col_12M)]), 
          na.rm = TRUE)
  } else {
    col_EOT = "VISITNUM_13_SCORE"
    pcr_visit_score$cum_score_EOT_M6[i] <- 
      sum(as.numeric(pcr_visit_score[i, which(names(pcr_visit_score) == col_EOT):which(names(pcr_visit_score) == col_6M)]), 
          na.rm = TRUE)
    
    pcr_visit_score$cum_score_EOT_M12[i] <- 
      sum(as.numeric(pcr_visit_score[i, which(names(pcr_visit_score) == col_EOT):which(names(pcr_visit_score) == col_12M)]), 
          na.rm = TRUE)
  }
}

pcr_visit_score %>% 
  select(USUBJID, ARMCD, cum_score_EOT_M6, everything()) %>% 
  filter(ARMCD == "300BZN8W", 
         !is.na(VISITNUM_16_SCORE)) %>% 
  View

pcr_visit_score %>% 
  group_by(ARMCD) %>% 
  summarise(n_NA_M6 = length(which(is.na(VISITNUM_16_SCORE))),
            n_SPC_M6 = length(which(cum_score_EOT_M6 == 0 & !is.na(VISITNUM_16_SCORE))),
            n_M6 = 30 - n_NA_M6,
            n_NA_M12 = length(which(is.na(VISITNUM_17_SCORE))),
            n_SPC_M12 = length(which(cum_score_EOT_M12 == 0 & !is.na(VISITNUM_17_SCORE))),
            n_M12 = 30 - n_NA_M12) 
