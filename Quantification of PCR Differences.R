library(readr)
library(tidyverse)
library(ggpubr)
library(hrbrthemes)

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

pcr_chagas$pcr_band <- cut(pcr_chagas$CT, c(19, 34, 36, 38, 40, 42),
                           right = FALSE)

table(pcr_chagas$pcr_band, useNA = "ifany")

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

table(pcr_chagas$score, useNA = "ifany")

pcr_score <- pcr_chagas %>% 
  filter(ID != 'CGTNWOV_347') %>%  #person who only took 3 days
  group_by(USUBJID, ARMCD, VISITNUM) %>% 
  summarise(total_score = sum(score),
            avg_score = mean(score))

EOT <- data.frame(
  ARMCD = c(
    "PLACEBO", "150BZN4W", "150BZN4WFOS",
    "300BZN2W", "300BZN4W", "300BZN8W", "300BZN8WFOS"
  ),
  EOT_VISITDY = c(
    3, 8, 8, 6, 8, 12, 12
  )
)

ggplot(pcr_score, aes(x = VISITNUM, y = total_score)) +
  geom_bar(stat = "identity", col = "#00b386", fill = "#00b386") +
  facet_wrap(~ARMCD) +
  theme_bw()+
  geom_vline(data = EOT, aes(xintercept = EOT_VISITDY + 0.5), 
             lwd = 1.12, col = "#b3002d") +
  geom_text(data = EOT, mapping = aes(x = EOT_VISITDY + 2.4, y = 460, 
                                      label = "EOT"), col = "#b3002d")


pt_dosing <- data.frame(
  VISITNUM = 1:17,
  x300BZN4W_BNZ = c(0, 0, 300, 600, 900, 2400, 4500, 6600, rep(8400, 9)),
  x150BZN4W_BNZ = c(0, 0, 300, 600, 900, 2400, 4500, 6600, rep(8400, 9))/2,
  x300BZN2W_BNZ = c(0, 0, 300, 600, 900, 2400, 4200, rep(4200, 10))
)

ggplot(pcr_score %>% 
         filter(ARMCD == "300BZN4W"), aes(x = VISITNUM, y = total_score)) +
  geom_bar(stat = "identity", col = "#00b386", fill = "#00b386") +
  geom_vline(xintercept = 9, col = 2, lwd = 1.5) +
  geom_line(pt_dosing, mapping = aes(x = VISITNUM, y = x300BZN4W_BNZ / 15), col = "blue", lwd = 1.2) + 
  scale_y_continuous(
    name = "Total Score",
    sec.axis = sec_axis(~.*15, name = "Cumulative Dose"),
    breaks = seq(100,600, 100)
  )  +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(colour = "black", size = 13),
    axis.title.y = element_text(color = "#00b386", size=13),
    axis.title.y.right = element_text(color = "blue", size=13)
  ) +
  annotate("text", x = 9.8, y = 620, label = "EOT", col = 2)+
  labs(
    title = "300mg BNZ Daily for 4 Weeks"
  )

ggplot(pcr_score %>% 
         filter(ARMCD == "150BZN4W"), aes(x = VISITNUM, y = total_score)) +
  geom_bar(stat = "identity", col = "#00b386", fill = "#00b386") +
  geom_vline(xintercept = 9, col = 2, lwd = 1.5) +
  geom_line(pt_dosing, mapping = aes(x = VISITNUM, y = x150BZN4W_BNZ / 15), col = "blue", lwd = 1.2) + 
  scale_y_continuous(
    name = "Total Score",
    sec.axis = sec_axis(~.*15, name = "Cumulative Dose"),
    breaks = seq(100,600, 100)
  )  +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(colour = "black", size = 13),
    axis.title.y = element_text(color = "#00b386", size = 13),
    axis.title.y.right = element_text(color = "blue", size = 13)
  ) +
  annotate("text", x = 9.8, y = 620, label = "EOT", col = 2) +
  labs(
    title = "150mg BNZ Daily for 4 Weeks"
  )

ggplot(pcr_score %>% 
         filter(ARMCD == "300BZN2W"), aes(x = VISITNUM, y = total_score)) +
  geom_bar(stat = "identity", col = "#00b386", fill = "#00b386") +
  geom_vline(xintercept = 7, col = 2, lwd = 1.5) +
  geom_line(pt_dosing, mapping = aes(x = VISITNUM, y = x300BZN2W_BNZ / 15), col = "blue", lwd = 1.2) + 
  scale_y_continuous(
    name = "Total Score",
    sec.axis = sec_axis(~.*15, name = "Cumulative Dose"),
    breaks = seq(100,600, 100)
  )  +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(colour = "black", size = 13),
    axis.title.y = element_text(color = "#00b386", size = 13),
    axis.title.y.right = element_text(color = "blue", size = 13)
  ) +
  annotate("text", x = 7.8, y = 620, label = "EOT", col = 2) +
  labs(
    title = "300mg BNZ Daily for 2 Weeks"
  )

ggplot(pcr_score %>% 
         filter(ARMCD == "300BZN4W"), aes(x = VISITNUM, y = total_score, group = USUBJID, col = factor(USUBJID))) +
  geom_line(stat = "identity") +
  geom_vline(xintercept = 9, col = 2, lwd = 1.5) +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(colour = "black", size = 13),
    axis.title.y = element_text(color = "#00b386", size = 13),
    axis.title.y.right = element_text(color = "blue", size = 13)
  ) +
  annotate("text", x = 9.8, y = 32, label = "EOT", col = 2)
  

ggplot(pcr_score %>% 
         filter(ARMCD == "300BZN2W"), aes(x = VISITNUM, y = total_score, group = USUBJID, col = factor(USUBJID))) +
  geom_line(stat = "identity")  +
  geom_vline(xintercept = 7, col = 2, lwd = 1.5) +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(colour = "black", size = 13),
    axis.title.y = element_text(color = "#00b386", size = 13),
    axis.title.y.right = element_text(color = "blue", size = 13)
  ) +
  annotate("text", x = 7.8, y = 32, label = "EOT", col = 2)


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
         DIFF = W2_AFTER_TRT - BASELINE_SCORE,
         DIFF_M6 = VISITNUM_16_SCORE - BASELINE_SCORE,
         DIFF_M12 = VISITNUM_17_SCORE - BASELINE_SCORE) 

round(table(is.na(pcr_visit_score$BASELINE_SCORE), pcr_visit_score$ARMCD)/30, 2)

pcr_visit_score %>% 
  group_by(ARMCD) %>% 
  summarise(MEAN_DIFF = mean(na.omit(DIFF)),
            N = length(na.omit(DIFF)),
            SD = sd(na.omit(DIFF)))

ggboxplot(pcr_visit_score, 
          x = "ARMCD", y = "DIFF", 
          color = "ARMCD",
          ylab = "Score Difference", xlab = "ARMs",
          title = "Score Difference between Baseline and Two Weeks Post Treatment") +
  rotate_x_text(angle = 45)

t.test(
  (pcr_visit_score %>% filter(ARMCD == "300BZN2W", !is.na(DIFF)))$DIFF, 
  (pcr_visit_score %>% filter(ARMCD == "300BZN8W", !is.na(DIFF)))$DIFF,
  var.equal = FALSE)

# wilcox.test(DIFF ~ ARMCD, 
#             data = pcr_visit_score %>% 
#               filter((ARMCD == "300BZN8W" | ARMCD == "300BZN2W"),
#                      !is.na(DIFF)) %>% 
#               mutate(DIFF = DIFF + runif(1, -0.01, 0.01))
#             )


ggplot(pcr_chagas, aes(x = score)) +
  geom_density() +
  facet_wrap(~ARMCD) +
  theme_bw()


kruskal.test(DIFF ~ ARMCD, 
                         data = pcr_visit_score %>%
                           filter(!is.na(DIFF),
                                  ARMCD != "PLACEBO"))

kruskal.test(DIFF ~ ARMCD, 
             data = pcr_visit_score %>%
               filter(!is.na(DIFF_M6),
                      ARMCD != "PLACEBO"))

kruskal.test(DIFF ~ ARMCD, 
             data = pcr_visit_score %>%
               filter(!is.na(DIFF_M12),
                      ARMCD != "PLACEBO"))
