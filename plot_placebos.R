load('RData/pcr_chagas.RData')
pcr_chagas_ss_input = pcr_chagas %>% 
  filter(ARM=='PLACEBO') %>%
  group_by(USUBJID, VISIT, MBREFID) %>%
  mutate(
    N_PCRs_sample = length(CT),
    N_pos_sample = sum(CT<40)) %>%
  group_by(USUBJID,VISIT) %>% arrange(USUBJID, Day_frm_rand, MBREFID, MBGRPID) %>%
  mutate(PCR_number = 1:length(CT)) %>%
  group_by(USUBJID) %>%
  mutate(Baseline_CT = mean(CT[which(VISIT=='Screening')])) %>% ungroup() %>%
  arrange(Baseline_CT, USUBJID, Day_frm_rand, MBREFID, MBGRPID) %>%
  mutate(x_val = case_when(
    MBREFID=='Sample 1' ~ -0.2,
    MBREFID=='Sample 2' ~ 0,
    MBREFID=='Sample 3' ~ 0.2
  ))

pdf('placebo_plots.pdf',width = 10, height = 10)
pcr_chagas_ss_input$Timepoint =
  as.numeric(as.factor(pcr_chagas_ss_input$VISIT_numeric))
pcr_chagas_ss_input$MBREFID = as.factor(pcr_chagas_ss_input$MBREFID)
par(mfrow=c(2,2),las=1, mar=c(3,3,3,1))
for(id in unique(pcr_chagas_ss_input$USUBJID)){
  pcr_id = pcr_chagas_ss_input%>%filter(USUBJID==id)
  plot(pcr_id$Timepoint+pcr_id$x_val, pcr_id$CT, col=pcr_id$MBREFID,
       panel.first = grid(), xaxt='n', ylim = range(pcr_chagas_ss_input$CT))
  title(paste(unique(pcr_id$STUDYID), unique(pcr_id$USUBJID),sep = '-'))
}
dev.off()



pdf('placebo_plots_time.pdf',width = 10, height = 10)
pcr_chagas_ss_input$Timepoint =
  as.numeric(as.factor(pcr_chagas_ss_input$VISIT_numeric))
pcr_chagas_ss_input$MBREFID = as.factor(pcr_chagas_ss_input$MBREFID)
par(mfrow=c(2,2),las=1, mar=c(3,3,3,1))
for(id in unique(pcr_chagas_ss_input$USUBJID)){
  pcr_id = pcr_chagas_ss_input%>%filter(USUBJID==id)
  plot(pcr_id$Day_frm_rand+pcr_id$x_val, pcr_id$CT, col=pcr_id$MBREFID,
       panel.first = grid(), ylim = range(pcr_chagas_ss_input$CT))
  title(paste(unique(pcr_id$STUDYID), unique(pcr_id$USUBJID),sep = '-'))
}
dev.off()