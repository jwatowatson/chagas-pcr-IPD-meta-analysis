make_matrix_pcr = function(pcr_dat, N_max_PCR=NA){
  
  require(RColorBrewer)
  if(is.na(N_max_PCR)) N_max_PCR=max(pcr_dat$PCR_number)
  
  K=length(unique(pcr_dat$VISIT)) * N_max_PCR
  N=length(unique(pcr_dat$ID))
  
  pcr_dat$PCR_name = apply(pcr_dat[, c('PCR_number','VISIT')], 1, function(x) paste(x[1],x[2], sep = '_'))
  pcr_mat = array(NA, dim=c(N,K))
  colnames(pcr_mat) = unique(pcr_dat$PCR_name)
  rownames(pcr_mat) = unique(pcr_dat$ID)
  
  for(id in unique(pcr_dat$ID)){
    ind = which(pcr_dat$ID==id)
    pcr_mat[id, pcr_dat$PCR_name[ind]] = pcr_dat$CT[ind]
  }
  return(pcr_mat)
  
}
plot_pcr_matrix = function(pcr_mat, my_breaks, my_break_legend, my_cols, 
                           visit_labels, arm_labels,
                           ind_y){  
  image(t(pcr_mat), x=1:ncol(pcr_mat), y = 1:nrow(pcr_mat),
        breaks = my_breaks,
        col = my_cols,
        xlab='', xaxt = 'n', yaxt='n', ylab='')
  legend('topright', inset= c(-0.1,0), title = 'CT',
         fill= my_cols, xpd=T,
         legend = my_break_legend)
  axis(1, at = (1:(ncol(pcr_mat)/9))*9 - 4.5, 
       labels = visit_labels, las=2)
  
  abline(h=ind_y[-1], lty=2)
  
  abline(v=(1:(ncol(pcr_mat)/9))*9 + 0.5, lty=2)

  axis(2, at = ind_y+15, labels = arm_labels,las=1,tick = F)
  
}