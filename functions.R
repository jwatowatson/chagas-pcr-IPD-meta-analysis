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
  
  image(t(pcr_mat), x=1:ncol(pcr_mat), y = 1:nrow(pcr_mat),
        breaks = c(10, 34, 36, 38, 39.99999,40),
        col = (c(brewer.pal(4, 'RdYlBu'),brewer.pal(9, 'Pastel1')[9])),
        xlab='', xaxt = 'n', yaxt='n', ylab='')
  legend('topright', inset= c(-0.2,0), title = 'Mean CT',
         fill= c(brewer.pal(4, 'RdYlBu'),'lightgrey'), xpd=T,
         legend = c('<34','34-36','36-38','38-40','>=40'))
  axis(1, at = (1:(ncol(pcr_mat)/9))*9 - 4.5, 
       labels = gsub(x = unique(pcr_dat$VISIT),
                     pattern = '(',replacement = '\n(',fixed = T), las=2)
  xx2 = pcr_dat %>% distinct(ID,.keep_all = T)
  ind_y = which(!duplicated(xx2$ARM))-0.5
  abline(h=ind_y[-1], lty=2)
  
  abline(v=(1:(ncol(pcr_mat)/9))*9 + 0.5, lty=2)
  
  my_labs=unique(pcr_dat$ARM)
  my_labs=gsub(pattern = 'BENZNIDAZOLE',replacement = 'BNZ',x = my_labs)
  my_labs=gsub(pattern = 'FOSRAVUCONAZOLE',replacement = '\nE1224',x = my_labs)
  my_labs=gsub(pattern = ' WEEKLY',replacement = '/WK',x = my_labs)
  
  axis(2, at = ind_y+15, labels = my_labs,las=1,tick = F)
  
}