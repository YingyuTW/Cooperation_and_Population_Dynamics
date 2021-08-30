wd = "D:/my/working/directory/"
setwd(wd)
repro_data_processing("repro_harsh")
repro_data_visualization("repro_harsh")

#-----------------------------------------------------------------------------

transNA = function( df ){
  for ( i in 1:ncol(df) ){
    df[df[,i]==-1,i] = NA
  }
  return(df)
}

repro_data_processing = function(fname){
  
  df = read.table(paste0(fname, ".txt"), sep = "\t")
  colnames(df) = c("Resource", "Benefit", "CostRate", "EndTime", "Round", "PatchID", "IndivDeg", "IndivProd", 
                   "PatchSize", "cNum", "cProp", "cDeg", "cDegTotal", "ProdTotal", "AverProd")
  df = transNA(df)
  
  s = 0
  r = numeric(0)
  for (i in 0:1999){
    if (any(df$Round == i)){
      s = s + 1
      r[s] = i
    }
  }
  
  df_new = as.data.frame(matrix(0, s*90, 15))
  colnames(df_new) = c("Resource", "Benefit", "CostRate", "EndTime", "Round", "PatchID", "IndivDeg", "IndivProd", 
                       "PatchSize", "cNum", "cProp", "cDeg", "cDegTotal", "ProdTotal", "AverProd")
  for (j in 1:s){
    for (p in 0:89){
      id = which(df$Round == r[j] & df$PatchID == p)
      df_new[(j-1)*90+p+1,] = df[id[1],]
    }
  }
  df_new = df_new[which(!is.na(df_new$Resource)),]
  write.csv(df_new, paste0(fname, ".csv"), row.names = FALSE)
}


repro_data_visualization = function(fname){
  
  df = read.csv(paste0(fname, ".csv"))
  df$CooBenefits = df$Benefit*df$cDegTotal
  cortest = cor.test(df$ProdTotal, df$CooBenefits)
  
  pdf(file = paste0(fname, ".pdf"), width = 4.0, height = 3.7)
  par(mar = c(4.5, 4.5, 2, 2))
  plot(df$CooBenefits, df$ProdTotal, pch = 16, col = "#E9736D", cex = 0.6, ylim=c(0, 35), xlab = "Cooperation benefit", ylab = "Total reproductive output")
  if (cortest$p.value < 0.01){
    text(max(df$CooBenefits)/2, 30, paste0("Correlation coefficient = ", round(cortest$estimate, 3)))
  }else{
    text(max(df$CooBenefits)/2, 30, paste0("Correlation coefficient = ",0))
  }
  dev.off()
}
