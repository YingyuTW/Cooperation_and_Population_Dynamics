wd = "D:/my/working/directory/"
setwd(wd)
ts_data_processing("ts_social")
ts_data_visualization("ts_social")
myfft("ts_social")

#-----------------------------------------------------------------------------

transNA = function( df ){
  for ( i in 1:ncol(df) ){
    df[df[,i]==-1,i] = NA
  }
  return(df)
}

ts_data_processing = function(fname){
  df = read.table(paste0(fname, ".txt"), sep = "\t")
  colnames(df) = c("Resource", "Benefit", "CostRate", "Round", "Time", "FlucAmp", "FlucPeriod", "pSize", "cProp", "cDeg", "IC" )
  df = transNA(df)
  dyn.type  = c("pSize", "cProp", "cDeg")
  round.num = length(which(df$Time == 0))
  time.max  = max(df$Time)
  span      = df$Time[2] - df$Time[1]
  ncut      = time.max/span + 1
  df_new    = as.data.frame(matrix(NA, ncut, 10))
  colnames(df_new) = c("Time", "Resource", "pSizeM", "pSizeU", "pSizeL", "cPropM", "cPropU", "cPropL", "cDegM", "cDegU", "cDegL")
  
  i = 1
  t = 0
  while (t <= time.max) {
    id = which(df$Time == t)
    df_new$Time[i] = t
    df_new$Resource[i] = df$Resource[id[1]]
    df_new[i, c("pSizeM", "cPropM", "cDegM")] = colMeans(df[id, c("pSize", "cProp", "cDeg")])
    tmpSD = apply(df[id, c("pSize", "cProp", "cDeg")], 2, sd)
    error = qnorm(0.975)*tmpSD/sqrt(round.num)
    df_new[i, c("pSizeU", "cPropU", "cDegU")] = df_new[i, c("pSizeM", "cPropM", "cDegM")] + error
    df_new[i, c("pSizeL", "cPropL", "cDegL")] = df_new[i, c("pSizeM", "cPropM", "cDegM")] - error
    i = i + 1
    t = t + span
  }
  
  write.csv(df_new, paste0(fname, ".csv"), row.names = FALSE)
}


ts_data_visualization = function(fname){
  
  df = read.csv(paste0(fname, ".csv"))
  
  for (type in 1:3){
    gname = paste0(fname, "_", dyn.type[type], ".pdf")
    pdf(file = gname, width = 4.3, height = 3.7)
    if (type == 1){
      par(mar = c(4.5, 4.5, 2, 4.5))
      plot(df$Time, df$pSizeM, xlab = "Time", ylab = "Population size", ylim = c(0,2000), type = "l")
      legend("topright", legend=c("Population size", "Environmental resource availability"), 
             lwd = 1, col=c("black", "#418CCB"), bty = 'n', cex = 0.75)
    }else if (type == 2){
      par(mar = c(4.5, 4.5, 2, 4.5))
      plot(df$Time, df$cPropM, xlab = "Time", ylab = "Proportion of cooperators", ylim = c(0,1.2), type = "l")
      legend("topright", legend=c("Proportion of cooperators", "Environmental resource availability"), 
             lwd = 1, col=c("black", "#418CCB"), bty = 'n', cex = 0.75)
    }else if (type == 3){
      par(mar = c(4.5, 4.5, 2, 4.5))
      plot(df$Time, df$cDegM, xlab = "Time", ylab = "Average degree of cooperation", ylim = c(0,1), type = "l")
      legend("topright", legend=c("Average degree of cooperation", "Environmental resource availability"), 
             lwd = 1, col=c("black", "#418CCB"), bty = 'n', cex = 0.75)
    }
    par(new = T)
    plot(df$Time, df$Resource, type = 'l', col = "#418CCB", axes = FALSE, xlab = NA, ylab = NA, ylim = c(0, 25))
    axis(side = 4)
    mtext(side = 4, line = 3, 'Environmental resource availability')
    dev.off()
  }
}

myfft = function(fname){
  df = read.table(paste0(fname, ".txt"), sep = "\t")
  colnames(df) = c("Resource", "Benefit", "CostRate", "Round", "Time", "FlucAmp", "FlucPeriod", "pSize", "cProp", "cDeg", "IC" )
  df = transNA(df)
  round.num = length(which(df$Time == 0))
  tmp = as.data.frame(matrix(NA, 900, round.num))
  for (i in 1:round.num){
    id = which(df$Round == i-1 & df$Time >= 1000)
    ds = (df$pSize[id]-mean(df$pSize[id]))/mean(df$pSize[id])
    X.k <- fft(ds) # find all harmonics with fft()
    tmp[,i] = X.k[1:900]
  }
  Frequency <- (1:900)/9000
  X.k.m = rowMeans(tmp)
  magn.l = Mod(X.k.m)/sqrt(9000)
  pdf(file = paste0(fname, "_fft.pdf"), width = 4, height = 3.7)
  par(mar = c(4.5, 4.5, 2, 2))
  plot(Frequency[1:90], magn.l[1:90], xlab = "Frequency", ylab = "Magnitude", type = 'l')
  dev.off()
}
