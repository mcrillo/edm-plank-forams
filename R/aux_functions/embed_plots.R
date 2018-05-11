# code_planktons.R do Rafa

embed_plots <- function(data_ts, emb_dim, trap_name, overwrite){    
  
 if(overwrite == TRUE | !file.exists(paste("output/",trap_name,"/embedding_plots",sep=""))){
    
  if (!file.exists(paste("output/",trap_name,"/embedding_plots",sep=""))){ dir.create(paste("output/",trap_name,"/embedding_plots/",sep="")) }

  Emb1<-data.frame()
  sp<-colnames(data_ts)
  
  for(i in 2:length(colnames(data_ts))){ # simplex for each variable
    nome<-sp[i]
    lib<-c(1, NROW(data_ts[,i]))
    pred<-lib
    simplex_output<-simplex(as.data.frame(data_ts)[,i], lib = lib, pred = pred, E = 1:emb_dim, silent = T)
    emax<-which.max(simplex_output$rho)
    
    Emb1[i,"variable"]<-nome
    Emb1[i,"emax"]<-emax
    
    savename<-paste("output/",trap_name,"/embedding_plots/",nome, "_E",emb_dim,".png",sep="")
    png(savename, width = 800, height = 600)
      par(mfrow = c(1,1),mar = c(4, 5, 1, 1), mgp = c(2.5, 1, 0))
      plot(simplex_output$E, simplex_output$rho, type = "b", xlab = "Embedding dimension (E)", 
           ylab = expression(paste("Forecast skill (",rho,")",sep = "")), cex.lab=1.5)
      points(rho[emax] ~ E[emax], data = simplex_output, type = "b", col  =  "red", lwd=2, pch=19)
    dev.off()
  
  } # for

  Emb1 <- Emb1[-1,]
  
  write.csv(Emb1, paste("output/",trap_name,"/embedding_plots_",trap_name,".csv",sep=""), row.names = F)
  return(Emb1)
  
  
}else{
  
  Emb1 <- read.csv(paste("output/",trap_name,"/embedding_plots_",trap_name,".csv",sep=""), header = TRUE, stringsAsFactors = FALSE)
  return (Emb1)
 }
  
  
} # function

