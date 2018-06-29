# code_planktons.R do Rafa

embed_dim <- function(data_ts, emb_dim, trap_name, overwrite){    
  
 if(overwrite == TRUE | !file.exists(paste("output/",trap_name,"/embedding_plots",sep=""))){
    
  if (!file.exists(paste("output/",trap_name,"/embedding_plots",sep=""))){ dir.create(paste("output/",trap_name,"/embedding_plots/",sep="")) }

  Emb1<-data.frame()
  sp<-colnames(data_ts)
  # Optimized embeding dimension by eye (looking at each species plots)
  emax_eye = c(3,3,5,2,2,2,1,2,2,5,3,10,3,5,5)

  for(i in 2:length(colnames(data_ts))){ # simplex for each variable

    lib <- c(1, NROW(data_ts[,i]))
    pred <- lib
    simplex_output <- simplex(as.data.frame(data_ts)[,i], lib = lib, pred = pred, E = 1:emb_dim, silent = T)

    Emb1[i,"species"]<-sp[i]
    Emb1[i,"emax_auto"]<-which.max(simplex_output$rho)
    Emb1[i,"rho_auto"]<-simplex_output[which.max(simplex_output$rho),"rho"]
   # Emb1[i,"emax_eye"]<-emax_eye[i-1]
   # Emb1[i,"rho_eye"]<-simplex_output[emax_eye[i-1],"rho"]
    
    # Creating a objetive criteria for embedding dimension:
    # If there are other rhos less than 0.02 from the max rho, take min Embedding dimension
    #"Problematic cases": sst, G. crassaformis, G. falconensis
    # Emb1[i,"emax_opt"]<-min(simplex_output[which(simplex_output$rho > max(simplex_output$rho) - 0.02),"E"]) # finding other maxima
    # Emb1[i,"rho_opt"]<-simplex_output[Emb1[i,"emax_opt"],"rho"]

    savename<-paste("output/",trap_name,"/embedding_plots/",sp[i], "_E",emb_dim,".png",sep="")
    png(savename, width = 800, height = 600)
      par(mfrow = c(1,1),mar = c(4, 5, 1, 1), mgp = c(2.5, 1, 0))
      plot(x=simplex_output$E, y=simplex_output$rho, type = "b", xlab = "Embedding dimension (E)", 
           ylab = expression(paste("Forecast skill (",rho,")",sep = "")), cex.lab=1.5)
      points(max(simplex_output$rho) ~ simplex_output$E[which.max(simplex_output$rho)],  type = "b", col  =  "red", lwd=2, pch=19)
    dev.off()

  } # for

  Emb1 <- Emb1[-1,]
  
  # Emb1[,"rho_diff"] <- Emb1[,"rho_auto"] - Emb1[,"rho_opt"]
  
  write.csv(Emb1, paste("output/",trap_name,"/embedding_plots_",trap_name,".csv",sep=""), row.names = F)
  return(Emb1)
  
  
}else{
  
  Emb1 <- read.csv(paste("output/",trap_name,"/embedding_plots_",trap_name,".csv",sep=""), header = TRUE, stringsAsFactors = FALSE)
  return (Emb1)
 }
  
  
} # function

