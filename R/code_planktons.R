library(rEDM)

data_plankton<-read.csv("~/Dropbox/plankton_marina/data/GOM_NA_only.csv", na = "NA")
data_plankton_original<-read.csv("~/Dropbox/plankton_marina/data/GOM_original_sst.csv")

data_plankton_ts<-ts(data_plankton[,4:19])
data_plankton_original_ts<-ts(data_plankton_original[8:22])

data_plankton_df<-as.data.frame(data_plankton)
data_plankton_original_df<-as.data.frame(data_plankton_original)

plot(data_plankton_ts[,c(2:11)])
plot(data_plankton_ts[,c(12:16)])

plot(data_plankton_original_ts[,c(1:10)])
plot(data_plankton_original_ts[,c(11:15)])


Emb1<-c()
sp<-colnames(data_plankton_ts)
par(mfrow = c(1,1),mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
for(i in 2:16){ ##para rodar o simplex em todas espécies##
  nome<-sp[i]
  lib<-c(1, NROW(data_plankton_ts[,i]))
  pred<-lib
  simplex_output<-simplex(as.data.frame(data_plankton_ts)[,i], lib = lib, pred = pred, E = 1:20, silent = TRUE)
  #par(mfrow = c(1,1),mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
  savename<-paste("~/Dropbox/plankton_marina/data/embeddings/embedding_E20",sep = "", nome, ".png")
  png(savename, width = 800, height = 600)
  plot(simplex_output$E, simplex_output$rho, type = "b", xlab = "Dimensão de Embedding (E)", 
       ylab = expression(paste(sep = "",rho)), main = nome, cex.lab=1.5)
  emax<-which.max(simplex_output$rho)
  lines(rho[emax] ~ E[emax], data = simplex_output, type = "b", col  =  "red", lwd=2, pch=19)
  Emb1[i]<-emax
  dev.off()
  
}


Emb2<-c()
sp<-colnames(data_plankton_original_ts)
par(mfrow = c(1,1),mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
for(i in 1:15){ ##para rodar o simplex em todas espécies##
  nome<-sp[i]
  lib<-c(1, NROW(data_plankton_original_ts[,i]))
  pred<-lib
  simplex_output<-simplex(as.data.frame(data_plankton_original_ts)[,i], lib = lib, pred = pred, E = 1:20, silent = TRUE)
  #par(mfrow = c(5,2),mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
  savename<-paste("~/Dropbox/plankton_marina/data/embeddings/embedding_original_E20",sep = "", nome, ".png")
  png(savename, width = 800, height = 600)
  plot(simplex_output$E, simplex_output$rho, type = "b", xlab = "Dimensão de Embedding (E)", 
       ylab = expression(paste(sep = "",rho)), main = nome, cex.lab=1.5)
  emax<-which.max(simplex_output$rho)
  lines(rho[emax] ~ E[emax], data = simplex_output, type = "b", col  =  "red", lwd=2, pch=19)
  Emb2[i]<-emax
  dev.off()
}

