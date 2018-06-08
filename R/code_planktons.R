library(rEDM)
library(corrplot)

setwd("./time-series-forams/")

data_plankton<-read.csv("data/GOM/GOM_edm_use.csv", header = T, na = "NA") #loading the data frames, placing any NA as a NA chr#
data_plankton_ts<-ts(data_plankton[,4:19]) #only taking the interests series, from planktons and setting them as timeseries#
data_plankton_df<-as.data.frame(data_plankton) #turning into data frame, better format to the algorithms#

# plot(data_plankton_ts[,c(2:11)], main = "plot of plankton series") #plottings#
# plot(data_plankton_ts[,c(12:16)], main = "plot of plankton series") #same as above#

# used<-c("everything", "all.obs", "complete.obs", "na.or.complete","pairwise.complete.obs")
# cor_plankton_ts<-cor(data_plankton_ts[,c(2:16)], use = used[4], method = "pearson")
# corrplot(cor_plankton_ts, method = "square", title = "data_plankton_ts")


###SIMPLEX MEHOD TO EMBEDDING DIMENSION DETERMINATION###

tp_simplex1<- "default"
Emb1<-c() #embedding vector#
sp<-colnames(data_plankton_ts) #picking the colnames#
par(mfrow = c(1,1),mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0)) #parameters to the plots#
for(i in 2:16){ ##to run over all species of planktons##
  nome<-sp[i] #name variable to set plot with it#
  lib<-c(1, NROW(data_plankton_ts[,i])) #defining the lib size used to simplex algorithmic#
  pred<-lib #same as above#
  simplex_output<-simplex(as.data.frame(data_plankton_ts)[,i], 
                          lib = lib, 
                          pred = pred, 
                          E = 1:20, 
                          silent = TRUE, 
                          tp = tp_simplex1
                          ) #simplex method function calling#
  #par(mfrow = c(1,1),mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
  #savename<-paste("~/Dropbox/plankton_marina/data/embeddings/embedding_",sep = "", nome, tp_simplex1, ".png") #saving name to be set on the archive and as title of the plot#
  #png(savename, width = 800, height = 600) #png ambient to save the figure as a png figure#
  #plot(simplex_output$E, simplex_output$rho, type = "b", xlab = "Dimensão de Embedding (E)", 
  #     ylab = expression(paste(sep = "",rho)), main = paste(nome, sep = "_tp_",tp_simplex1), cex.lab=1.5) #plot function with extra parameters to plotting#
  emax<-which.max(simplex_output$rho) #picking the maximum of embedding dimension#
  #lines(rho[emax] ~ E[emax], data = simplex_output, type = "b", col  =  "red", lwd=2, pch=19) #setting the maximum emebedding dimension on the plot#
  Emb1[i]<-emax #placing the maximum into the embedding vector for the record and future uses#
  #dev.off() #shutting down the png ambient#
}

data.frame(sp,Emb1)

# Marina (pelos graficos de E 1:20, no olho):
Emb1 = c(NA, 3,3,5,2,2,2,1,2,2,5,3,10,3,5,5)

###S-MAP METHOD TO EMBEDING DIMENSIONS DETERMINATION##
##Same as for the simplex method but by s-map, which considers the whole attractor points##
if (!file.exists("output/GOM/theta_plots")){ dir.create("output/GOM/theta_plots/") } # creating output folder for plots
theta1<-c()
sp<-colnames(data_plankton_ts)
par(mfrow = c(1,1), mar = c(4,4,1,1), mgp =  c(2.5,1,0))
for (i in 2:16) {
  nome<-sp[i]
  lib<-c(1, NROW(data_plankton_ts[,i]))
  pred<-lib
  s_map_output<-s_map(as.data.frame(data_plankton_ts)[,i], 
                   lib = lib, 
                   pred = pred, 
                   E = Emb1[i], 
                   theta = c(0, 1e-04, 3e-04, 0.001,0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8), 
                   silent = TRUE, 
                   save_smap_coefficients = TRUE)
  savename<-paste("output/GOM/theta_plots/theta_",sep = "", nome, ".png") #saving name to be set on the archive and as title of the plot#
  png(savename, width = 800, height = 600) #png ambient to save the figure as a png figure#
  plot(s_map_output$theta, s_map_output$rho, type = "b", xlab = expression(paste("Theta (", sep = "", theta, ")")),ylab = expression(paste(sep = "",rho)), main = nome , cex.lab = 1.5)
  thetamax<-which.max(s_map_output$rho)
  lines(rho[thetamax]~theta[thetamax], data = s_map_output, type = "b", col = "red", lwd = 2, pch = 19)
  theta1[i]<-thetamax
  dev.off()
  
}
data.frame(sp,theta1)



###Embedding by CCM, which optimatize the embedding dimensions by CCM function###
if (!file.exists("output/GOM/ccm_embed_plots")){ dir.create("output/GOM/ccm_embed_plots/") } # creating output folder for plots

E_n1<-as.data.frame(matrix(data = NA, nrow = 14, ncol = 15))
sp<-colnames(data_plankton_ts)
colnames(E_n1)<-sp[2:16]

for (i in 2:15) { ##To run over each species##
  for (j in 2:16) { ##To run over each species differing from the first one##
    if(i == j)next
    E_max_ccm<-data.frame()
    #E_max_ccm<-colnames(sp[2:16])
    E.test.1=NULL
    cmxy.t<-list()
    for(E.t in 1:20){ ##to run over different embedding dimensions##
      cmxy.t <- ccm(data_plankton_ts, E = E.t, lib_column = sp[i], target_column = sp[j],
                  lib_sizes = 203, num_samples = 1, tp=-1,random_libs = F, silent = TRUE)
      E.test.1=rbind(E.test.1,cmxy.t)
  }
    E_n1[i-1, j-1] <- which.max(E.test.1$rho) # the optimal E
    
    savename<-paste("output/GOM/ccm_embed_plots/ccm_emb",sp[i], "xmap", sp[j],".png",sep = "_") #saving name to be set on the archive and as title of the plot#
    png(savename, width = 800, height = 600) #png ambient to save the figure as a png figure#
    par(mfrow = c(1,1), mar = c(4,4,1,1), mgp =  c(2.5,1,0))
    plot(E.test.1$E, E.test.1$rho, type = "b", 
         ylab = expression(paste(sep = "", rho)), 
         xlab = "Embedding dimension (E)",
         main = paste(sp[i], sep = "_", "xmap", sp[j]))
    lines(rho[E_n1[i-1, j-1]]~E[E_n1[i-1, j-1]], data = E.test.1, type = "b", col =  "red", lwd = 2, pch = 19)
    dev.off()
  }
}



# CCM analysis with varying library size (L)
    #libs<-c(1, NROW(data_plankton_ts[,i]))

##CCM from each species to SST###
par(mfrow = c(1,1), mar = c(4,4,1,1), mgp =  c(2.5,1,0))
for (i in 2:15) { ##To run over each species##
    spi_xmap_spj <- ccm(data_plankton_ts, E=E_n1[i-1, j-1],lib_column=sp[i], target_column="sst",
                        lib_sizes=seq(1,203, by = 10), replace=T, silent = TRUE)
    
    spi_xmap_spj_means<-ccm_means(spi_xmap_spj)
    
    plot(spi_xmap_spj_means$lib_size, pmax(0, spi_xmap_spj_means$rho), 
         data = spi_xmap_spj_means, 
         type = "l", 
         xlab = "Library size (L)", 
         ylab = expression(paste("Mapping skil (", sep = "", rho, ")")),
         col = "red", 
         main = paste(sp[i], sep = "_", "xmap_sst"),
         lwd = 2, 
         pch = 19, 
         cex.lab = 1.5)
    #abline()
}


###CCM between species###
par(mfrow = c(1,1), mar = c(4,4,1,1), mgp =  c(2.5,1,0))
for (i in 2:15) { ##To run over each species##
  for (j in 2:15) {
    if(i == j)next
    spi_xmap_spj <- ccm(data_plankton_ts, E=E_n1[i-1, j-1],lib_column=sp[i], target_column=sp[j],
                  lib_sizes=seq(1,203, by = 10), replace=T, silent = TRUE)
    
    spi_xmap_spj_means<-ccm_means(spi_xmap_spj)
    
    plot(spi_xmap_spj_means$lib_size, pmax(0, spi_xmap_spj_means$rho), 
         data = spi_xmap_spj_means, 
         type = "l", 
         xlab = "Library size (L)", 
         ylab = expression(paste("Mapping skil (", sep = "", rho, ")")),
         col = "red", 
         main = paste(sp[i], sep = "_", "xmap", sp[j]),
         lwd = 2, 
         pch = 19, 
         cex.lab = 1.5)
    #abline()

  }
}


