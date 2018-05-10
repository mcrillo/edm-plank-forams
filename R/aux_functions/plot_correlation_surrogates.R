
plot_correlation_surrogates <- function(corr_surrog, trap_name, overwrite){     

  
if(overwrite == TRUE | !file.exists(paste("output/",trap_name,"/correlation_surrogates",sep=""))){
    
  if (!file.exists(paste("output/",trap_name,"/correlation_surrogates",sep=""))){ dir.create(paste("output/",trap_name,"/correlation_surrogates",sep="")) }
  
  for(i in unique(c(corr_surrog$var1,corr_surrog$var2))){ # i = focus variable
      
    corr_subset <- corr_surrog[c(which(corr_surrog$var1 == i), which(corr_surrog$var2 == i)),]
    # corr_subset[1:14,1:10]
    
    # Generating vector with variables names that are being compared with the focus variable i
    corr_subset$group <- c(corr_subset$var2[which(corr_subset$var2 != i)],corr_subset$var1[which(corr_subset$var1 != i)])
    # corr_subset[,c("var1","var2", "group")]

    
    # Plot
    subset_melt <- melt(corr_subset[, (which(names(corr_subset)=="corr_resid_p")+1):ncol(corr_subset)], id = "group")
    # (which(names(corr_subset)=="corr_resid_p")+1) is the same as which(names(corr_subset)=="V1"), but the name "V1" can be "V1.tau" depending on the correlation test used, so this 'which()+1'is to avoid this problem when plotting
    
    p <- ggplot() + 
         geom_boxplot(data = subset_melt, aes(factor(group), value)) +
         labs(y = "Correlation", x = element_blank()) +
         theme_bw() + ylim(-1, 1.1) +
         geom_point(data = corr_subset, aes(x = factor(group), y = corr_series), color = 'red',size = 2) +
         geom_point(data = corr_subset, aes(x = factor(group), y = corr_splines), color = 'blue',shape = 115,size = 2) + # shape: s
         geom_point(data = corr_subset, aes(x = factor(group), y = corr_resid), color = 'blue',shape = 114,size = 2) + # shape: r
         theme(axis.title=element_text(size=18, face="bold"))

    if(any(corr_subset$corr_series_p_surrog=="signif")){
      p <- p + geom_point(data = corr_subset[which(corr_subset$corr_series_p_surrog =="signif"),], aes(x = factor(group), y = 1.05), color = 'red',shape = 8, size = 3)
    }

    pdf(file =  paste("output/",trap_name,"/correlation_surrogates/",i,"_GOM.pdf",sep=""), width=length(corr_subset$group), height=5, paper = "special")
      print(p)
    dev.off()  
    
    rm(corr_subset)
    rm(subset_melt)
    
  }
}# if
}
  
