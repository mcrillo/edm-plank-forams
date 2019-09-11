# function_cormtest

cormatrix <- function(mat, cormethod, conf.level = 0.95){
  mat <- as.matrix(mat)
  n <- ncol(mat)
  cor.mat <- p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(cor.mat) <- 1
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level, 
                      method=cormethod, alternative="two.sided")
      cor.mat[i,j] <- cor.mat[j,i] <- tmp$estimate
      p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
      
      if(cormethod=="pearson"){
        lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$conf.int[1]
        uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$conf.int[2]
      }
    }
  }
  return(list(corr = cor.mat, p_value = p.mat, lower_ci = lowCI.mat, upper_ci = uppCI.mat))
}
