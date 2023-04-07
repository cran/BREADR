utils::globalVariables(c("deam", "pop", "chr",
                         "keep", "site", "nsnps",
                         "pmr", "mismatch", "ave_rel",
                         "Same_Twins", "First_Degree", "Second_Degree",
                         "Unrelated", "normConst", "relationship", "pair",
                         "n", "sd", "everything", "col_names_tibble",
                         "meanPMR", "ymin", "ymax", "model", "p.sim", "y",
                         "posterior", "label"))
getfilter <- function(x,gap){
  i <- 1
  res <- c(1)
  while(T){
    newMin <- which((x-x[i])>=gap)
    if(length(newMin)==0){
      break
    }else{
      res <- c(res,newMin%>%min())
      i <- newMin%>%min()
    }
  }
  return(1:length(x)%in%res)
}
weightedBinom <- function(x,N,p,lambda=10,uppern=1e3,log=F){
  pVec <- (1-0.5^(1:(uppern+1)))*p
  logP <- dtruncatedPoisson(0:uppern,d=2,lambda=lambda,log=T)+stats::dbinom(x,N,pVec,log=T)
  if(log){
    return(matrixStats::logSumExp(logP))
  }else{
    return(sum(exp(logP)))
  }
}
dtruncatedPoisson <- function(k,d,lambda,log=F){
  numerator <- k*log(lambda)
  denominator <- lfactorial(k)+log(exp(lambda)-sum(lambda^(0:d)/(factorial(0:d))))
  if(log){
    return(numerator-denominator+log(k%%1==0)+log(k>d))
  }else{
    return(exp(numerator-denominator)*(k%%1==0)*(k>d))
  }
}
ggcolorhue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}
convertP <- function(p,pbar){
  suppressWarnings(ifelse(p<pbar,log(1-p/pbar)/log(0.5)-1,Inf)) %>% return()
}
label_formatter <- function(x){
  label <- dplyr::case_when(
    x >= 0.01 ~ as.character(round(x, 4)),
    x < 1e-6 ~ "{}<1 %*% 10^-6",
    TRUE ~ stringr::str_replace(
      formatC(x, digits = 2, format = "e"),
      "e", " %*% 10^"
    )
  )
  return(label)
}

