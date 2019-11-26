

MVApowerstat.null <- function(stats.null){
  stats.null  <-suppressMessages(suppressWarnings(anova(suppressMessages(suppressWarnings(extend.manyglm(fit,N=N,coeffs = coeff.null,use.design=use.design,newdata=newdata))), nBoot=1, test=test,show.time = "none")))$table[pred,3]
}

MVApowerstat.alt <- function(stats){
  stats     <-suppressMessages(suppressWarnings(anova(suppressMessages(suppressWarnings(extend.manyglm(fit,N=N,coeffs = coeffs,use.design=use.design,newdata=newdata))), nBoot=1, test=test,show.time = "none")))$table[pred,3]
}


powersim.manyglm_nested <- function(fit,N,coeffs = coef(fit),
                                    pred,nsim = 100,nboot=100,test="score",alpha=0.05,use.design=TRUE,
                                    newdata,cl = makeCluster(detectCores()-1),coeff.null= effect.null(fit,pred)){
  stats <- rep(NA,nsim)
  stats.null <- matrix(ncol=nsim,nrow=nboot)
  old <- Sys.time()
  clusterExport(cl,objects(envir = .GlobalEnv),envir = .GlobalEnv)
  clusterExport(cl,objects(envir = environment()),envir = environment())
  libraries <- clusterEvalQ(cl, {
    library(mvabund)
    library(MASS)
    library(psych)
    library(matrixcalc)
    library(parallel)
  })
  #find nsim by nboot statistics under the null
  stats.null <- parApply(cl,stats.null,c(1,2),MVApowerstat.null)
  #find nsim statistics under the alternative
  stats <- parSapply(cl,stats,MVApowerstat.alt)
  stopCluster(cl)
  #find p_values by comparing each nsim test statistics under alternative to those under the null
  p_values <- rep(NA,length=nsim)
  for(i in c(1:nsim)){
    p_values[i] <- mean(stats[i]<stats.null[,i])
  }
  new <- Sys.time() - old
  #Estimate power as the proportion of times these p-values are less then alpha
  print(c(mean(p_values[!is.na(p_values)]<alpha),new))

  #print(binom.confint(table((stats[!is.na(stats)]>criticalStat))[2],n=nsim,method="wilson")[c(5,6)])
}





