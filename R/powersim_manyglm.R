#' Provide power estimates for multivariate abundance models
#'
#' @description
#' `powersim.manyglm` returns a power estimate for a `manyglm` object for a given sample size `N`
#' and effect size of interest.
#'
#' @details
#' `powersim.manyglm` takes a `manyglm` object, sample size `N` and coefficient matrix `coef` which
#' specifies an effect size of interest and returns a power estimate.
#'
#' The power estimate is obtained by first parsing the inputed `manyglm` object into \link[ecopower]{extend.manyglm},
#' `nsim` times with an effect size specified by `coeff`. Next, the `manyglm` object is parsed into
#' `extend.manyglm` an additional `nsim` times with a null effect, which is defined by default by
#'  \link[ecopower]{effect.null}. This effectively simulates `nsim` `manyglm` models under both the null
#'  and alternative hypothesis.
#'
#'  For each simulated `manyglm` object a test statistic `test` is obtained. A critical test statistic
#'  is then obtained as the upper 1 - `alpha` quantile of simulated test statistics under the null
#'  hypothesis. Power is then estimated as the proportion of times the test statistics simulated under
#'  the alternative hypothesis exceed the critical test statistic under the null.
#'
#'  To improve computation time, simulations are computed in parrellel using the 'socket' approach, which
#'  be default uses all available cores for clustering. If the function returns an error relating to
#'  clustering or nodes, it is recommended to use 1 less than the number of available cores for your
#'  machine; makeCluster(detectCores()-1).
#'@param fit objects of class `manyglm`, typically the result of a call to \link[mvabund]{manyglm}.
#'@param N Number of samples for power estimate.
#'@param coeffs Coefficient matrix for a `manyglm` object that characterises the size of effects to be simulated.
#'See `effect.alt` for help in producing this matrix. Defaults to the coefficient matrix from the inputed `manyglm`
#'object `coef(fit)`.
#'@param pred Name of predictor of interest in quotes.
#'@param nsim Number of simulations for power estimate to be based upon. Defaults to 1000.
#'@param test Test statistic for power estimate to based upon. Defaults to "score", however "wald" is also allowed.
#'@param alpha Type I error rate for power estimate, defaults to 0.05.
#'@param use.design Logical. Wether to utilise the design of the inputed `manyglm` object or the design specified by
#'the data frame `newdata`.
#'@param newdata Data frame of same size as the original data frame from the inputed `manyglm` fit, that specifies
#'a different design of interest.
#'@param cl Number of clusters for parrelel computing. Defaults to the total number of clusters available on the
#'machine.
#'@param coeff.null Coefficient matrix under the null hypothesis. Defaults to being specified by \link[ecopower]{effect.null}.
#'@examples
#'library(mvabund)
#'data(spider)
#' spiddat <- mvabund(spider$abund)
#' X <- data.frame(spider$x)
#'
#' #Specify 'increasers' and 'decreasers'
#' increasers <- c("Alopacce","Arctlute" ,"Arctperi","Pardnigr", "Pardpull")
#' decreasers <- c("Alopcune","Alopfabr" ,"Zoraspin")
#'
#' #Find power for continuous predictor, N=20 and effect.size=3
#' glm.spid <- manyglm(spiddat~soil.dry, family="negative.binomial",data=X)
#' effect.mat <- effect.alt(glm.spid,effect.size=3,pred="soil.dry",increasers,decreasers)
#' extend.fit <- extend.manyglm(glm.spid,N=10,
#'                              coeffs=effect.mat) #not needed to be executed for power estimate
#' powersim.manyglm(glm.spid,N=20,pred="soil.dry",coeffs=effect.mat,cl=makeCluster(1))
#'
#' #Find power for categorical predictor with 4 levels, N=10, effect.size=1.5
#' X$Treatment <- rep(c("A","B","C","D"),each=7)
#' glm.spid <- manyglm(spiddat~Treatment, family="negative.binomial",data=X)
#' effect.mat <- effect.alt(glm.spid,effect.size=1.5,pred="Treatment",increasers,decreasers)
#' extend.fit <- extend.manyglm(glm.spid,N=20,
#'                              coeffs=effect.mat) #not needed to be executed for power estimate
#' powersim.manyglm(glm.spid,N=20,pred="Treatment",coeffs=effect.mat,cl=makeCluster(1))
#'
#' #change effect size parameterisation
#' effect.mat <- effect.alt(glm.spid,effect.size=1.5,
#'                          pred="Treatment",increasers,decreasers,
#'                          K=c(3,1,2),OrderedLevels = FALSE)
#' powersim.manyglm(glm.spid,N=20,pred="Treatment",
#'                  coeffs=effect.mat,use.design = FALSE,newdata=X_new,cl=makeCluster(1))
#'
#'#change sampling design
#' X_new <- X
#' X_new$Treatment[6:7] <- c("B","B")
#' extend.fit <- extend.manyglm(glm.spid,N=20,
#'                              coeffs=effect.mat,use.design = FALSE,newdata=X_new) #not needed to be executed for power estimate
#' powersim.manyglm(glm.spid,N=20,pred="Treatment",
#' coeffs=effect.mat,use.design = FALSE,newdata=X_new,cl=makeCluster(1))
powersim.manyglm <- function(fit,N,coeffs = coef(fit),
                             pred,nsim = 1000,test="score",alpha=0.05,use.design=TRUE,
                             newdata=NULL,cl = makeCluster(detectCores()),coeff.null= effect.null(fit,pred)){
  coeff.null <- effect.null(fit,pred)
  stats.null <- stats <- rep(NA,nsim)
  old <- Sys.time()
  #find critical value under the null
  clusterExport(cl,objects(envir = .GlobalEnv),envir = .GlobalEnv)
  clusterExport(cl,objects(envir = environment()),envir = environment())
  libraries <- clusterEvalQ(cl, {
    library(mvabund)
    library(MASS)
    library(psych)
    library(matrixcalc)
    library(parallel)
  })
  stats.null <- parSapply(cl,stats.null,MVApowerstat.null)
  criticalStat <- quantile(unlist(stats.null[!is.na(unlist(stats.null))]),1-alpha,na.rm=TRUE)
  criticalStat_lower <- quantile(unlist(stats.null[!is.na(unlist(stats.null))]),1-alpha+0.025,na.rm=TRUE)
  criticalStat_upper <- quantile(unlist(stats.null[!is.na(unlist(stats.null))]),1-alpha-0.025,na.rm=TRUE)
  #obseve the proportion of times our test statistics exceed this value
  stats <- unlist(parSapply(cl,stats,MVApowerstat.alt))
  stopCluster(cl)
  new <- Sys.time() - old
  #print(c(mean(stats[!is.na(stats)]>criticalStat),mean(stats[!is.na(stats)]>criticalStat_lower),mean(stats[!is.na(stats)]>criticalStat_upper),new))
  p <- c(mean(stats[!is.na(stats)]>criticalStat),new)
  names(p) <- c("Power", "Comp time")
  print(p)
  #print(binom.confint(table((stats[!is.na(stats)]>criticalStat))[2],n=nsim,method="wilson")[c(5,6)])
}



