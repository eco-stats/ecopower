#' @name effect.alt
#'
#' @title Specify multivariate effect sizes.
#'
#'
#' @description
#' `effect.alt` returns a coefficient matrix to be parsed to `extend.manyglm` and `powersim.manyglm`
#' to specify an effect size of interest.
#'
#' @details
#' `effect.alt` helps users to create interpretable multivariate effect sizes to be parsed into `extend.manyglm`
#' and `powersim.manyglm`, so that researchers can investigate the relaitionship between effect size, power and
#' sample size in a complicated multivariate abundance setting.
#'
#' `effect.alt` creates an effect of size `log(effect.size)` for a predictor of interest (`pred`), for responses
#'  who have been specified to increase (`increasers`) and `-log(effect.size)` for responses who have been secified
#'  to decrease (`decreasers`). Responses that have not been specified in the `increasers` or `decreasers` vectors
#'  are specified to have no effect with a coefficient of 0. The effect has been logged to make the effect size
#'  interpretable within the coefficient matrix.
#'
#'  For poisson regression `family=poisson()` and negative binomial regression `family="negative.binomial"`the effect
#'  size is interpreted for a categorical variable as the multiplicative change in mean abundance in the treatment
#'  group relative to the control group, whilst for a continuous variable it is interpreted as the multiplicative
#'  change in abundance for a 1 unit increase in the predictor of interest.
#'
#'  For logit regression `family=binomial("logit")` the effect size is interpreted as an odds ratio. For a categorical
#'  variable this is the change in odds of obtaining outcome `1` when being in the treatment group relative to the
#'  control group. Whilst for continuous variables, this is interpreted as the change in odds of obtaining outcome `1`
#'  with a 1 unit increase in the predictor of interest.
#'
#'  For cloglog regression `family=binomial("cloglog")` the effect size is interpreted similarly to poisson and
#'  negative binomial regression. For a categorical variable it is interpreted as the multiplicative change in
#'  the mean of the underlying count in the treatment group relative to the control. Whilst for a continuous
#'  variable it is interpreted as the multiplicative change in the mean of the underlying count for a 1 unit
#'  increase in the predictor of interest.
#'
#'  For categorical variables, the intercept is also changed to be the group mean intercept by taking the
#'  intercept of a model without the categorical predictor of interest. This is done to avoid messy comparisons
#'  of null control groups.
#'
#'  For categorical variables with more than two levels, effect size is changed to `effect.size^K[i]` where K defaults
#'  to  be `c(1,2,...,nlevels - 1)`, where `nlevels` are the number of levels of the categorical variable and is
#'  specified along the order of the levels. To change this, specify `OrderedLevels = FALSE` and provide a vector
#'  `K` with length of `nlevels - 1`. To change the control group, this must be done prior to specifying the
#'  `manyglm` object using `relevel` (which can also change the order of the levels).
#'
#'  Note that if the predictor of interest is a categorical variable it must be classed either as a factor or
#'  character otherwise results may be misleading.
#'
#'@param fit objects of class manyglm, typically the result of a call to \link[mvabund]{manyglm}.
#'@param effect.size An effect size of interest, see description for interpretation.
#'@param increasers A vector list of responses which increase relative to the control group/intercept.
#'@param decreasers A vector list of responses which decrease relative to the control group/intercept.
#'@param pred Name of predictor of interest in quotes.
#'@param OrderedLevels Logical. Whether or not to default the effect size to increase it's exponent by the
#'order of factor variables. See description for further details.
#'@param K If `OrderedLevels=FALSE`, the exponent of the effect.size for each level of a factor variable.
#'@import mvabund
#'@export
#'@examples
#' library(mvabund)
#' data(spider)
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
#' powersim.manyglm(glm.spid,N=20,pred="soil.dry",coeffs=effect.mat)
#'
#' #Find power for categorical predictor with 4 levels, N=10, effect.size=1.5
#' X$Treatment <- rep(c("A","B","C","D"),each=7)
#' glm.spid <- manyglm(spiddat~Treatment, family="negative.binomial",data=X)
#' effect.mat <- effect.alt(glm.spid,effect.size=1.5,pred="Treatment",increasers,decreasers)
#' extend.fit <- extend.manyglm(glm.spid,N=20,
#'                              coeffs=effect.mat) #not needed to be executed for power estimate
#' powersim.manyglm(glm.spid,N=20,pred="Treatment",coeffs=effect.mat)
#'
#' #change effect size parameterisation
#' effect.mat <- effect.alt(glm.spid,effect.size=1.5,
#'                          pred="Treatment",increasers,decreasers,
#'                          K=c(3,1,2),OrderedLevels = FALSE)
#' powersim.manyglm(glm.spid,N=20,pred="Treatment",
#'                  coeffs=effect.mat)
effect.alt <- function(fit,effect.size,increasers,decreasers,pred,OrderedLevels = TRUE,K){
  coeff  <- fit$coefficients

  #create a function with the appropriate link
  if (fit$family=="negative.binomial"){
    inv.func = function(x){negative.binomial(theta=1)$linkfun(x)}
  }else{if(fit$family=="poisson"){
    inv.func = function(x){poisson(link="log")$linkfun(x)}
  }else{if(fit$family=="binomial(link=logit)"){
    inv.func = function(x){binomial(link="logit")$linkfun(x)}
  }else{if(fit$family=="binomial(link=cloglog)"){
    inv.func = function(x){binomial(link="cloglog")$linkfun(x)}
  }else{
    stop("This package does not currently support this family")
  }
  }
  }
  }



  #specify pred based on who
  if(class(fit$data[,pred])=="integer"|class(fit$data[,pred])=="numeric"){
    coeff[pred,]<- 0
    coeff[pred,][increasers] <- log(effect.size)
    coeff[pred,][decreasers] <- -log(effect.size)
  }else{if(class(fit$data[,pred])=="factor"|class(fit$data[,pred])=="character"){
    #set mean abundance as referance group to avoid negative comparisons of null abundances across groups
    fit$call[[2]] <- gsub(pred,1,fit$call[2])
    fit_noeffect <- eval(fit$call)
    #coeff[1,] <- fit_noeffect$coefficients[1,]
    coeff[1,]<- inv.func(apply(fit$y,2,mean))
    if(OrderedLevels == TRUE){
      K=seq(length(row.names(coeff)[startsWith(row.names(coeff),pred)]))
    }else{if(length(K)==length(row.names(coeff)[startsWith(row.names(coeff),pred)])){
      K=K
    }else{
      stop("Length of K is not equal to the number of pred levels - 1")
    }
    }
    for (i in c(1:length(row.names(coeff)[startsWith(row.names(coeff),pred)]))){
      coeff[row.names(coeff)[startsWith(row.names(coeff),pred)][i],] <- 0
      coeff[row.names(coeff)[startsWith(row.names(coeff),pred)][i],][increasers] <- log(effect.size^K[i])
      coeff[row.names(coeff)[startsWith(row.names(coeff),pred)][i],][decreasers] <- -log(effect.size^K[i])
    }
  }else{
    stop("This function does not recognise this type of predictor")
  }
  }
  return(coeff)
}
