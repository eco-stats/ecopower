#' Specify null effects for multivariate abundance data
#'
#' @description
#' `effect.null` returns a coefficient matrix to be parsed to `powersim.manyglm` by default
#' to specify a null effect.
#'
#' @details
#' `effect.null` produces a coefficient matrix with a null effect that is specified by setting the parameter
#' estimates of a predictor of interest `pred` to 0. This function is used by default in \link[ecopower]{powersim.manyglm}.
#' Note that intercept values are parameterised as in \link[ecopower]{effect.alt}.
#'@param fit objects of class manyglm, typically the result of a call to \link[mvabund]{manyglm}.
#'@param pred Name of predictor of interest in quotes.
#'@examples
#'library(mvabund)
#'data(spider)
#' spiddat <- mvabund(spider$abund)
#' X <- data.frame(spider$x)
#'
#' #Find null effect size for continuous predictor
#' glm.spid <- manyglm(spiddat~soil.dry, family="negative.binomial",data=X)
#' coeff.null <- effect.null(glm.spid,pred="soil.dry")
effect.null <- function(fit,pred){
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
  }else{if(class(fit$data[,pred])=="factor"|class(fit$data[,pred])=="character"){
    #set mean abundance as referance group/intercept to avoid negative comparisons of null abundances across groups
    #coeff[1,]<- inv.func(apply(fit$y,2,mean))
    fit$call[[2]] <- gsub(pred,1,fit$call[2])
    fit_noeffect <- eval(fit$call)
    coeff[1,]<- inv.func(apply(fit$y,2,mean))
    for (i in c(1:length(row.names(coeff)[startsWith(row.names(coeff),pred)]))){
      coeff[row.names(coeff)[startsWith(row.names(coeff),pred)][i],] <- 0
    }
  }else{
    stop("This function does not recognise this type of predictor")
  }
  }
  return(coeff)
}
