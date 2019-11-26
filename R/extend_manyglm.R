#' Produce 'extended' manyglm object with simulated abundances
#'
#' @description
#' `extend.manyglm` returns a manyglm object with `N` observations and simulated response matrix,
#' that utilises the existing correlation structure of the data.
#'
#' @details
#' `extend.manyglm` takes a manyglm object and returns back an 'extended' manyglm object with `N`
#' observations and a new simulated response matrix. Response abundances are simulated through a Gaussian
#' copula model that utilises a coefficient matrix `coeff`, the specified `manyglm` model and the joint
#' correlation structure exhibited between the response variables. To help with the specification of
#' `coeff`, see `effect.alt` which simplifies this process.
#'
#' Response variables are simulated through a copula model by first extracting Gaussian copular scores
#' as Dunn-Smyth residuals, which are obtained from abundances \eqn{y_{ij}} with marginal distributions \eqn{F_j}
#' which have been specified via the original `manyglm` model `fit`;
#'
#' \deqn{z_{ij} = \Phi^{-1}{F_{j}(y_{ij}^-) + u_{ij} F_{j}(y_{ij})}}
#'
#'  These scores then follow a multivariate Gaussian distribution with zero mean and covariance structure \eqn{\Sigma},
#'
#'  \deqn{z_{ij} ~ N_p(0,\Sigma)}
#'
#' To avoid estimating a large number \eqn{p(p-1)/2} pairwise correlations within \eqn{\Sigma}, factor analysis is utilised
#' with one latent factor variable, which can be interpreted as an unobserved environmental covariate.
#'
#' Thus, in order to simulate new multivariate abundances we simulate new copula scores and back transform them to
#' abundances as \eqn{y_{ij}= {F^*}_j^{-1}(\Phi(z_{ij}))}, where the inputed coefficient matrix `coeff` specifies the
#' effect size within the new marginal distributions \eqn{{F^*}_j}.\\
#'
#' The data frame is also extended in a manner that preserves the original design structure. This is done by first
#' repeating the design matrix until the number of samples exceeds `N`, then randomly removing rows from the last
#' repeated data frame until the number of samples equals `N`.
#'
#' `use.design=False` and `newdata` can be utilised if a different data frame is wanted for simulation.
#'
#' If users are interested in obtaining simulated multivariate abundances alone from a `manyglm` model, these can be
#' extracted as the response matrix from the returned `manyglm` object.
#'
#'@param fit objects of class `manyglm`, typically the result of a call to \link[mvabund]{manyglm}.
#'@param N Number of samples for the returned `manyglm` object to be extended.
#'@param coeffs Coefficient matrix for a `manyglm` object that characterises the size of effects to be simulated.
#'See `effect.alt` for help in producing this matrix. Defaults to the coefficient matrix from the inputed `manyglm`
#'object `coef(fit)`.
#'@param use.design Logical. Wether to utilise the design of the inputed `manyglm` object or the design specified by
#'the data frame `newdata`.
#'@param newdata Data frame of same size as the original data frame from the inputed `manyglm` fit, that specifies
#'a different design of interest.
#'@import mvabund
#'@import psych
#'@import matrixcalc
#'@import MASS
#'@export
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
#'
#'#change sampling design
#' X_new <- X
#' X_new$Treatment[6:7] <- c("B","B")
#' extend.fit <- extend.manyglm(glm.spid,N=20,
#'                              coeffs=effect.mat,use.design = FALSE,newdata=X_new) #not needed to be executed for power estimate
#' powersim.manyglm(glm.spid,N=20,pred="Treatment",coeffs=effect.mat,use.design = FALSE,newdata=X_new)
extend.manyglm <- function(fit,N, coeffs = coef(fit),use.design=TRUE,newdata){
  #extract the dun-smyth residuals
  Z_ds <- residuals(fit)

  #check for infinite values in residuals
  if(sum(apply(Z_ds,c(1,2),is.infinite))>0){
    stop("Infinite Dunn-smyth residuals, consider using negative binomial model or turning data into presence/absences")
  }


  factor      <- suppressMessages(suppressWarnings(fa(Z_ds,warnings=FALSE)))
  unique      <- factor$uniquenesses + 10^(-4)
  loadings    <- factor$loadings
  corr        <- loadings%*%t(loadings) + diag(unique)

  #Ensure correlation matrix is positive definite
  while (is.positive.definite(corr)==FALSE){
    Z_ds <- residuals(fit)
    factor      <- suppressMessages(suppressWarnings(fa(Z_ds,warnings=FALSE)))
    unique      <- factor$uniquenesses
    loadings    <- factor$loadings
    corr        <- loadings%*%t(loadings) + diag(unique)
  }


  #simulate N multivariate normal random variables
  mu      <- rep(0,times=ncol(Z_ds))
  sim     <- suppressWarnings(mvrnorm(N,mu,corr))


  #extend the design matrix appropriately

  #find the number of observations
  Nobs<-nrow(fit$fitted.values)
  #choose the number of data sets to simulate
  if (N%%Nobs > 0){
    nDatasets <- N%/%Nobs+1
  } else{
    nDatasets <- N%/%Nobs
  }



  #Repeat the data set nDataset times
  if (use.design==TRUE){
    if (dim(fit$data[-1])[2]>1){
      Xnew        <- fit$data[rep( 1:Nobs , nDatasets ),]
    }else{
      Xnew        <- data.frame(fit$data[rep( 1:Nobs , nDatasets ),])
      names(Xnew) <- names(fit$data)
    }
  }else{
    if (dim(fit$data[-1])[2]>1){
      Xnew        <- newdata[rep( 1:Nobs , nDatasets ),]
    }else{
      Xnew        <- data.frame(newdata[rep( 1:Nobs , nDatasets ),])
      names(Xnew) <- names(fit$data)
    }
  }


  #remove the excess values from the model frame
  if (N%%Nobs > 0){
    whichRemove <- sample(Nobs,(Nobs- N%%Nobs)) + Nobs*N%/%Nobs
    Xnew        <- Xnew[-whichRemove,]
  } else {
    Xnew        <- Xnew
  }

  #fix to column if model frame becomes a vector
  if (dim(fit$data[-1])[2]>1){
    Xnew        <- Xnew
  }else{
    Xnew        <- data.frame(Xnew)
    names(Xnew) <- names(fit$data)
  }


  #Repeat family nuisance parameters

  if (use.design==TRUE){
    design.matrix = model.matrix(fit$formula,data=fit$data)
  }else{
    design.matrix = model.matrix(fit$formula,data=newdata)
  }


  if (fit$family=="negative.binomial"){
    Xlam      = design.matrix[rep( 1:Nobs , nDatasets ),]
    lambdanew = negative.binomial(theta=1)$linkinv(Xlam%*% coeffs)
    size      <- fit$theta
    #remove excess values
    if (N%%Nobs > 0){
      lambdanew   <- lambdanew[-whichRemove,]
    } else {
      lambdanew   <- lambdanew
    }
  }else{if(fit$family=="poisson"){
    Xlam      = design.matrix[rep( 1:Nobs , nDatasets ),]
    lambdanew = poisson(link="log")$linkinv(Xlam%*% coeffs)
    #remove excess values
    if (N%%Nobs > 0){
      lambdanew   <- lambdanew[-whichRemove,]
    } else {
      lambdanew   <- lambdanew
    }
  }else{if(fit$family=="binomial(link=logit)"){
    Xlam      = design.matrix[rep( 1:Nobs , nDatasets ),]
    probnew = binomial(link="logit")$linkinv(Xlam%*% coeffs)
    #remove excess values
    if (N%%Nobs > 0){
      probnew   <- probnew[-whichRemove,]
    } else {
      probnew   <- probnew
    }
  }else{if(fit$family=="binomial(link=cloglog)"){
    Xlam      = design.matrix[rep( 1:Nobs , nDatasets ),]
    probnew = binomial(link="cloglog")$linkinv(Xlam%*% coeffs)
    #remove excess values
    if (N%%Nobs > 0){
      probnew   <- probnew[-whichRemove,]
    } else {
      probnew   <- probnew
    }
  }else{
    stop("This package does not currently support this family")
  }
  }
  }
  }





  #turn these simulated MVN residuals back to abundances
  ys <- matrix(nrow=N,ncol=ncol(fit$y))

  if (fit$family=="negative.binomial"){
    for (i in c(1:ncol(fit$y))){
      ys[,i]    <- qnbinom(pnorm(sim[,i]),size=size[i],mu=lambdanew[,i])
    }
  }else {if(fit$family=="poisson"){
    for (i in c(1:ncol(fit$y))){
      ys[,i]    <- qpois(pnorm(sim[,i]),lambda = lambdanew[,i])
    }
  }else{
    for (i in c(1:ncol(fit$y))){
      ys[,i]    <- qbinom(pnorm(sim[,i]),size=1,prob = probnew[,i])
    }
  }
  }


  #voila we have simulated multivariate abundance data!!
  #give the column names back to our simulated values
  colnames(ys) <- colnames(fit$y)



  #now join back our data frame
  data.sim <- data.frame(ys,Xnew,row.names = NULL)
  ##CHange below line in legit code to be able
  #colnames(data.sim)[[dim(data.sim)[2]]] <- names(fit$data)[1]
  ys <- as.matrix(ys)


  #so now we have our data frame and we want to refit the model:
  fit$call[[4]] <- quote(data.sim)
  fit$call[[2]][[2]] <- quote(ys)


  return(eval(fit$call))
}
