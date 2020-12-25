#' Simulate or extend multivariate abundance data
#'
#' @description
#' \code{extend} returns a simulated response matrix or a \code{\link[mvabund]{manyglm}} object with \code{N} observations
#' and simulated response matrix that utilises the existing correlation structure of the data.
#'
#' @details
#' \code{extend} takes a \code{\link[ecoCopula]{cord}} object and returns a new simulated response matrix or an "extended" \code{\link[mvabund]{manyglm}} object
#' with \code{N} observations and the new simulated response matrix. Response abundances are simulated through a Gaussian
#' copula model that utilises a coefficient matrix \code{coeffs}, the specified \code{cord} model and the joint
#' correlation structure exhibited between the response variables. To help with the specification of
#' \code{coeffs}, see \code{\link{effect_alt}} which simplifies this process.
#'
#' Response variables are simulated through a copula model by first extracting Gaussian copular scores
#' as Dunn-Smyth residuals, which are obtained from abundances \eqn{y_{ij}} with marginal distributions \eqn{F_j}
#' which have been specified via the original \code{manyglm} model \code{fit.manyglm};
#'
#' \deqn{z_{ij} = \Phi^{-1}{F_{j}(y_{ij}^-) + u_{ij} F_{j}(y_{ij})}}
#'
#'  These scores then follow a multivariate Gaussian distribution with zero mean and covariance structure \eqn{\Sigma},
#'
#'  \deqn{z_{ij} ~ N_p(0,\Sigma)}
#'
#' To avoid estimating a large number \eqn{p(p-1)/2} pairwise correlations within \eqn{\Sigma}, factor analysis is utilised
#' with two latent factor variables, which can be interpreted as an unobserved environmental covariate.
#'
#' Thus, in order to simulate new multivariate abundances we simulate new copula scores and back transform them to
#' abundances as \eqn{y_{ij}= {F^*}_j^{-1}(\Phi(z_{ij}))}, where the inputed coefficient matrix \code{coeffs} specifies the
#' effect size within the new marginal distributions \eqn{{F^*}_j}.
#'
#' The data frame is also extended in a manner that preserves the original design structure. This is done by first
#' repeating the design matrix until the number of samples exceeds \code{N}, then randomly removing rows from the last
#' repeated data frame until the number of samples equals \code{N}. Alternatively, a balanced design structure can be
#' obtained by specifying the number of replicates.
#'
#' \code{newdata} can be utilised if a different data frame is wanted for simulation.
#'
#' If users are interested in obtaining a \code{manyglm} model, \code{do.fit=TRUE} can be used to obtain a \code{\link[mvabund]{manyglm}}
#' object from the simulated responses.
#'
#' @param object objects of class \code{cord}, typically the result of a call to \code{\link[ecoCopula]{cord}}.
#' @param N Number of samples to be extended. Defaults to the number of observations in the original sample.
#' @param coeffs Coefficient matrix for a \code{\link[mvabund]{manyglm}} object that characterises the size of effects to be simulated.
#' See \code{\link{effect_alt}} for help in producing this matrix. Defaults to the coefficient matrix from the inputed \code{\link[mvabund]{manyglm}}
#' object \code{coef(object$obj)}.
#' @param newdata Data frame of same size as the original data frame from the inputed \code{manyglm} fit, that specifies
#' a different design of interest.
#' @param n_replicate Number of unique replicates of the original data frame. Defaults to \code{NULL}, overwrites \code{N} if specified.
#' @param do.fit Logical. If \code{TRUE}, fits a \code{\link[mvabund]{manyglm}} object from the simulated data. Defaults to \code{FALSE}.
#' @param seed Random number seed, defaults to a random seed number.
#' @return Simulated data or \code{manyglm} object.
#' @seealso \code{\link{effect_alt}}
#' @import ecoCopula
#' @import mvabund
#' @importFrom stats simulate
#' @importFrom stats reshape
#' @importFrom stats coef
#' @importFrom stats terms
#' @importFrom stats formula
#' @rdname extend
#' @examples
#' library(ecoCopula)
#' library(mvabund)
#' data(spider)
#' spiddat = mvabund(spider$abund)
#' X = data.frame(spider$x)
#'
#' # Specify increasers and decreasers
#' increasers = c("Alopacce", "Arctlute", "Arctperi", "Pardnigr", "Pardpull")
#' decreasers = c("Alopcune", "Alopfabr", "Zoraspin")
#' 
#' # Simulate data
#' fit.glm = manyglm(spiddat~1, family="negative.binomial")
#' fit.cord = cord(fit.glm)
#' simData = extend(fit.cord)
#'
#' # Simulate data with N=20
#' fit.glm = manyglm(spiddat~soil.dry, family="negative.binomial", data=X)
#' fit.cord = cord(fit.glm)
#' simData = extend(fit.cord, N=20)
#'
#' # Obtain a manyglm fit from simulated data with N=10 and effect_size=1.5
#' X$Treatment = rep(c("A","B","C","D"),each=7)
#' fit_factors.glm = manyglm(spiddat~Treatment, family="negative.binomial", data=X)
#' effect_mat = effect_alt(fit_factors.glm, effect_size=1.5,
#'      increasers, decreasers, term="Treatment")
#' fit_factors.cord = cord(fit_factors.glm)
#' newFit.glm = extend(fit_factors.cord, N=10,
#'      coeffs=effect_mat, do.fit=TRUE)
#'
#' # Change sampling design
#' X_new = X
#' X_new$Treatment[6:7] = c("B","B")
#' simData = extend(fit_factors.cord, N=NULL,
#'    coeffs=effect_mat, newdata=X_new, n_replicate=5)
#' @export

extend.cord = function(object, N=nrow(object$obj$data), coeffs=coef(object$obj), newdata=NULL, 
  n_replicate=NULL, do.fit=FALSE, seed=NULL) {

  #find the number of observations
  Nobs = nrow(object$obj$fitted.values)
  N = compute_N(object, N, n_replicate)
  #choose the number of data sets to simulate
  nDatasets = compute_nDatasets(N, Nobs)

  Xnew = get_Xnew(object, n_replicate, newdata, Nobs, nDatasets, N)
  Xnew = remove_excess_rows(N, Nobs, Xnew)

  object$obj$coefficients = coeffs
  Ynew = simulate(object=object, nsim=1, seed=seed, newdata=Xnew)

  extended_data = data.frame(Ynew, Xnew, row.names = NULL)
  if (do.fit == TRUE) {
    out = get_new_fit(object, Ynew, extended_data, Xnew)
  } else {
    extended_data = drop_ones(extended_data, Xnew)
    out = list(data = extended_data)
  }
  return (out)
}

get_Xnew = function(object, n_replicate, newdata, Nobs, nDatasets, N) {
  if (formula(object$obj)[-2] == ~1) {
    Xnew = data.frame(array(1, N))
    names(Xnew) = "ones.intercept"
  } else {
    if (!is.null(n_replicate)) {
      if (is.null(newdata)) {
        Xnew = rep_unique_Xnew(object$obj$data, n_replicate)
      } else {
        Xnew = rep_unique_Xnew(newdata, n_replicate)
      }
    } else {
      #extend the design matrix appropriately

      # Repeat the data set nDataset times
      if (is.null(newdata)) {
        Xnew = rep_Xnew(object$obj$data, Nobs, nDatasets)
      } else {
        Xnew = rep_Xnew(newdata, Nobs, nDatasets)
      }
    }
  }

  Xnew = Xnew[!colnames(Xnew) %in% colnames(object$obj$y)]
  return (Xnew)
}

rep_unique_Xnew = function(frame, n_replicate) {
  if (ncol(frame) > 1) {
    Xnew = do.call("rbind", replicate(n_replicate, unique(frame), simplify = FALSE))
  } else {
    name = colnames(frame)
    colnames(frame) = "V1"
    Xnew = data.frame(rep(unique(frame), each=n_replicate))

    Xnew = vectorize_Xnew(Xnew, name)
  }
  return (Xnew)
}

rep_Xnew = function(frame, Nobs, nDatasets) {
  if (ncol(frame) > 1) {
    Xnew = frame[rep( 1:Nobs , nDatasets ),]
  } else {
    name = colnames(frame)
    colnames(frame) = "V1"
    Xnew = data.frame(frame[rep( 1:Nobs , nDatasets ),])

    Xnew = vectorize_Xnew(Xnew, name)
  }
  return (Xnew)
}

vectorize_Xnew = function(Xnew, name) {
  if (ncol(Xnew) > 1) {
    Xnew = reshape(Xnew, direction = "long", varying=1:ncol(Xnew))["V1"]
  }
  colnames(Xnew) = name
  return (Xnew)
}

get_new_fit = function(object, Ynew, extended_data, Xnew) {
  Ynew <- as.matrix(Ynew)
  # data = extended_data
  object$obj$call[[2]][[2]] <- quote(Ynew)

  if (!"ones.intercept" %in% names(Xnew)) {
    object$obj$call[[4]] <- quote(extended_data)
  }

  new_fit = eval(object$obj$call)
  return (new_fit)
}

compute_N = function(object, N, n_replicate) {
  if (!is.null(n_replicate)) {
    N_cal = nrow(unique(object$obj$data)) * n_replicate
    if (!is.null(N) && N_cal!=N) {
      warning("N has been overwritten by n_replicate.")        
    }
    N = N_cal
  } else {
    N = N
  }
  return (N)
}

compute_nDatasets = function(N, Nobs) {
  if (N%%Nobs > 0) {
    nDatasets = N%/%Nobs+1
  } else {
    nDatasets = N%/%Nobs
  }
  return (nDatasets)
}

remove_excess_rows = function(N, Nobs, frame) {
  if (N%%Nobs > 0 && N!=nrow(frame)) {
    whichRemove = sample(Nobs,(Nobs- N%%Nobs)) + Nobs*N%/%Nobs
    frame   = frame[-whichRemove,]
  } else {
    frame   = frame
  }
  return (frame)
}

drop_ones = function(extended_data, Xnew) {
  if ("ones.intercept" %in% names(Xnew)) {
    extended_data = extended_data[, !(names(extended_data) %in% "ones.intercept")]
  } else {
    extended_data = extended_data
  }
  return (extended_data)
}
