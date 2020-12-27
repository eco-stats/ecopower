#' Specify null effects for multivariate abundance data
#'
#' @description
#' \code{effect_null} returns a coefficient matrix to be parsed to \code{\link{powersim}} by default
#' to specify a null effect.
#'
#' @details
#' \code{effect_null} produces a coefficient matrix with a null effect that is specified by setting the parameter
#' estimates of a predictor of interest \code{term} to 0. This function is used by default in \code{\link{powersim}}.
#' Note that intercept values are parameterised as in \code{\link{effect_alt}}.
#' @param object objects of class \code{manyglm}, typically the result of a call to \code{\link[mvabund]{manyglm}}.
#' @param term Name of predictor of interest in quotes.
#' @return A coefficient matrix with the null effect.
#' @seealso \code{\link{effect_alt}}, \code{\link{powersim}}
#' @import mvabund
#' @importFrom stats update
#' @importFrom stats drop.terms
#' @importFrom stats binomial
#' @importFrom stats poisson
#' @rdname effect_null
#' @examples
#' library(mvabund)
#' data(spider)
#' spiddat = mvabund(spider$abund)
#' X = data.frame(spider$x)
#'
#' # Find null effect size for continuous predictor
#' spid.glm = manyglm(spiddat~soil.dry, family="negative.binomial", data=X)
#' coeff_null = effect_null(spid.glm, term="soil.dry")
#' @export

effect_null.manyglm = function(object, term) {

  coeff = object$coefficients
  rowNames = row.names(coeff)[startsWith(row.names(coeff), term)]
  nRow = length(rowNames)
  nCol = ncol(coeff)

  # Specify term based on class
  if (class(object$data[,term]) %in% c("integer", "numeric")) {
    coeff[term,] = 0
  } else if (class(object$data[,term]) %in% c("factor", "character")) {
    # Set mean abundance as referance group/intercept to avoid negative comparisons of null abundances across groups
    #coeff[1,]<- inv.func(apply(fit$y,2,mean))
    fit_null = get_object0(object, term, NULL)
    coeff_null = fit_null$coefficients

    term_matrix = matrix(data=0, nrow=nRow, ncol=nCol)
    rownames(term_matrix) = rowNames

    coeff = rbind(coeff_null, term_matrix)

  } else {
    stop("This function does not recognise this type of term")
  }

  return(coeff)
}
