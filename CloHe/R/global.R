#' Print a CloHe S4 class to standard output.
#'
#' @param x a CloHe object: a \code{\linkS4class{GaussianMutSigmatModel}} or
#' a \code{\linkS4class{FuncSpectraModel}}.
#' @param ... further arguments passed to or from other methods
#'
#' @return NULL. Prints to standard out.
#'
#' @name print
#' @rdname print-methods
#' @docType methods
#' @seealso \code{\link{print}}
#'
NULL

#' Extract parts of a CloHe S4 class
#'
#' @param x object from which to extract element(s) or in which to replace element(s).
#' @param i the name of the element we want to extract or replace.
#' @param j if the element designing by i is complex, j specifying elements to extract or replace.
#' @param drop For matrices and arrays.  If TRUE the result is coerced to the lowest
#' possible dimension (see the examples).  This only works for extracting elements,
#' not for the replacement.  See drop for further details.
#' @param value	typically an array-like R object of a similar class as the element
#' of x we want to replace.
#' @name [
#' @aliases Extract
#' @docType methods
#' @rdname extract-methods
#'
NULL

#' Show description of a CloHe S4 class to standard output.
#'
#' @param object a MixAll object: a \code{\linkS4class{GaussianMutSigmatModel}},
#' or a \code{\linkS4class{FuncSpectraModel}}.
#'
#' @return NULL. Prints to standard out.
#'
#' @name show
#' @docType methods
#' @rdname show-methods
#' @seealso \code{\link{show}}
NULL


#' Produce summary of a CloHe S4 class.
#'
#' @param object any cluster model deriving from a \code{\linkS4class{ICloheModel}} object.
#' @param ... further arguments passed to or from other methods
#'
#' @return NULL. Summaries to standard out.
#'
#' @name summary
#' @docType methods
#' @rdname summary-methods
#'
# @seealso \code{\link{base::summary}}
#'
NULL

