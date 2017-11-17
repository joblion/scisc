#-----------------------------------------------------------------------
#     Copyright (C) 2012-2016  Inria
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as
#    published by the Free Software Foundation; either version 2 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with this program; if not, write to the
#    Free Software Foundation, Inc.,
#    59 Temple Place,
#    Suite 330,
#    Boston, MA 02111-1307
#    USA
#
#    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
#
#-----------------------------------------------------------------------
NULL

learnGaussian <- function(data)
{
  classNumbers <- as.integer(levels(factor(data$labels)))
  if (length(classNumbers) < 2)
  { stop("in learnGaussian, not enough class")}
  # create vector of results
  res <- lapply( rep("GaussianMutSigmatModel", length(classNumbers)), new)
  # launch computations
  resLearn <- .Call( "launchGaussian_mut_sigmat"
                   , data$labels, data$times
                   , data$xb, data$xg, data$xr, data$xi
                   , data$clouds
                   , res
                   , PACKAGE = "CloHe")
  resLearn
}


#-----------------------------------------------------------------------
#' Definition of the [\code{\linkS4class{GaussianMutSigmatModel}}] class
#'
#' This class defines a Gaussian model for the spectrum with mean
#' and variance varying along the time. It inherits
#' from [\code{\linkS4class{ICloHeModel}}].
#'
#' @slot mut A list with the vectors of the means.
#' @slot sigmat A list with the variances matrices.
#'
#' @seealso [\code{\linkS4class{ICloHeModel}}] class
#'
#' @examples
#' getSlots("GaussianMutSigmatModel")
#'
#' @author Serge Iovleff
#'
#' @name GaussianMutSigmatModel
#' @rdname GaussianMutSigmatModel-class
#' @aliases GaussianMutSigmatModel-class
#'
setClass(
    Class="GaussianMutSigmatModel",
    representation( mut = "list"
                  , sigmat = "list"
                  ),
    contains=c("ICloHeModel"),
    validity=function(object)
    {
      if (length(object@mut)!=length(object@sigmat))
      {stop("mut and sigmat must have the same length.")}
      return(TRUE)
    }
)

#' Initialize an instance of a MixAll S4 class.
#'
#' Initialization method of the [\code{\linkS4class{GaussianMutSigmatModel}}] class.
#' Used internally in the 'MixAll' package.
#'
#' @rdname initialize-methods
#' @keywords internal
setMethod(
    f="initialize",
    signature=c("GaussianMutSigmatModel"),
    definition=function(.Object)
    {
      # validate
      .Object <- callNextMethod(.Object)
      validObject(.Object)
      return(.Object)
    }
)

#' @rdname print-methods
#' @aliases print GaussianMutSigmatModel-method
#'
setMethod(
    f="print",
    signature=c("GaussianMutSigmatModel"),
    function(x,...){
      cat("****************************************\n")
      callNextMethod();
      cat("****************************************\n")
      cat("* mut    =\n",)
      print( (matrix(unlist(x@mut), nrow=4, byrow = F)), quote=FALSE)
      cat("****************************************\n")
      cat("* sigmat =\n")
      print( (matrix(unlist(x@sigmat), nrow=4, byrow = F)), quote=FALSE)
      cat("****************************************\n")
    }
)

#' @rdname show-methods
#' @aliases show-GaussianMutSigmatModel,GaussianMutSigmatModel,GaussianMutSigmatModel-method
setMethod(
    f="show",
    signature=c("GaussianMutSigmatModel"),
    function(object)
    {
      cat("****************************************\n")
      callNextMethod();
      cat("****************************************\n")
      if(length(object@mut) != 0)
      {
        ncolShow <- min(10,length(object@mut));
        cat("* mut (limited to 10 dates) =\n")
        print( format(matrix(unlist(object@mut), nrow = 4, byrow = F)[,1:ncolShow]), quote = FALSE)
      }
      cat("* ... ...\n")
      cat("****************************************\n")
      if(length(object@sigmat) != 0)
      {
        ncolShow <- min(10,length(object@sigmat))*4;
        cat("* sigmat (limited to 10 dates) =\n")
        print( format(matrix(unlist(object@sigmat), nrow=4, byrow=F)[,1:ncolShow]), quote=FALSE)
      }
      cat("****************************************\n")
    }
)

#' @rdname summary-methods
#' @aliases summary summary,GaussianMutSigmatModel-method
#'
setMethod(
    f = "summary",
    signature=c("GaussianMutSigmatModel"),
    function(object, ...)
    {
      cat("**************************************************************\n")
      callNextMethod()
      cat("**************************************************************\n")
    }
)

