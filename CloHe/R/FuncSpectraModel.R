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

#-----------------------------------------------------------------------
#' Create an instance of the [\code{\linkS4class{FuncSpectraModel}}] class
#'
#' This function ajust the parameters of a functionnal gaussian curve model.
#'
#' @param data list containing the data.
#' @param kernelName [\code{string}] with the kernel to use for the covariance matrix.
#' Available kernels are "gaussian", "exponential", "rationalQuadratic".
#' Default is "gaussian".
#' @param width the width of the kernel to use. Default is 50.
#' @param dim  [\code{\link{integer}}] number of base to use. Default is 5.
#' @param basisName [\code{string}] with the name of the basis to use. Available
#' basis are "sines", "cosines", "trigonometric", "chebyshev", "bspline". Default
#' is "bspline".
#' @param posKnots position for the knots of the B-spline. Available positions
#' are "periodic", "density", "uniform". Default is "density". Only used if
#' \code{basisName} is "bspline".
#' @param degree degree of the B-splines. Default is 3. Only used if
#' \code{basisName} is "bspline".
#' @param criterion character defining the criterion to select the best model.
#' The best model is the one with the lowest criterion value.
#' Possible values: "BIC", "AIC", "ML". Default is "BIC".
#' @param test list with test data to classify
#'
#' @examples
#' ## the famous formosat data set
#' #d <-readSplitFormosatData(path="./data/Formosat/",firstYear=2008, lastYear=2012)
#' data("multiFormosat")
#' model <- learnFuncSpectra(multiFormosat)
#' ## get summary
#' summary(model)
#'
#' ## use graphics functions
#' \dontrun{
#' plot(model)
#' plot(model$models[[1]])
#' }
#'
#' ## print model
#' \dontrun{
#' print(model)
#' }
#'
#' @return An instance of the [\code{\linkS4class{FuncModel}}] class.
#' @author Serge Iovleff
#'
learnFuncSpectra <- function( data
                            , kernelName = "gaussian", width = 50
                            , dim = 5
                            , basisName = "Bspline", posKnots = "density", degree = 3
                            , criterion = "BIC", maxValue = NULL
                            , test = NULL)
{
  # check number of class
  nbClass <- 0
  for (i in 1:length(data$labels))
  { nbClass <- max(nbClass, nlevels(factor( (data$labels)[[i]] ))) }
  if (nbClass < 2) { stop("in learnFuncSpectra, not enough class")}

  # check maxValue
  if (!is.null(maxValue))
  {
    if (!is.double(maxValue) || maxValue<= 0)
    stop("maxValue must be a positive real number\n")
  }

  # check criterion
  if (sum(criterion %in% c("BIC","AIC", "ML")) !=  1)
  { stop("criterion is not valid. See ?learnFuncSpectra for the list of valid criterion\n")}

  # check kernel name
  flag = .Call("checkKernelNames", kernelName, PACKAGE = "CloHe")
  if (!flag) { stop("in learnFuncSpectra, kernelName is incorrect")}

  # check basis name
  flag = .Call("checkBasisNames", basisName, PACKAGE = "CloHe")
  if (!flag) { stop("in learnFuncSpectra, basisName  is incorrect")}

  # additional tests for BSPLINE basis
  if (toupper(basisName) == "BSPLINE")
  {
    # check knots position
    flag = .Call("checkKnotsPositionNames", posKnots, PACKAGE = "CloHe")
    if (!flag) { stop("in learnFuncSpectra, posKnots is incorrect")}
    # check degree
    if (!is.numeric(degree)) { stop("degree must be an integer\n")}
  }
  else
  { # set arbitrary values
    degree = 0
    posKnots = "unknown"
  }

  # create vector of results
  res <- lapply( rep("FuncSpectraModel", nbClass), new)
  # launch computations
  if (length(data$spectra) == 4)
  {
    resLearn <- .Call( "launchFuncSpectra4"
                     , data
                     , list( kernelName = kernelName, width = width
                           , basisName = basisName, dim = dim, posKnots = posKnots
                           , degree = degree
                           , criterion=criterion
                           , maxValue = maxValue)
                     , res
                     , test
                     , PACKAGE = "CloHe"
                 )
  }
  else
  if (length(data$spectra) == 10)
  {
    resLearn <- .Call( "launchFuncSpectra10"
                     , data
                     , list( kernelName = kernelName, width = width
                           , basisName = basisName, dim = dim, posKnots = posKnots
                           , degree = degree
                           , criterion=criterion
                           , maxValue = maxValue)
                     , res
                     , test
                     , PACKAGE = "CloHe"
                     )
  }
  else
  {
    stop(paste("Only 4 or 10 spectra are implemented and you have ", length(data$spectra), "spectra"))
  }

  class(resLearn) <- "FuncModel"
  resLearn
}

#-----------------------------------------------------------------------
#' Definition of the [\code{\linkS4class{FuncSpectraModel}}] class
#'
#' This class defines a Gaussian model for the spectrum with mean
#' and variance varying along the time. It inherits
#' from [\code{\linkS4class{ICloHeModel}}].
#'
#' @slot mu matrix with the sampled (100 times) mean spectrum.
#' @slot sigma2 vector with the variance of each spectrum
#'
#' @seealso [\code{\linkS4class{ICloHeModel}}] class
#'
#' @examples
#' getSlots("FuncSpectraModel")
#'
#' @author Serge Iovleff
#'
#' @name FuncSpectraModel
#' @aliases FuncSpectraModel-class
#' @rdname FuncSpectraModel-class
#'
setClass(
    Class = "FuncSpectraModel",
    representation( mu           = "matrix"
                  , sigma2       = "vector"
                  , alpha        = "matrix"
                  , eigenvalues  = "vector"
                  , knots        = "vector"
                  , coefficients = "vector"
                  , tmin         = "numeric"
                  , tmax         = "numeric"
              ),
    contains = c("ICloHeModel"),
    validity = function(object)
    {
      return(TRUE)
    }
)

#' Initialize an instance of a CloHe S4 class.
#'
#' Initialization method of the [\code{\linkS4class{FuncSpectraModel}}] class.
#' Used internally in the 'CloHe' package.
#'
#' @rdname initialize-methods
#' @keywords internal
setMethod(
    f = "initialize",
    signature = c("FuncSpectraModel"),
    definition = function(.Object)
    {
      # validate
      .Object <- callNextMethod(.Object)
      validObject(.Object)
      return(.Object)
    }

)

#' @rdname print-methods
#' @aliases print FuncSpectraModel-method
#'
setMethod(
    f = "print",
    signature = c("FuncSpectraModel"),
    function(x,...){
      cat("****************************************\n")
      callNextMethod();
      cat("****************************************\n")
      cat("* mu        =\n")
      print( format(x@mu), quote = FALSE)
      cat("****************************************\n")
      cat("* sigma2    =\n")
      print( format(x@sigma2), quote = FALSE)
      cat("****************************************\n")
      cat("* alpha    =\n")
      print( format(x@alpha), quote = FALSE)
    }
)

#' @rdname show-methods
#' @aliases show-FuncSpectraModel,FuncSpectraModel,FuncSpectraModel-method
setMethod(
    f = "show",
    signature = c("FuncSpectraModel"),
    function(object)
    {
      cat("****************************************\n")
      callNextMethod();
      cat("****************************************\n")
      if (length(object@mu) != 0)
      {
        cat("* mu of spectra =\n")
        print( format(object@mu), quote = FALSE)
      }
      cat("****************************************\n")
      if (length(object@sigma2) != 0)
      {
        cat("* sigma2 of spectra =\n")
        print( format(object@sigma2), quote = FALSE)
      }
      cat("****************************************\n")
      if (length(object@sigma2) != 0)
      {
        cat("* alpha of spectra =\n")
        print( format(object@alpha), quote = FALSE)
      }
      cat("****************************************\n")
    }
)

#' @rdname summary-methods
#' @aliases summary summary,FuncSpectraModel-method
#'
setMethod(
    f = "summary",
    signature = c("FuncSpectraModel"),
    function(object, ...)
    {
      cat("**************************************************************\n")
      callNextMethod()
      cat("**************************************************************\n")
    }
)

#' Plotting of a class [\code{\linkS4class{FuncSpectraModel}}]
#'
#' Plotting data from a [\code{\linkS4class{FuncSpectraModel}}] object
#' using the estimated parameters.
#'
#' @param x an object of class [\code{\linkS4class{FuncSpectraModel}}]
#' @param ... further arguments passed to or from other methods
#'
#' @aliases plot-FuncSpectraModel
#' @docType methods
#' @rdname plot-FuncSpectraModel-method
#'
#' @seealso \code{\link{plot}}
#' @examples
#'  ## the famous formosat data set
#' \dontrun{
#'   d <-readFiles(firstYear=2012, lastYear=2012)
#'   model <- learnFuncSpectra(d)
#'   plot(r$models[[1]])
#'   }
#'
setMethod(
    f = "plot",
    signature = c(x = "FuncSpectraModel", y = "missing"),
    function(x,y,...)
    {
      nbSpectrum <- x@nbSpectrum
      nbSampling <- nrow(x@mu)
      t <- seq(from = 1, to = 366, length.out = nbSampling)
      nbrow = 1
      nbcol = nbSpectrum
      if ((nbSpectrum %% 2) == 0)
      {
        nbrow = 2
        nbcol = nbSpectrum %/% 2
      }
      if ((nbSpectrum %% 3) == 0)
      {
        nbrow = 3
        nbcol = nbSpectrum %/% 3
      }
      # get old par
      op <- par(no.readonly = TRUE) # the whole list of settable par's.
      op.palette <- colorRampPalette(c('dark blue', 'green', 'dark red'), space = "Lab")
      palette( op.palette(nbSpectrum) )
      par(mfrow = c(nbrow, nbcol))
      for (i in 1:nbSpectrum)
      {
        plot(x = t, y = x@mu[,i], type = 'l', col = i, ylab = "Values", xlab = "dates")
      }
      # restore plotting parameters
      par(op)
      palette("default")
    }
)

#' Plotting of a class [\code{FuncModel}]
#'
#' Plotting data from a [\code{FuncModel}] S3 object
#' using the estimated parameters.
#'
#' @param x an object of class [\code{FuncModel}]
#' @param y not used. Will be silently ignored
#' @param ... further arguments passed to or from other methods
#'
#' @aliases plot-FuncModel
#' @docType methods
#' @rdname plot-FuncModel-method
#'
#' @seealso \code{\link{plot}}
#' @examples
#' ## the famous formosat data set
#' \dontrun{
#'   d <-readFiles(firstYear=2012, lastYear=2012)
#'   r <- learnFuncSpectra(d)
#'   plot(r)
#'   }
#'
plot.FuncModel <- function(x,...)
{
  # get dimensions
  nbSpectrum <- x$nbSpectrum
  nbClass    <- x$nbClass
  if ((nbSpectrum %% 2) == 0)
  {
    nbRow = 2
    nbCol = nbSpectrum %/% 2
  }
  else
  {
    nbRow = 1
    nbCol = nbSpectrum
  }
  # get old par
  op <- par(no.readonly = TRUE) # the whole list of settable par's.
  fmin <- rep(1e18,nbSpectrum)
  fmax <- rep(0,nbSpectrum)
  classLabels <- c()
  for (k in 1:nbClass)
  {
    classLabels <- c(classLabels, x$models[[k]]@classLabel)
    for (i in 1:nbSpectrum)
    {
      fmi <- min(x$models[[k]]@mu[,i])
      fma <- max(x$models[[k]]@mu[,i])
      fmin[i] <- min(fmin[i], fmi)
      fmax[i] <- max(fmax[i], fma)
    }
  }
  palette(rainbow(nbClass+3))
  # cluster parameters
  par(cex = .75, oma = c(4, 1, 1, 1)) # font size and margin. Let a big bottom margin for the legend
  par(mfrow = c(nbRow, nbCol))
  # one plot per spectrum
  for (i in 1:nbSpectrum)
  {
    # plot spectrum i, class 1
    nbTimes <- nrow(x$models[[1]]@mu)
    t <- seq(from = x$tMin, to = x$tMax, length.out = nbTimes)
    plot( x = t, y = x$models[[1]]@mu[,i], type = 'l', col = 1
        , xlim = c(x$tMin, x$tMax), ylim = c(fmin[i],fmax[i])
        , ylab = "Values", xlab = "dates")
    title(main = paste("Class Means, Spectrum ", i))
    # plot spectrum i, other class
    for (k in 2:nbClass)
    {
      nbTimes <- nrow(x$models[[k]]@mu)
      t <- seq(from = x$models[[k]]@tmin, to = x$models[[k]]@tmax, length.out = nbTimes)
      lines( x = t, y = x$models[[k]]@mu[,i], type = 'l', col = k
           , ylab = "Values", xlab = "dates")
    }
  }
  #close.screen(all.screens = TRUE)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend( "bottom"
        , cex = 0.8
        , horiz = TRUE
        , inset = c(0, 0), bty = "n"
        , pch = rep(19, nbClass), col = 1:nbClass, legend = as.character(classLabels)
        )
  # restore plotting parameters
  par(op)
  palette("default")
  invisible()
}

#' summary of a class [\code{FuncModel}]
#'
#' summary data from a [\code{FuncModel}] S3 object.
#'
#' @param x an object of class [\code{FuncModel}]
#'
#' @aliases summary-FuncModel
#' @docType methods
#' @rdname summary-FuncModel-method
#'
#' @examples
#' ## the famous formosat data set
#' \dontrun{
#'   d <-readFiles(firstYear=2012, lastYear=2012)
#'   r <- learnFuncSpectra(d)
#'   summary(r)
#'   }
#'
summary.FuncModel <- function(x)
{
  cat("**************************************************************\n")
  cat("* nbSample     = ", x$nbSample,"\n")
  cat("* nbClass      = ", x$nbClass,"\n")
  cat("* nbSpectrum   = ", x$nbSpectrum,"\n")
  for( k in (1:(x$nbClass)) )
  { summary(x$models[[k]])}
  cat("**************************************************************\n")
}

#' print of a class [\code{FuncModel}]
#'
#' print data from a [\code{FuncModel}] S3 object.
#'
#' @param x an object of class [\code{FuncModel}]
#'
#' @aliases print-FuncModel
#' @docType methods
#' @rdname print-FuncModel-method
#'
#' @examples
#' ## the famous formosat data set
#' \dontrun{
#'   d <-readFiles(firstYear=2012, lastYear=2012)
#'   r <- learnFuncSpectra(d)
#'   print(r)
#'   }
#'
print.FuncModel <- function(x)
{
  cat("**************************************************************\n")
  cat("* nbSample     = ", x$nbSample,"\n")
  cat("* nbClass      = ", x$nbClass,"\n")
  cat("* nbSpectrum   = ", x$nbSpectrum,"\n")
  cat("* labels       =\n")
  print( format(x$labels), quote = FALSE)
  cat("* zi           =\n")
  print( format(x$zi), quote = FALSE)
  for( k in (1:(x$nbClass)) )
  { print(x$models[[k]])}
  cat("**************************************************************\n")
}

#' @title get parts of a class [\code{FuncModel}]
#'
#' get log-likelihoods from a [\code{FuncModel}] S3 object.
#'
#' @param x an object of class [\code{FuncModel}]
#'
#' @aliases likelihood-FuncModel
#' @docType methods
#' @rdname FuncModel-methods
#'
#' @examples
#' ## the famous formosat data set
#' \dontrun{
#'   d <-readFiles(firstYear=2012, lastYear=2012)
#'   r <- learnFuncSpectra(d)
#'   likelihood(r)
#'   }
#'
lnLikelihood <- function(x)
{
  likelihood <- (1:x$nbClass)
  for( k in (1:(x$nbClass)) )
  {likelihood[k] <- x$models[[k]]@lnLikelihood}
  likelihood
}

#' @title get parts of a class [\code{FuncModel}]
#'
#' get criteria from a [\code{FuncModel}] S3 object.
#'
#' @aliases criteria-FuncModel
#' @docType methods
#' @rdname FuncModel-methods
#'
#' @examples
#' ## the famous formosat data set
#' \dontrun{
#'   d <-readFiles(firstYear=2012, lastYear=2012)
#'   r <- learnFuncSpectra(d)
#'   criteria(r)
#'   }
#'
criteria <- function(x)
{
  criteria <- (1:x$nbClass)
  for( k in (1:(x$nbClass)) )
  {criteria[k] <- x$models[[k]]@criterion}
  criteria
}


