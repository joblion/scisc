#
# Author: iovleff
###############################################################################

#' Read the csv formosat files by concatenating all the years in a single series.
#'
#' Each file is in the format i_year.csv, with i=1  for blue spectrum, i=2 for green spectrum
#' i=3 for red spectrum and i=4 for infra-red spectrum.
#' @param path a string with the path to the files
#' @param firstYear first year of observations
#' @param lastYear last year of observations
#' @return a list with the labels of each observations, the days of each observation (starting from
#' day 1 of firstYear), the four spectrums (b,g,r,ir), the clouds presence.
#' @rdname  readFormosatFiles.rd
readFormosatData <- function(path="./data/Formosat/", firstYear=2008, lastYear=2014)
{
  if (lastYear < firstYear) {stop("lastYear must be greater or equal to lastYear")}
  # create years
  years <- firstYear:lastYear
  nameYears <- as.character(years)
  names(years) <- nameYears
  # blue values
  bList <- vector(mode = "list", length = length(nameYears))
  names(bList) <- nameYears
  # green values
  gList <- vector(mode = "list", length = length(nameYears))
  names(gList) <- nameYears
  # red values
  rList <- vector(mode = "list", length = length(nameYears))
  names(rList) <- nameYears
  # infrared values
  iList <- vector(mode = "list", length = length(nameYears))
  names(iList) <- nameYears
  # date values
  dList <- vector(mode="list", length = length(nameYears))
  names(dList) <- nameYears
  # read files
  days <- 0; times <- NULL; xb <- NULL; xg <- NULL; xr <- NULL; xi <- NULL; mb <- NULL;
  for(y in nameYears)
  {
    bList[[y]] <- read.csv( file=paste(path, "1_", y, ".csv", sep=""), header = FALSE, sep = ",")
    gList[[y]] <- read.csv( file=paste(path, "2_", y, ".csv", sep=""), header = FALSE, sep = ",")
    rList[[y]] <- read.csv( file=paste(path, "3_", y, ".csv", sep=""), header = FALSE, sep = ",")
    iList[[y]] <- read.csv( file=paste(path, "4_", y, ".csv", sep=""), header = FALSE, sep = ",")
    # dates
    dList[[y]] <- days+read.csv( file = paste(path, "sample_time_", y, ".csv", sep=""), header = FALSE, sep = ",")
    nbSample <- length(dList[[y]][,1]);
    days = days + 365
    # handle bissextile years
    if (years[y] == 4 * ( years[y] %/% 4 )) { days = days+1 }
    # create clouds
    mb <- cbind( mb, as.matrix(bList[[y]][(2):(1+nbSample)]))
    # create spectrum
    xb <- cbind( xb, as.matrix(bList[[y]][(2+nbSample):(1+2*nbSample)]))
    xg <- cbind( xg, as.matrix(gList[[y]][(2+nbSample):(1+2*nbSample)]))
    xr <- cbind( xr, as.matrix(rList[[y]][(2+nbSample):(1+2*nbSample)]))
    xi <- cbind( xi, as.matrix(iList[[y]][(2+nbSample):(1+2*nbSample)]))
    # create times
    times <- rbind(times, dList[[y]])
  }
  # create labels
  labels <- as.integer(bList[[nameYears[1]]][,1])
  # return data in a list
  list( labels  = labels
      , times   = as.vector(as.matrix(times))
      , xb      = as.matrix(xb)
      , xg      = as.matrix(xg)
      , xr      = as.matrix(xr)
      , xi      = as.matrix(xi)
      , clouds  = as.matrix(mb)
      )
}

#' Read the csv formosat files by adding each year as a new series.
#'
#' Each file is in the format i_year.csv, with i=1  for blue spectrum, i=2 for green spectrum
#' i=3 for red spectrum and i=4 for infra-red spectrum.
#' @param path a string with the path to the files
#' @param firstYear first year of observations
#' @param lastYear last year of observations
#' @return a list with the labels of each observations, the days of each observation (starting from
#' day 1 of firstYear), the four spectrums (b,g,r,ir), the clouds presence.
#'
#' @rdname  readFormosatFiles.rd
readMultiFormosatData <- function(path="./data/Formosat/", firstYear=2008, lastYear=2014)
{
  if (lastYear<firstYear) {stop("lastYear must be greater or equal to lastYear")}
  # create years
  years <- firstYear:lastYear
  nameYears <- as.character(years)
  names(years) <- nameYears
  # blue values
  bList <- vector(mode = "list", length = length(nameYears))
  names(bList) <- nameYears
  xb <- vector(mode = "list", length = length(nameYears))
  names(xb) <- nameYears
  # green values
  gList <- vector(mode = "list", length = length(nameYears))
  names(gList) <- nameYears
  xg <- vector(mode = "list", length = length(nameYears))
  names(xg) <- nameYears
  # red values
  rList <- vector(mode = "list", length = length(nameYears))
  names(rList) <- nameYears
  xr <- vector(mode = "list", length = length(nameYears))
  names(xr) <- nameYears
  # infrared values
  iList <- vector(mode = "list", length = length(nameYears))
  names(iList) <- nameYears
  xi <- vector(mode = "list", length = length(nameYears))
  names(xi) <- nameYears
  # date values
  dList <- vector(mode = "list", length = length(nameYears))
  names(dList) <- nameYears
  times <- vector(mode = "list", length = length(nameYears))
  names(times) <- nameYears
  # clouds values
  mb <- vector(mode = "list", length = length(nameYears))
  names(mb) <- nameYears
  # labels
  labels <- vector(mode = "list", length = length(nameYears))
  names(labels) <- nameYears
  # read files
  for (i in 1:length(years))
  {
    y <- years[i]
    # spectra
    bList[[y]] <- read.csv( file=paste(path, "1_", y, ".csv", sep = ""), header = FALSE, sep = ",")
    gList[[y]] <- read.csv( file=paste(path, "2_", y, ".csv", sep = ""), header = FALSE, sep = ",")
    rList[[y]] <- read.csv( file=paste(path, "3_", y, ".csv", sep = ""), header = FALSE, sep = ",")
    iList[[y]] <- read.csv( file=paste(path, "4_", y, ".csv", sep = ""), header = FALSE, sep = ",")
    # dates
    dList[[y]] <- read.csv( file = paste(path, "sample_time_", y, ".csv", sep=""), header = FALSE, sep = ",")

    nbSample <- length(dList[[y]][,1]);
    # create clouds
    mb[[i]] <- as.matrix(bList[[y]][(2):(1+nbSample)])
    # create spectrum
    xb[[i]] <- as.matrix(bList[[y]][(2+nbSample):(1+2*nbSample)])
    xg[[i]] <- as.matrix(gList[[y]][(2+nbSample):(1+2*nbSample)])
    xr[[i]] <- as.matrix(rList[[y]][(2+nbSample):(1+2*nbSample)])
    xi[[i]] <- as.matrix(iList[[y]][(2+nbSample):(1+2*nbSample)])
    # create times
    times[[i]] <- dList[[y]][,1] - 1
    # create labels
    labels[[i]] <- as.integer(bList[[y]][,1]) # could be any spectrum
  }
  # return data in a list
  list( labels  = labels
      , times   = times
      , spectra = list(xb = xb
                      ,xg = xg
                      ,xr = xr
                      ,xi = xi)
      , clouds  = mb
      )
}
