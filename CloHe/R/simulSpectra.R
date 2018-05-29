#-----------------------------------------------------------------------
#' Create a list with a simulated data set of spectra
#'
#' @param nbPixel number of pixels belonging to class k
#' @param nbCluster number of cluster
#' @param nbSpectrum number of spectra
#' @param nbSampling number of sampling
#' @param mean a function returning the mean value of the spectra. mean(t) returns
#' a matrix of size nbSpectrum x nbSampling.
#' @param sigma a vector of size nbSpectrum giving the variance level of
#' the spectrum
#' @param kernelName [\code{string}] with the kernel to use for the covariance matrix.
#' Available kernels are "gaussian", "exponential", "rationalQuadratic".
#' Default is "gaussian".
#' @param width the width of the kernel to use. Default is 50.
#'
#' @examples
#'
#'
#' @return A list with the spectra
#' @author Serge Iovleff & Asmita Poddar
#'

library(mvtnorm)
mean=function(t, nbSpectrum, nbCluster)
{
  #matrix(rep(cos(t*(2*pi)/365), nbSpectrum), nrow=nbSpectrum, byrow=T)}
  res <- array(0, c(nbCluster, nbSpectrum, length(t)));
  a0 = 2

  #for(i in 1:nbSpectrum)
  #{
  #  if (i%%2==0)
  #   res[i,]=cos(t*i*(2*pi)/365)
  #  else
  #    res[i,]=sin(t*i*(2*pi)/365)
  #}

  for(i in 1:nbCluster)
  {
    for(j in 1:nbSpectrum)
    {
      s = rep(0, length(t))
      for(k in 1:5)
      {
        a = 1/sqrt(2*pi)*exp(-k*k)
        lambda = 1/sqrt(2*pi)*exp(-(k-1)*(k-1))
        s = s+a0*a*cos((2*pi*lambda*t*j*i)/365)
      }
      res[i,j,]=s
    }
  }
  res
}

simulSpectra<-function( nbPixel, nbCluster, nbSpectrum, nbSampling = 20
                      , sigma = dnorm(rnorm(nbSpectrum)) , kernelName = "gaussian"
                      , width = 50 )
{
  times <- seq(from=0, to=365, length.out=nbSampling)
  means <- mean(times, nbSpectrum, nbCluster)
  data <- array(0, dim = c(nbPixel, nbSpectrum, nbSampling))
  process <- matrix(0, nrow=nbSpectrum, ncol= nbSampling)

  #creating a vector of size nbPixel containing the labels (number of labels = nbCluster)
  #the probablilty of each cluster being between 0 and 1
  labels <- sample(1:nbCluster, nbPixel , prob = dnorm(rnorm(nbCluster)) , replace = T)
  ##prob = rep(1, nbCluster)

   for (i in 1:nbPixel)
   {
     k <- labels[i]
     for ( s in 1:nbSpectrum)
     {
       covariance <- gaussianKernelCov(times, width, sigma[s])
       process[s,] <- rmvnorm(1, mean = means[k,s,], sigma = covariance )
     }
     data[i,,] <- process
   }
  list(labels=labels , times = times
       , means = means, sigma = sigma, process = process);
  #labels, times, spectra, clouds


}

gaussianKernelCov <- function(times, width, sigma)
{
  tLength <- length(times)
  res <- matrix(0, nrow=tLength, ncol=tLength)
  for(i in 1:tLength)
  {
    for(j in 1:tLength)
    {
      s <- times[i]
      t <- times[j]
      res[i,j] <- exp(-abs(s-t)^2/width)
    }
  }
  sigma * res
}
