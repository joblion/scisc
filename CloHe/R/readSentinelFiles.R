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
readMultiSentinelSmoothedData <- function(path="./data/", size=NULL)
{
  library(h5)
  f <- h5file(name = paste(path, "gp.hdf5", sep = ""))
  y <- f["y"][]
  x <- f["x"][]
  t <- f["dates"][]
  nbSample <- length(as.vector(y))
  if (is.null(size))
  { samples <- 1:nbSample }
  else
  { samples <- sample(1:nbSample, size = size)}

  labels <- list(years1 = as.vector(y)[samples])
  times  <- list(years1 = t)
  spect1 <- list(years1 = as.matrix(x[samples,1:33]))
  spect2 <- list(years1 = as.matrix(x[samples,34:66]))
  spect3 <- list(years1 = as.matrix(x[samples,67:99]))
  spect4 <- list(years1 = as.matrix(x[samples,100:132]))
  spect5 <- list(years1 = as.matrix(x[samples,133:165]))
  spect6 <- list(years1 = as.matrix(x[samples,166:198]))
  spect7 <- list(years1 = as.matrix(x[samples,199:231]))
  spect8 <- list(years1 = as.matrix(x[samples,232:264]))
  spect9 <- list(years1 = as.matrix(x[samples,265:297]))
  spect10 <- list(years1 = as.matrix(x[samples,298:330]))
  clouds <- list(years1 = matrix(0, nrow = length(samples), ncol = length(t) ))

  # return data in a list
  list( labels  = labels
      , times   = times
      , spectra = list( spect1      = (spect1)
                      , spect2      = (spect2)
                      , spect3      = (spect3)
                      , spect4      = (spect4)
                      , spect5      = (spect5)
                      , spect6      = (spect6)
                      , spect7      = (spect7)
                      , spect8      = (spect8)
                      , spect9      = (spect9)
                      , spect10     = (spect10)
                      )
      , clouds = clouds
      )
}

# raws files
readMultiSentinelRawData <- function(path="./data/", size=NULL)
{
  library(h5)

  f <- h5file(name = paste(path, "raw.hdf5", sep = ""))
  times <- f["dates"][]
  labels <- f["y"][]
  spectra <- f["x"][]

  f <- h5file(name = paste(path, "mask.hdf5", sep = ""))
  clouds <- f["mask"][]

  nbSample <- length(as.vector(labels))
  if (is.null(size))
  {
    size <- nbSample
    samples <- 1:nbSample
  }
  else
  { samples <- sample(1:nbSample, size = size)}

  times  <- list(years1 = times)
  labels <- list(years1 = as.vector(labels)[samples])
  clouds <- list(years1 = (as.matrix(clouds)[samples,]))
  spectra<- as.matrix(spectra[samples,])

  spect1 <- list(years1 = spectra[,1:30])
  for (i in 1:size) { for (j in 1:30) { if(spect1$years1[i,j] == -10000) { clouds$years1[i,j] <- 1} } }

  spect2 <- list(years1 = spectra[,31:60])
  for (i in 1:size) { for (j in 1:30) { if(spect2$years1[i,j] == -10000) { clouds$years1[i,j] <- 1} } }

  spect3 <- list(years1 = spectra[,61:90])
  for (i in 1:size) { for (j in 1:30) { if(spect3$years1[i,j] == -10000) { clouds$years1[i,j] <- 1} } }

  spect4 <- list(years1 = spectra[,91:120])
  for (i in 1:size) { for (j in 1:30) { if(spect4$years1[i,j] == -10000) { clouds$years1[i,j] <- 1} } }

  spect5 <- list(years1 = spectra[,121:150])
  for (i in 1:size) { for (j in 1:30) { if(spect5$years1[i,j] == -10000) { clouds$years1[i,j] <- 1} } }

  spect6 <- list(years1 = spectra[,151:180])
  for (i in 1:size) { for (j in 1:30) { if(spect6$years1[i,j] == -10000) { clouds$years1[i,j] <- 1} } }

  spect7 <- list(years1 = spectra[,181:210])
  for (i in 1:size) { for (j in 1:30) { if(spect7$years1[i,j] == -10000) { clouds$years1[i,j] <- 1} } }

  spect8 <- list(years1 = spectra[,211:240])
  for (i in 1:size) { for (j in 1:30) { if(spect8$years1[i,j] == -10000) { clouds$years1[i,j] <- 1} } }

  spect9 <- list(years1 = spectra[,241:270])
  for (i in 1:size) { for (j in 1:30) { if(spect9$years1[i,j] == -10000) { clouds$years1[i,j] <- 1} } }

  spect10 <- list(years1 = spectra[,271:300])
  for (i in 1:size) { for (j in 1:30) { if(spect10$years1[i,j] == -10000) { clouds$years1[i,j] <- 1} } }

  # return data in a list
  list( labels  = labels
      , times   = times
      , spectra = list( spect1      = (spect1)
                      , spect2      = (spect2)
                      , spect3      = (spect3)
                      , spect4      = (spect4)
                      , spect5      = (spect5)
                      , spect6      = (spect6)
                      , spect7      = (spect7)
                      , spect8      = (spect8)
                      , spect9      = (spect9)
                      , spect10     = (spect10)
                      )
      , clouds = clouds
      )
}
