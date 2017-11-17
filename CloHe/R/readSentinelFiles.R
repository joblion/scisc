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
readMultiSentinelData <- function(path="./data/")
{
  library(h5)
  f <- h5file(name = paste(path, "gp.hdf5", sep = ""))
  y <- f["y"][]
  x <- f["x"][]
  t <- f["dates"][]

  labels <- list(years1 = as.vector(y))
  times  <- list(years1 = t)
  spect1 <- list(years1 = as.matrix(x[,1:33]))
  spect2 <- list(years1 = as.matrix(x[,34:66]))
  spect3 <- list(years1 = as.matrix(x[,67:99]))
  spect4 <- list(years1 = as.matrix(x[,100:132]))
  spect5 <- list(years1 = as.matrix(x[,133:165]))
  spect6 <- list(years1 = as.matrix(x[,166:198]))
  spect7 <- list(years1 = as.matrix(x[,199:231]))
  spect8 <- list(years1 = as.matrix(x[,232:264]))
  spect9 <- list(years1 = as.matrix(x[,265:297]))
  spect10 <- list(years1 = as.matrix(x[,298:330]))

  clouds <- list(years1 = matrix(0, nrow = length(y), ncol = length(t) ))

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
