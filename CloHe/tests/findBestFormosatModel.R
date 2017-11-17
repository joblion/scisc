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

findBestFormosatModel <- function( widths = seq(from = 5, to = 150, by = 5)
                                 , dims = 6:40
                                 , degree = 3
                                 )
{
  library("CloHe")
  multiFormosat <- readMultiFormosatData()
  widths <- widths
  dims  <- dims
  degree <- degree

  nbLabels <- 0;
  for (y in 1:length(multiFormosat$labels))
  { nbLabels = nbLabels + length(multiFormosat$labels[[y]])}

  ## Launch Gaussian tests
  predictGaussian <- matrix(0,nrow=length(dims), ncol= length(widths))
  for (i in 1:length(dims))
  {
    for (j in 1:length(widths))
    {
      rtest <- learnFuncSpectra(multiFormosat, kernelName = "gaussian", width = widths[j], posKnots = "density", degree = degree, dim = dims[i])
      predictGaussian[i,j] <- sum(diag(buildConfusionMatrix(rtest$labels, rtest$zi)))/nbLabels
    }
  }
  resultsGaussian <- list(widths = widths, dims = dims, predicts = predictGaussian)

  ## Launch Laplace
  predictLaplace <- matrix(0,nrow = length(dims), ncol = length(widths))
  for (i in 1:length(dims))
  {
    for (j in 1:length(widths))
    {
      rtest <- learnFuncSpectra(multiFormosat, kernelName = "laplace", width = widths[j], posKnots = "density", degree = degree, dim = dims[i])
      predictLaplace[i,j] <- sum(diag(buildConfusionMatrix(rtest$labels, rtest$zi)))/nbLabels
    }
  }
  resultsLaplace <- list(widths = widths, dims = dims, predicts = predictLaplace)

  ## Launch rationalQuadratic
  predictRational <- matrix(0,nrow=length(dims), ncol= length(widths))
  for (i in 1:length(dims))
  {
    for (j in 1:length(widths))
    {
      rtest <- learnFuncSpectra(multiFormosat, kernelName = "rationalQuadratic", width = widths[j], posKnots = "density", degree = degree, dim = dims[i])
      predictRational[i,j] <- sum(diag(buildConfusionMatrix(rtest$labels, rtest$zi)))/nbLabels
    }
  }
  resultsRational <- list(widths = widths, dims= dims, predicts = predictRational)
  list(resultsGaussian, resultsLaplace, resultsRational)
}

  # plot results
  # plot(widths, resultsGaussian6$predicts, type='l', col = "blue", xlab = "widths", ylab = "predicted", ylim = c(0.65,0.77))
  # lines(widths, resultsGaussian7$predicts, type='l', col = "green")
  # lines(widths, resultsGaussian8$predicts, type='l', col = "orange")
  # lines(widths, resultsGaussian9$predicts, type='l', col = "red")
  #
  # lines(widths, resultsLaplace6$predicts, type='l', col = "blue", lty = 2)
  # lines(widths, resultsLaplace7$predicts, type='l', col = "green", lty = 2)
  # lines(widths, resultsLaplace8$predicts, type='l', col = "orange", lty = 2)
  # lines(widths, resultsLaplace9$predicts, type='l', col = "red", lty = 2)
  #
  # lines(widths, resultsRational6$predicts, type='l', col = "blue", lty = 4)
  # lines(widths, resultsRational7$predicts, type='l', col = "green", lty = 4)
  # lines(widths, resultsRational8$predicts, type='l', col = "orange", lty = 4)
  # lines(widths, resultsRational9$predicts, type='l', col = "red", lty = 4)
  #
  # legend("bottomleft", legend=c("6 Control Points", "7 Control Points", "8 Control Points", "9 Control Points"),
  #        col=c("blue", "green","orange","red"), lty=1, cex=0.8)
  # legend("bottomright", legend=c("Gaussian Kernel", "Laplace Kernel", "Rational quadra."),
  #         lty=c(1,2,4), cex=0.8)

  # compute and plot best model
  # rtest <- learnFuncSpectra(multiFormosat, kernelName = "gaussian", width = 60, posKnots = "density", degree = 2, dim = 8)
  # mat <- buildConfusionMatrix(rtest$labels, rtest$zi)
  # bestPred <- sum(diag(mat))/nbLabels
  # res <- rep(0, 7)
  # nbDiff <- rep(0, 1029)
  # for (j in 1:1029)
  # {
  #   # count
  #   for (i in 0:6) { res[i + 1] = rtest$zi[1029 * i + j] }
  #   nbDiff[j] <- nlevels(factor(res))
  # }
  # table(factor(nbDiff))
  # mat
  # bestPred
  # # barplot
  # plot(factor(nbDiff))
