
predictGaussianMutSigmat <- function(xb, xg, xr, xi, clouds, labels, models)
{
  labels <- as.vector(labels)
  # dimensions
  nbSample <- length(labels)
  nbTimes  <- ncol(xb)
  nbClass  <- length(models)
  # predicted labels
  predict <- integer(nbSample);

  for (i in 1:nbSample)
  {
    smax <- -Inf
    kmax <- 1L
    for (k in 1:nbClass)
    {
      sum = 0;
      for (t in 1:nbTimes)
      {
        x <- c(xb[i,t],xg[i,t],xr[i,t],xi[i,t]);
        if (clouds[i,t] == 0)
        { sum <- sum + dmvnorm(x, mean = models[[k]]@mut[[t]], sigma = models[[k]]@sigmat[[t]], log = TRUE) }
      }
      if (smax < sum)
      {
        kmax <- k
        smax <- sum
      }
    }
    predict[i] <- kmax;
  }
  buildConfusionMatrix(labels, predict)
}

buildConfusionMatrix <- function(trueLabels, predictLabels)
{
  trueLabels <- as.vector(trueLabels)
  predictLabels <- as.vector(predictLabels)
  if (length(trueLabels) != length(predictLabels))
  { stop("trueLabels and predictLabels are not of the same length")}

  classNumbers <- as.integer(levels(factor(trueLabels)))
  nbClass <- length(classNumbers);
  nbSample <- length(trueLabels);

  confusionMatrix <- matrix(0, nrow = nbClass, ncol = nbClass);
  rownames(confusionMatrix) <- paste0(rep("T", length.out = nbClass), classNumbers)
  colnames(confusionMatrix) <- paste0(rep("P", length.out = nbClass), classNumbers)

  for (i in 1:nbSample)
  {
    confusionMatrix[trueLabels[i], predictLabels[i]] =   confusionMatrix[trueLabels[i], predictLabels[i]] + 1
  }
  confusionMatrix
}
