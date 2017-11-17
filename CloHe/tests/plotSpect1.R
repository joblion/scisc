spect1 <- multiSentinel$spectra$spect1$years1
clouds <- multiSentinel$clouds$years1
times  <- multiSentinel$times$years1

nbSample <- dim(clouds)[1]
nbTimes <- dim(clouds)[2]

y <- spect1[1, clouds[1,]==0]
x <- times[ clouds[1,]==0]
plot(x, y, type = 'l', col = "blue", xlim=c(0,355), ylim=c(0,6000))

for (i in 2:nbSample)
{
  y <- spect1[i, clouds[i,]==0]
  x <- times[ clouds[i,]==0]
  lines(x, y, col = "blue")
}
