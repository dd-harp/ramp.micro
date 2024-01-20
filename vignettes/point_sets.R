## ----suppressMessages=T-------------------------------------------------------
library(viridisLite)
library(knitr)
library(viridis)
library(motrap.micro)

## ----fig.height=5, fig.width=5------------------------------------------------
par(bty = "o")
xy = runif(50, -5, 5)
plot(xy, xaxt = "n", yaxt = "n", xlab = "", ylab="")

