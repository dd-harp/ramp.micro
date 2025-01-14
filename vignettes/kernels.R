## -----------------------------------------------------------------------------
library(ramp.micro)

## -----------------------------------------------------------------------------
set.seed(24328)
bb = unif_xy(256, -17, 17) 
qq = unif_xy(289, -17, 17) 

## -----------------------------------------------------------------------------
kFb = make_kF_exp(k=2, s=1, gamma=1.5)
kFq = make_kF_exp(k=2, s=2, gamma=2)

## ----fig.height=4, fig.width=6------------------------------------------------
dd = seq(0, 5, by = 0.01)
plot(dd, kFb(dd), type = "l", col = "darkred", xlab = "Distance", ylab = "Weight")
lines(dd, kFq(dd), type = "l", col = "darkblue")

