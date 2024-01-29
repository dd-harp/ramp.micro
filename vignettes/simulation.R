## ----suppressMessages=T, echo=F-----------------------------------------------
library(viridisLite)
library(knitr)
library(viridis)
library(motrap.micro)
#devtools::load_all()

## -----------------------------------------------------------------------------
set.seed(24328)
bb = unif_xy(256, -17, 17) 
qq = unif_xy(289, -17, 17) 

## ----fig.height=7, fig.width=7------------------------------------------------
plot_points_bq(bb, qq)

## -----------------------------------------------------------------------------
kFb = make_kF_exp(k=2, s=1, gamma=1.5)
kFq = make_kF_exp(k=2, s=2, gamma=2)

## ----fig.height=4, fig.width=6------------------------------------------------
dd = seq(0, 5, by = 0.01)
plot(dd, kFb(dd), type = "l", col = "darkred", xlab = "Distance", ylab = "Weight")
lines(dd, kFq(dd), type = "l", col = "darkblue")

## -----------------------------------------------------------------------------
Psi_bb = make_Psi_xx(bb, kFb)
Psi_qb = make_Psi_xy(bb, qq, kFq)
Psi_bq = make_Psi_xy(qq, bb, kFb)
Psi_qq = make_Psi_xx(qq, kFq)

## ----fig.height=8, fig.width=8------------------------------------------------
par(mfcol = c(2,2), mar = c(1,2,1,2))
plot_Psi_BQ(bb, qq, Psi_bb, Psi_qb, Psi_bq, Psi_qq)

## ----eval=F-------------------------------------------------------------------
#  #devtools::load_all()

## -----------------------------------------------------------------------------
opts_d = list(kFb = kFb, kFq = kFq)
model = setup_model(bb, qq, dispersal_opts=opts_d)

## -----------------------------------------------------------------------------
names(model)

## -----------------------------------------------------------------------------
names(model$Mpar)

## -----------------------------------------------------------------------------
names(model$Lpar)

