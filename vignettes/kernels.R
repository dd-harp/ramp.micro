## -----------------------------------------------------------------------------
library(ramp.micro)

## -----------------------------------------------------------------------------
set.seed(24328)
bb = unif_xy(96, -17, 17) 
qq = unif_xy(89, -17, 17) 
ss = unif_xy(72, -17, 17) 

## -----------------------------------------------------------------------------
kFb = make_kF_exp(k=2, s=1, gamma=1.5)
kFq = make_kF_exp(k=2, s=2, gamma=2)
kFs = make_kF_exp(k=3, s=2, gamma=2)

## ----fig.height=4, fig.width=6------------------------------------------------
dd = seq(0, 5, by = 0.01)
plot(dd, kFb(dd), type = "l", col = "darkred", xlab = "Distance", ylab = "Weight")
lines(dd, kFq(dd), type = "l", col = "darkblue")
lines(dd, kFs(dd), type = "l", col = "olivedrab4")

## ----fig.height=7, fig.width=7------------------------------------------------
par(mar = c(1,1,2,2))
Psi_bq = make_Psi_xy(qq, bb, kFb, w=1)
plot_Psi_bq(bb, qq, Psi_bq)

## ----fig.height=7, fig.width=7------------------------------------------------
par(mar = c(1,1,2,2))
Psi_bb = make_Psi_xx(bb, kFb, w=1, stay=0.5)
plot_Psi_bb(bb, qq, Psi_bb)

## ----fig.height=7, fig.width=7------------------------------------------------
par(mar = c(1,1,2,2))
Psi_bs = make_Psi_xy(ss, bb, kFb, w=1)
plot_Psi_bs(bb, qq, ss, Psi_bs)

## ----fig.height=7, fig.width=7------------------------------------------------
par(mar = c(1,1,2,2))
Psi_qb = make_Psi_xy(bb, qq, kFb, w=1)
plot_Psi_qb(bb, qq, Psi_qb)

## ----fig.height=7, fig.width=7------------------------------------------------
par(mar = c(1,1,2,2))
Psi_qq = make_Psi_xx(qq, kFq, w=1, stay=0.3)
plot_Psi_qq(bb, qq, Psi_qq)

## ----fig.height=7, fig.width=7------------------------------------------------
par(mar = c(1,1,2,2))
Psi_sb = make_Psi_xy(bb, ss, kFs, w=1)
plot_Psi_sb(bb, qq, ss, Psi_sb)

## ----fig.height=7, fig.width=7------------------------------------------------
par(mar = c(1,1,2,2))
Psi_sq = make_Psi_xy(qq, ss, kFs, w=1)
plot_Psi_sq(bb, qq, ss, Psi_sq)

## ----fig.height=7, fig.width=7------------------------------------------------
par(mar = c(1,1,2,2))
Psi_ss = make_Psi_xx(ss, kFs, w=1, stay=0.05)
plot_Psi_ss(bb, qq, ss, Psi_ss)

