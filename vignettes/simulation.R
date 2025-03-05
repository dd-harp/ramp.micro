## ----suppressMessages=T, echo=F-----------------------------------------------
library(viridisLite)
library(knitr)
library(viridis)
library(ramp.micro)

## -----------------------------------------------------------------------------
set.seed(24328)
bb = unif_xy(256, -17, 17) 
qq = unif_xy(289, -17, 17) 

## ----fig.height=6, fig.width=6------------------------------------------------
par(mar = c(1,1,1,1))
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
plot_Psi_BQmod(bb, qq, Psi_bb, Psi_qb, Psi_bq, Psi_qq)

## -----------------------------------------------------------------------------
model = setup_model(bb, qq, kFb = kFb, kFq = kFq)

## -----------------------------------------------------------------------------
names(model)

## -----------------------------------------------------------------------------
names(model$Mpar)

## -----------------------------------------------------------------------------
names(model$Lpar)

## -----------------------------------------------------------------------------
model <- SIM(model)

## ----fig.height=6, fig.width=6------------------------------------------------
par(mar = c(1,1,1,1))
plot_points(model, 
            bwts=model$states$M$B_t[[201]], 
            qwts= model$states$M$B_t[[201]], 
            max_pt_sz=2)

## -----------------------------------------------------------------------------
model <- steady_state(model)

## ----fig.height=6, fig.width=6------------------------------------------------
par(mar = c(1,1,1,1))
plot_points(model, 
            bwts=model$steady$M$B,
            qwts=model$steady$M$Q,
            max_pt_sz=2)

