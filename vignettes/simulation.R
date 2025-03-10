## ----suppressMessages=T, echo=F-----------------------------------------------
library(viridisLite)
library(knitr)
library(viridis)
library(ramp.micro)

## ----echo=F-------------------------------------------------------------------
#devtools::load_all()

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
Psi_bb = make_Psi_xx(bb, kFb, stay=0.5)
Psi_qb = make_Psi_xy(bb, qq, kFq)
Psi_bq = make_Psi_xy(qq, bb, kFb)
Psi_qq = make_Psi_xx(qq, kFq, stay=0.5)

## ----fig.height=8, fig.width=8------------------------------------------------
par(mfcol = c(2,2), mar = c(1,2,1,2))
plot_all_Psi_BQ(bb, qq, Psi_bb, Psi_qb, Psi_bq, Psi_qq, r=0.5, cx_D=1.5, cx_S=0.5)

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
            wts_b = model$states$M$B_t[[60]], 
            wts_q = model$states$M$B_t[[60]], 
            cx_b=2, cx_q=2)

## -----------------------------------------------------------------------------
model <- steady_state(model)

## ----fig.height=6, fig.width=6------------------------------------------------
par(mar = c(1,1,1,1))
plot_points(model, 
            wts_b=model$steady$M$B,
            wts_q=model$steady$M$Q,
            cx_b=2, cx_q=2)

