## ----suppressMessages=T-------------------------------------------------------
library(viridisLite)
library(knitr)
library(viridis)
library(ramp.micro)
#suppressWarnings(devtools::load_all() )
set.seed(25)

## -----------------------------------------------------------------------------
nPv = 15 
nPw = 7
dd = 10
Pv = unif_xy(nPv,-dd, dd) 
Pw = unif_xy(nPw,-dd, dd) 

## -----------------------------------------------------------------------------
M_wv = make_Psi_xy(Pv, Pw)

## ----fig.height=5, fig.width=5, echo=F----------------------------------------
par(mar = c(1,1,1,1))
plot_matrix_xy(Pv, Pw, M_wv, r=0.25) -> wts

## -----------------------------------------------------------------------------
M_vw = make_Psi_xy(Pw, Pv)
Mvv = M_vw %*% M_wv
Mww = M_wv %*% M_vw

## ----fig.height=5, fig.width=10, echo=F---------------------------------------
par(mar = c(1,1,1,1), mfcol = c(1,2))
plot_matrix_xx(Pv, Mvv, min_edge_frac=1e-3, lwd=23, r=0.25) -> wts
plot_matrix_xx(Pw, Mww, min_edge_frac=1e-3, lwd=23, r=0.25) -> wts

