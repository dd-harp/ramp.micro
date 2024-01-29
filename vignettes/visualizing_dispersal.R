## ----suppressMessages=T-------------------------------------------------------
library(viridisLite)
library(knitr)
library(viridis)
library(motrap.micro)
set.seed(25)

## ----echo=F, eval=F-----------------------------------------------------------
#  devtools::load_all()

## -----------------------------------------------------------------------------
nPv = 37 
nPw = 73 
dd = 10
Pv = unif_xy(nPv,-dd, dd) 
Pw = unif_xy(nPw,-dd, dd) 

## -----------------------------------------------------------------------------
M_wv = make_Psi_xy(Pv, Pw, make_kF_exp())
M_vw = make_Psi_xy(Pw, Pv, make_kF_exp())

## ----fig.height=4, fig.width=8, echo=F----------------------------------------
par(mar = c(1,1,1,1), mfcol = c(1,2))
plot_matrix_xy(Pv, Pw, M_wv, r=0.25) -> wts
plot_matrix_xy(Pw, Pv, M_vw, r=0.25, arw_clr = "cornflowerblue") -> wts

## ----fig.height=3, fig.width=5------------------------------------------------
hist(rowSums(M_wv), xlab = "Mosquito Density at Destination", main = expression(sum(M[list(j,i)], i)))

## ----fig.height=3, fig.width=5------------------------------------------------
hist(log10(M_wv), 40, xaxt = "n")
axis(1, -5:-1, c(expression(10^-5), expression(10^-4), "0.001", "0.01", "0.1"))
segments(-2, 0, -2, 50, col = "red", lwd=2)

## ----fig.height=4, fig.width=8, echo=F----------------------------------------
par(mar = c(1,1,1,1), mfcol = c(1,2))
plot_matrix_xy(Pv, Pw, M_wv, min_edge_frac = 0) 
plot_matrix_xy(Pv, Pw, M_wv, min_edge_frac = .01) 


## ----fig.height=6, fig.width=5------------------------------------------------
par(mfrow = c(2,1), mar = c(4,2,1,1))
ot = order(as.vector(M_wv), decreasing=T)
Mo = as.vector(M_wv)[ot]
plot(Mo, type = "l", ylab = "PMF", xlab = "")
ix1 = max(which(Mo>0.01))
points(ix1, Mo[ix1], pch = 15, col = "darkblue")

ecdf = cumsum(Mo)/sum(M_wv)
plot(ecdf, type = "l", ylim = c(0,1), xlab = "", ylab = "eCMF")
ix2 = min(which(ecdf >0.9))
points(ix2, ecdf[ix2], pch=19, col = "darkred")

## -----------------------------------------------------------------------------
Mvv = M_vw %*% M_wv
Mww = M_wv %*% M_vw

## ----fig.height=4, fig.width=8, echo=F----------------------------------------
par(mar = c(1,1,1,1), mfcol = c(1,2))
plot_matrix_xx(Pv, Mvv/mean(Mvv), lwd=5, r=0.5) -> wts
plot_matrix_xx(Pw, Mww/mean(Mww), lwd=5, r=0.45, pt_clr = "lightblue", arw_clr = "cornflowerblue") -> wts

## -----------------------------------------------------------------------------
nEdges = function(nS, nD, frac=0.01){
  dd = 10
  S = unif_xy(nS,-dd, dd) 
  D = unif_xy(nD,-dd, dd) 
  M = make_Psi_xy(S, D, make_kF_exp())
  ot = order(as.vector(M), decreasing=T)
  Mo = as.vector(M)[ot]
  ix1 = max(which(Mo>frac))
  c(ix1, nS*nD, ix1/nS/nD)
}

## -----------------------------------------------------------------------------
t(replicate(20, nEdges(27, 43))) -> out
c(nEdges = mean(out[,1]))
c(fracEdges =  mean(out[,3]))

## -----------------------------------------------------------------------------
t(replicate(20, nEdges(270, 430))) -> out
c(nEdges = mean(out[,1]))
c(fracEdges =  mean(out[,3]))

## ----fig.height=10, fig.width=10, eval=F--------------------------------------
#  nPv = 270
#  nPw = 430
#  dd=10
#  Pv = unif_xy(nPv,-dd, dd)
#  Pw = unif_xy(nPw,-dd, dd)
#  M = make_Psi_xy(Pv, Pw, make_kF_exp())
#  plot_matrix_xy(Pv, Pw, M, min_edge_frac = 0.02)

## -----------------------------------------------------------------------------
t(replicate(20, nEdges(270*5, 430*5))) -> out
c(nEdges = mean(out[,1]))
c(fracEdges =  mean(out[,3]))

## ----fig.height=10, fig.width=10, eval=F--------------------------------------
#  nPv = 270*5
#  nPw = 430*5
#  dd=10
#  Pv = unif_xy(nPv,-dd, dd)
#  Pw = unif_xy(nPw,-dd, dd)
#  M = make_Psi_xy(Pv, Pw, make_kF_exp())
#  plot_matrix_xy(Pv, Pw, M, min_edge_frac = 0.01)

