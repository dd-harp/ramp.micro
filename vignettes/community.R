## -----------------------------------------------------------------------------
library(ramp.micro) 

## -----------------------------------------------------------------------------
load("bq_mod1.rda")

## -----------------------------------------------------------------------------
devtools::load_all()

## ----fig.height=4, fig.width=12-----------------------------------------------

par(mar = c(1,1,4,4), mfrow=c(1,3))
plot_subgraph(bq_mod1, bq_mod1$graphs$GG_net, stretch=0.1, min_edge_frac=0.1, cut=15, mtl = "Egg Dispersal")
with(bq_mod1, frame_bq(b,q))
add_hulls(bq_mod1, bq_mod1$graphs$GG_net, stretch=-0.01, min_edge_frac=0.1, cut=15)
add_hulls(bq_mod1, bq_mod1$graphs$VC_net, stretch=-0.01, min_edge_frac=0.1, cut=15,)
plot_subgraph(bq_mod1, bq_mod1$graphs$VC_net, stretch=0.1, min_edge_frac=0.1, cut=15, mtl = "Potential Transmission")

