## -----------------------------------------------------------------------------
library(ramp.micro) 

## -----------------------------------------------------------------------------
set.seed(24328)
bb = unif_xy(256, -17, 17) 
qq = unif_xy(289, -17, 17) 
ker_b = make_kF_exp(k=2, s=1, gamma=1.5)
ker_q = make_kF_exp(k=2, s=2, gamma=2)
dispO = list(kFb=ker_b, kFq =ker_q)
bq_mod1 = setup_model(b=bb, q=qq, dispersal_opts = dispO) 
bq_mod1 <- basic_analysis(bq_mod1)

## ----fig.height=4, fig.width=12-----------------------------------------------
par(mar = c(1,1,4,4), mfrow=c(1,3))
plot_subgraph(bq_mod1, bq_mod1$graphs$GG_graph, stretch=0.1, min_edge_frac=0.1, cut=15, mtl = "Egg Dispersal")
with(bq_mod1, frame_bq(b,q))
add_hulls(bq_mod1, bq_mod1$graphs$GG_graph, stretch=-0.01, min_edge_frac=0.1, cut=15)
add_hulls(bq_mod1, bq_mod1$graphs$VC_graph, stretch=-0.01, min_edge_frac=0.1, cut=15,)
plot_subgraph(bq_mod1, bq_mod1$graphs$VC_graph, stretch=0.1, min_edge_frac=0.1, cut=15, mtl = "Potential Transmission")

