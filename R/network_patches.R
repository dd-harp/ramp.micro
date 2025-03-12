
#' Use communities to define patches
#'
#' @param model a ramp.micro model object
#' @param i which graph (see get_net)
#' @param cut cutat argument, # communities (optional)
#' @importFrom igraph membership cut_at
#' @returns a ramp.micro model object
#' @export
net2patches =function(model, i, cut=NULL){
  net = get_net(model, i)
  clusters = net$walktrap_clusters
  if(net$type == "b") xy = model$b
  if(net$type == "q") xy = model$q
  M = get_matrix(model, i)
  if(is.null(cut)){
    memix = membership(clusters)
  } else{
    memix = cut_at(clusters, cut)
  }
  dm = max(memix)
  pop = matrix(0, dm, dm)
  for(i in 1:dm){
    for(j in 1:dm){
      ixi = which(memix==i)
      ixj = which(memix==j)
      pop[i,j] = sum(M[ixi,ixj])
    }
  }
  cntr = c()
  chull = list()
  bq = list()
  for(i in 1:dm){
    ix = which(memix==i)
    cx = mean(xy[ix,1])
    cy = mean(xy[ix,2])
    cxy = cbind(cx,cy)
    chull[[i]] = make_convex_hull_xy(xy[ix,])
    bq[[i]] = cxy
    cntr = rbind(cntr, cxy)
  }
  self = diag(pop)
  p.self = self/colSums(pop)
  pop1 = pop
  diag(pop1)<-0
  p.connect = pop1%*% diag(1/colSums(pop))
  if(net$type == "b")
   model$b_patches = list(i=1:dm, pop=pop, chull=chull, cxy=bq, self=self, p.self=p.self, p.connect=p.connect, connect=pop1, centers=cntr)
  if(net$type == "q")
   model$q_patches = list(i=1:dm, pop=pop, chull=chull, cxy=bq, self=self, p.self=p.self, p.connect=p.connect, connect=pop1, centers=cntr)
  return(model)
}

#' Visualize the patchespopulation approximation
#'
#' @param model a ramp.micro model object
#' @param i the graph
#' @param cut a cutat argument, the # of communities
#' @param f_color a function that returns a list of colors (e.g. viridis::turbo)
#' @param stretch make the hull slightly larger or slightly smaller
#' @param lwd a plotting argument
#' @param bbend the bend in the arrows
#' @param mtl a title for the graph
#'
#' @returns a ramp.micro model object
#' @export
plot_patches = function(model, i, cut=NULL, f_color = viridis::turbo, stretch=0.1, lwd=2, bbend=3, mtl = NULL){
  net = get_net(model, i)
  model = net2patches(model, i, cut)
  if(net$type == "b") patches <- model$b_patches
  if(net$type == "q") patches <- model$q_patches
  n = dim(patches$centers)[1]
  clrs = viridis::turbo(n)
  if(net$type == "b") model$b_patches$clrs = clrs
  if(net$type == "q") model$q_patches$clrs = clrs
  with(model, frame_bq(b,q, mtl))
  graphics::points(patches$centers, pch=21, bg = "white", col=clrs, cex=2)
  for(j in 1:n){
    sxy <- patches$chull[[j]]
    sxy = stretch_convex_hull(sxy, 1+stretch)
    graphics::polygon(sxy[,1], sxy[,2], border=clrs[j], lwd=lwd)
  }
  add_bent_arrows_xx(patches$centers, patches$p.connect, bbend=bbend, endd=0.65, clr=clrs)
  return(model)
}

#' Plot the ouptuts of a graph using
#'
#' @param model a `ramp.micro` model object
#' @param graphs a graphswork object
#' @param cut optional arguent for cut_at
#' @param alg walktrap = "wt" or greedy = "gr"
#' @param f_color a function that returns a list of colors (e.g. viridis::turbo)
#' @param min_edge_frac the fraction of the mass to plot
#' @param cx cex for plotting points
#' @param pw power relationship for scaling point size: pw=1 is linear
#' @param mtl a plot title
#' @param stretch make the hull slightly larger or slightly smaller
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return invisible(NULL)
#' @export
plot_meta = function(model, graph, cut=NULL, alg="wt",
                      f_color = viridis::turbo,
                      min_edge_frac = 0.01, cx=2,
                      pw=1, mtl = "",
                      stretch=0.1,
                      r=.02, arw_lng=0.05, lwd=2){with(model,{
  if(alg == "wt") clusters = graph$walktrap_clusters
  if(alg == "gr") clusters = graph$greedy_clusters
  memix = if(is.null(cut)){
    membership(clusters)
  } else {
    cut_at(clusters, cut)
  }
  nC = max(memix)
  clrs = f_color(nC)
  frame_bq(b,q,mtl)
  for(i in 1:nC){
    mx <- max(rowSums(graph$M))
    ix = which(memix == i)
    M = graph$M[ix,][,ix]
    sx <- max(rowSums(M))
    if(graph$type == "b"){
      add_arrows_xx(b[ix,], M, arw_clr = clrs[i], lwd=lwd*sx/mx, min_edge_frac=min_edge_frac)
      add_points_b(b[ix,], wts_b=rowSums(M), clr=clrs[i], cx_b=cx*sx/mx, pw=pw)
      sxy = make_convex_hull_xy(b[ix,])
      sxy = stretch_convex_hull(sxy, 1+stretch)
      graphics::polygon(sxy[,1], sxy[,2], border=clrs[i], lwd=lwd)
    }
    if(graph$type == "q"){
      add_arrows_xx(q[ix,], M, arw_clr = clrs[i], lwd=lwd*sx/mx, min_edge_frac=min_edge_frac)
      add_points_q(q[ix,], wts_q=rowSums(M), clr=clrs[i], cx_q=cx*sx/mx, pw=pw)
      sxy = make_convex_hull_xy(q[ix,])
      sxy = stretch_convex_hull(sxy, 1+stretch)
      graphics::polygon(sxy[,1], sxy[,2], border=clrs[i], lwd=lwd)
    }
  }
  return(invisible())
                      })}
