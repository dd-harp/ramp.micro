
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
  clusters = net$clusters_walktrap
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
  for( i in 1:dm){
    ix = which(memix==i)
    cx = mean(xy[ix,1])
    cy = mean(xy[ix,2])
    cxy = cbind(cx,cy)
    cntr = rbind(cntr, cxy)
  }
  self = diag(pop)
  p.self = self/colSums(pop)
  pop1 = pop
  diag(pop1)<-0
  p.connect = pop1%*% diag(1/colSums(pop))
  model$patches = list(pop=pop, self=self, p.self=p.self, p.connect=p.connect, connect=pop1, centers=cntr)
  return(model)
}

#' Visualize the patchespopulation approximation
#'
#' @param model a ramp.micro model object
#' @param i the graph
#' @param cut a cutat argument, the # of communities
#' @param lwd a plotting argument
#' @param bbend the bend in the arrows
#' @param mtl a title for the graph
#'
#' @returns invisible(NULL)
#' @export
plot_patches = function(model, i, cut=NULL, lwd=2, bbend=3, mtl = NULL){
  net = get_net(model, i)
  model = net2patches(model, i, cut)
  n = dim(model$patches$centers)[1]
  clrs = viridis::turbo(n)[sample(1:n)]
  with(model, frame_bq(b,q, mtl))
  graphics::points(model$patches$centers, col = clrs)
  net <- make_convex_hulls(model, net, cut, f_color = viridis::turbo)
  plot_convex_hulls(net, lwd=lwd)
  with(model, add_bent_arrows_xx(patches$centers, patches$p.connect, bbend=bbend, clr=clrs))
  return(invisible())
}
