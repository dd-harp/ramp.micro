# Visualize communities as convex hulls

#' Make a set of convex hulls to visualize the community
#'
#' @param model a model defined as a compound [list]
#' @param net a network object
#' @param cut if !null, the number of communities for igraph::cut_at
#' @param clrs a set of colors
#' @param f_color a function that returns a list of colors (e.g. viridis::turbo)
#'
#' @return the network object with convex hulls attached
#' @export
make_convex_hulls = function(model, net, cut=NULL, clrs=NULL, f_color = viridis::turbo){
  clusters = net$walktrap_clusters
  memix = if(is.null(cut)){
    igraph::membership(clusters)
  } else {
    igraph::cut_at(clusters, cut)
  }
  if(net$type == "b"){
    xy = model$b
  }
  if(net$type == "q"){
    xy = model$q
  }
  if(net$type == "bq"){
    xy = with(model,rbind(b, q))
  }
  net$convex_hulls = list()
  if(is.null(clrs)) clrs = f_color(max(memix))
  for(i in 1:max(memix)){
    net$convex_hulls[[i]] = list()
    net$convex_hulls[[i]]$xy = make_convex_hull_i(i, memix, xy)
    net$convex_hulls[[i]]$clr = clrs[i]
  }
  return(net)
}

#' Make the convex hull for the i^th community
#'
#' @param i the i^th community
#' @param memix the community membership indices
#' @param xy the point set
#'
#' @return the xy points that define a convex hull
#'
#' @return a [list] with a hull
#' @export
make_convex_hull_i = function(i, memix, xy){
  ixj = which(memix == i)
  hpts <- grDevices::chull(xy[ixj,])
  ixk = c(hpts, hpts[1])
  hxy = xy[ixj[ixk],]
  return(hxy)
}

#' Make the convex hull for the i^th community
#'
#' @param i the i^th community
#' @param memix the community membership indices
#' @param xy the point set
#'
#' @return the xy points that define a convex hull
#'
#' @return a [list] with a hull
#' @export
make_convex_hull_xy = function(xy){
  hpts <- grDevices::chull(xy)
  ixk = c(hpts, hpts[1])
  hxy = xy[ixk,]
  return(hxy)
}


#' Add the convex hulls to a framed plot
#'
#' @param net a network object
#' @param stretch make the hull slightly larger or slightly smaller
#' @param lwd hull line width
#'
#' @return invisible(NULL)
#' @export
plot_convex_hulls = function(net, stretch=1.1, lwd=2){
  n = length(net$convex_hulls)
  for(i in 1:n){
    with(net$convex_hulls[[i]],{
      sxy = stretch_convex_hull(xy, stretch)
      graphics::polygon(sxy[,1], sxy[,2], border=clr, lwd=lwd)
    })
  }
  return(invisible())
}

#' Add convex hulls around the points in a community
#'
#' @param memix the membership index
#' @param xy the xy locasions
#' @param clrs colors
#' @param stretch a stretch factor for the hulls
#' @param lwd the line width for plotting
#' @param llty the line type
#'
#' @returns invisible(NULL)
#' @export
add_convex_hulls = function(memix, xy, clrs, stretch=0.1, lwd=2, llty=1){
  for(i in 1:max(memix)){
    hxy = make_convex_hull_i(i, memix, xy)
    sxy = stretch_convex_hull(hxy, 1+stretch)
    graphics::polygon(sxy[,1], sxy[,2], border=clrs[i], lwd=lwd, llty=llty)
  }
  return(invisible())
}


#' Stretch (or shrink) the convex hull for plotting
#'
#' @param xy a convex hull
#' @param fac a factor to stretch (>1) or shrink (<1) the hull for plotting
#'
#' @return a new point set
#' @export
stretch_convex_hull = function(xy, fac){
  cx = mean(xy[,1])
  cy = mean(xy[,2])
  xn = xy[,1] - cx
  yn = xy[,2] - cy

  r = sqrt(xn^2 + yn^2)
  theta = atan2(yn,xn)
  sxy =  cbind(
    r*fac*cos(theta) + cx,
    r*fac*sin(theta) + cy
  )
  return(sxy)
}


#' Make a set of convex hulls to visualize the community
#'
#' @param model a model defined as a compound [list]
#' @param net a network object
#' @param cut if !null, the number of communities for igraph::cut_at
#' @param f_color a function that returns a list of colors (e.g. viridis::turbo)
#'
#' @return the network object with convex hulls attached
#' @export
mod_convex_hulls.i = function(model, i, cut=NULL, f_color=turbo, stretch=0.1, lwd=2){
  net = getNet.i(model, i)
  mod_convex_hulls.i(model, net, cut, clr_scheme, stretch, lwd)
}

#' Make a set of convex hulls to visualize the community
#'
#' @param model a model defined as a compound [list]
#' @param net a network object
#' @param cut if !null, the number of communities for igraph::cut_at
#' @param clrs a list of colors
#' @param f_color a function that returns a list of colors (e.g. viridis::turbo)
#' @param stretch a stretch factor for the hulls
#' @param lwd the line width for plotting
#'
#' @return the network object with convex hulls attached
#' @export
mod_convex_hulls = function(model, net, cut=NULL, clrs=NULL, f_color=viridis::turbo, stretch=0.1, lwd=2){
  clusters = net$walktrap_clusters
  memix = if(is.null(cut)){
    membership(clusters)
  } else {
    cut_at(clusters, cut)
  }
  if(net$type == "b"){
    xy = model$b
  }
  if(net$type == "q"){
    xy = model$q
  }
  if(net$type == "bq"){
    xy = with(model,rbind(b, q))
  }
  if(is.null(clrs)) clrs = f_color(max(memix))
  add_convex_hulls(memix, xy, clrs, stretch, lwd)
}
