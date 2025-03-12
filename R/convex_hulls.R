# Visualize communities as convex hulls

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

# Visualize communities as convex hulls

#' Make a set of convex hulls to visualize the community
#'
#' @param model a model defined as a compound [list]
#' @param graph a graphwork object
#' @param cut if !null, the number of communities for igraph::cut_at
#' @param clrs a set of colors
#' @param f_color a function that returns a list of colors (e.g. viridis::turbo)
#'
#' @return the graphwork object with convex hulls attached
#' @export
make_convex_hulls = function(model, graph, cut=NULL, clrs=NULL, f_color = viridis::turbo){
  clusters = graph$walktrap_clusters
  memix = if(is.null(cut)){
    igraph::membership(clusters)
  } else {
    igraph::cut_at(clusters, cut)
  }
  if(graph$type == "b"){
    xy = model$b
  }
  if(graph$type == "q"){
    xy = model$q
  }
  if(graph$type == "bq"){
    xy = with(model,rbind(b, q))
  }
  graph$convex_hulls = list()
  if(is.null(clrs)) clrs = f_color(max(memix))
  for(i in 1:max(memix)){
    graph$convex_hulls[[i]] = list()
    graph$convex_hulls[[i]]$xy = make_convex_hull_i(i, memix, xy)
    graph$convex_hulls[[i]]$clr = clrs[i]
  }
  return(graph)
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
  ix = which(memix == i)
  return(make_convex_hull(xy[ix,]))
}

#' Make the convex hull for the i^th community
#'
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
#' @param graph a graphwork object
#' @param stretch make the hull slightly larger or slightly smaller
#' @param lwd hull line width
#'
#' @return invisible(NULL)
#' @export
plot_convex_hulls = function(graph, stretch=1.1, lwd=2){
  n = length(graph$convex_hulls)
  for(i in 1:n){
    with(graph$convex_hulls[[i]],{
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

