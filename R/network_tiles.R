
#' Make a set of tiles that tesselate space
#'
#' @param model a model, defined as a compound [list]
#'
#' @return the model, a compound [list]
#' @export
make_tiles = function(model){
  UseMethod("make_tiles", model$Mpar)
}

#' Make a set of tiles that tesselate space for the BQ model
#'
#'
#' @param model a model, defined as a compound [list]
#'
#' @importFrom deldir deldir tile.list
#'
#' @return the model, a compound [list]
#' @export
make_tiles.BQ = function(model){with(model,{
  model$tilesB = deldir::tile.list(deldir::deldir(b))
  model$tilesQ = deldir::tile.list(deldir::deldir(q))
  return(model)
})}

#' Make a set of tiles that tesselate space for the BQS model
#'
#' @param model a model, defined as a compound [list]
#'
#' @importFrom deldir deldir tile.list
#'
#' @return the model, a compound [list]
#' @export
make_tiles.BQS = function(model){with(model,{
  model$tilesB = tile.list(deldir(b))
  model$tilesQ = tile.list(deldir(q))
  model$tilesS = tile.list(deldir(s))
  return(model)
})}
