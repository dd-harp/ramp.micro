
makeTiles = function(mod){
  UseMethod("makeTiles", mod)
}

makeTiles.BQ = function(mod){with(mod,{
  mod$tilesB = tile.list(deldir(b))
  mod$tilesQ = tile.list(deldir(q))
  return(mod)
})}

makeTiles.BQS = function(mod){with(mod,{
  mod$tilesB = tile.list(deldir(b))
  mod$tilesQ = tile.list(deldir(q))
  mod$tilesS = tile.list(deldir(s))
  return(mod)
})}
