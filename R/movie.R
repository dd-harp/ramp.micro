
#' Make a movie
#'
#' @param model a model, defined as a [list]
#' @param moviestemx the stem name
#'
#' @return invisible(NULL)
#' @export
make_movie = function(model, moviestemx){
  UseMethod("make_movie", model$Mpar)
}

#' Make a movie of a simulation with a BQ model
#'
#' @param model a model, defined as a [list]
#' @param moviestemx the stem name
#'
#' @return invisible(NULL)
#' @export
make_movie.BQ = function(model, moviestemx){
  Bt = model$states$Bt; Qt=model$states$Qt
  Tmax=dim(Bt)[2]
  sclB = max(sqrt(Bt))
  sclQ = max(sqrt(Qt))
  framenames = c()
  with(model,{
    for(i in 1:Tmax){
      fname = paste("./tmp/frame", 100+i, ".png", sep="")
      framenames = c(framenames, fname)
      png(fname, 480, 480)
      frame_bq(b, q, mtl = paste("t =", i))
      points(b, pch=15, col= "red", cex=2*sqrt(Bt[,i])/sclB)
      points(q, pch=19, col= "blue", cex=2*sqrt(Qt[,i])/sclQ)
      dev.off(dev.cur())
    }
    moviename = paste(moviestem, ".mp4", sep = "")
    av::av_encode_video(framenames, moviename, verbose=F) -> out
  })
  return(invisible())
}

#' Make a movie of a simulation with a BQS model
#'
#' @param model a model, defined as a [list]
#' @param moviestemx the stem name
#'
#' @return invisible(NULL)
#' @export
make_movie.BQS = function(model,  moviestemx="BQSmovie"){
  Bt = model$states$Bt; Qt=model$states$Qt; St=model$states$St
  Tmax=dim(Bt)[2]
  sclB = max(sqrt(Bt))
  sclQ = max(sqrt(Qt))
  sclS = max(sqrt(St))
  framenames = c()
  with(model,{
    for(i in 1:Tmax){
      fname = paste("./tmp/frame", 100+i, ".png", sep="")
      framenames = c(framenames, fname)
      png(fname, 480, 480)
      frame_bqs(b, q, s, mtl = paste("t =", i))
      points(b, pch=15, col= "red", cex=2*sqrt(Bt[,i])/sclB)
      points(q, pch=19, col= "blue", cex=2*sqrt(Qt[,i])/sclQ)
      points(s, pch=19, col= "orange", cex=2*sqrt(St[,i])/sclS)
      dev.off(dev.cur())
    }
    moviename = paste(moviestem, ".mp4", sep = "")
    av::av_encode_video(framenames, moviename, verbose=F) -> out
  })
  return(invisible())
}
