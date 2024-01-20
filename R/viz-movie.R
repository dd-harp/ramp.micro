
makeMovie = function(mod, moviestemx){
  UseMethod("makeMovie", mod)
}

makeMovie.BQ = function(mod, moviestem="BQmovie"){
  Bt = mod$states$Bt; Qt=mod$states$Qt
  Tmax=dim(Bt)[2]
  sclB = max(sqrt(Bt))
  sclQ = max(sqrt(Qt))
  framenames = c()
  with(mod,{
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
    av_encode_video(framenames, moviename, verbose=F) -> out
  })
}



makeMovie.BQS = function(mod,  moviestem="BQSmovie"){
  Bt = mod$states$Bt; Qt=mod$states$Qt; St=mod$states$St
  Tmax=dim(Bt)[2]
  sclB = max(sqrt(Bt))
  sclQ = max(sqrt(Qt))
  sclS = max(sqrt(St))
  framenames = c()
  with(mod,{
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
    av_encode_video(framenames, moviename, verbose=F) -> out
  })
}
