
frame_bq = function(b,q,mtl=NULL){
  plot(rbind(b,q), type = "n", xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", frame.plot=T,
       main = mtl)
}

frame_bqs = function(b,q,s, mtl=NULL){
  plot(rbind(b,q,s), type = "n", xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", frame.plot=T,
       main = mtl)
}

addP.b = function(b, wts=1, pw=1, adj=2, clr="#fe5f55CC"){
  points(b, col = clr, pch = 15, cex=adj*wts^pw/max(wts^pw))
}

addP.q = function(q, wts=1, pw=1, adj=2, clr="#4361eeCC"){
  points(q, col=clr, pch=19, cex=adj*wts^pw/max(wts^pw))
}

addP.s = function(s, wts=1, pw=1, adj=2, clr="#ff4800CC"){
  points(s, col=clr, pch=17, cex=adj*wts^pw/max(wts^pw))
}

addP.bb = function(b, M, pw=1, adj=2, colA="#cc444b66", colB="#cc444bCC"){
  addP.b(b, as.vector(rowSums(M)), pw, adj, colA)
  diag(M) <- 0
  addP.b(b, as.vector(rowSums(M)), pw, adj, colB)
}

addP.qq = function(q, M, pw=1, adj=2, colA="#4e148c66", colB="#4e148cCC"){
  addP.q(q, as.vector(rowSums(M)), pw, adj, colA)
  diag(M) <- 0
  addP.q(q, as.vector(rowSums(M)), pw, adj, colB)
}

addP.ss = function(s, M, pw=1, adj=2, colA="#8bc74266", colB="#8bc742CC"){
  addP.s(s, as.vector(rowSums(M)), pw, adj, colA)
  diag(M) <- 0
  addP.s(s, as.vector(rowSums(M)), pw, adj, colB)
}

plotPoints = function(model, ...){
  UseMethod("plotPoints", model)
}

plotPoints.BQ = function(model, ...){with(model,{
  plotPoints_bq(b, q, ...)
})}

plotPoints.BQS = function(model, ...){with(model,{
  plotPoints_bqs(b, q, s, ...)
})}

plotPoints_bq = function(b, q,
                         bwts=1, qwts=1, pw=1, adj=2,
                         mtl=NULL){
  frame_bq(b,q,mtl)
  addP.b(b, bwts, pw, adj)
  addP.q(q, qwts, pw, adj)
}

plotPoints_bqs = function(b, q, s,
                          bwts=1, qwts=1, swts=1, pw=1, adj=2,
                          mtl=NULL){

  frame_bqs(b,q,s,mtl)
  addP.s(s, swts, pw, adj)
  addP.b(b, bwts, pw, adj)
  addP.q(q, qwts, pw, adj)
}


