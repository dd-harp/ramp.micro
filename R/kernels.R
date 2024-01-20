
kFmake.exp = function(k=1, s=2, gamma=1,  w=1){
  return(function(dd, w=1){
    wij = w*(exp(-k*(dd/s)^gamma))
    wij/max(wij)
  })
}

kFmake.pwr = function(delta=1, s=1, w=1){
  return(function(dd, w=1){
    wij = w/(dd+s)^delta
    wij/max(wij)
  })
}

kFmake.mix = function(p=0.001, k=1, s=2, gamma=1, delta=1, w=1){
  return(function(dd, w=1){
    wij = (1-p)*w*(exp(-k*(dd/s)^gamma)) + p*w/(dd+s)^delta
    wij/max(wij)
  })
}

makePsi = function(S, D, kF, w=1){
  lS = length(S[,1])
  lD = length(D[,1])
  K = matrix(0, lD, lS)
  if(length(w)==1) w=rep(w, lD)
  for(i in 1:lS){
    K[,i] = kF(sqrt((S[i,1]-D[,1])^2 + (S[i,2]-D[,2])^2), w)
    K[,i] = K[,i]/sum(K[,i])
  }
  K
}

makePsi_stay = function(S, kF, w=1, stay=0){
  lS = dim(S)[1]
  K = matrix(0, lS, lS)
  if(length(w)==1) w=rep(w, lS)
  for(i in 1:lS){
    K[,i] = kF(sqrt((S[i,1]-S[,1])^2 + (S[i,2]-S[,2])^2), w)
    K[i,i] = 0
    K[,i] = (1-stay)*K[,i] /sum(K[,i])
    K[i,i] = stay
  }
  return(K)
}
