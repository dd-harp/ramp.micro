
KxyProfiles = function(
    mod,
    lng=0.002, cutat=0.99, pw=0.5,
    adj1=1.5, adj2=2.5,
    clrLA="#4361eeCC", clrLB="#fe5f55CC",
    clrA='red', clrB='blue',
    clrNB='darkred', clrNQ='darkblue'
){
  with(mod,{
    par (mfrow = c(1,2), mar = c(2,2,2,2))
    frame_bq(b, q, mtl = expression(K[q%->%b]))
    if(exists("s")) addP.s(s, adj=0.5)
    arrowsX2Y(q, b, Kqb, clr=clrLA, lng=lng, cutat=cutat)
    addP.q(q, adj=adj1, clr=clrA)
    addP.b(b, rowSums(Kqb), pw=pw, adj=adj2, clr=clrNB)

    frame_bq(b, q, mtl = expression(K[b%->%q]))
    if(exists("s")) addP.s(s, adj=0.5)
    arrowsX2Y(b, q, Kbq, clr=clrLB, lng=lng, cutat=cutat)
    addP.b(b, adj=adj1, clr=clrB)
    addP.q(q, rowSums(Kbq), pw=pw, adj=adj2, clr=clrNQ)
  }
  )}













KxxProfiles = function(
    mod, lng=0.02, wd=.1, withHist=FALSE,
    cutat=0.99,
    adj1=0.5, adj2=2,
    clrLA="#e2739655", clrLB="#abc4ff55", clrLS='#00000022',
    clrA="#fe5f5599", clrB="#858ae399",
    clrNB="#cc444bCC", clrNQ="#4e148cCC"
){
  with(mod,
       {
         if(withHist == TRUE) par(mfcol = c(2,2), mar = c(2,2,2,2))
         if(withHist == FALSE) par(mfcol = c(1,2), mar = c(2,2,2,2))

         frame_bq(b, q, mtl = expression(K[b %->%b]))
         if(exists("s")) addP.s(s, adj=adj1)
         arrowsX2X(b, Kbb, lng=lng, wd=wd, cutat=cutat, clrA=clrLA, clrS=clrLS)
         addP.q(q, adj=adj1, clr=clrB)
         addP.bb(b, Kbb, adj=adj2, colA=clrNB, colB=clrNB)

         if(withHist == TRUE) hist(rowSums(Kbb), col="darkred")

         frame_bq(b, q, mtl = expression(K[q %->%q]))
         if(exists("s")) addP.s(s, adj=adj1)
         arrowsX2X(q, Kqq, lng=lng, wd=wd, cutat=cutat, clrA=clrLB, clrS=clrLS)
         addP.b(b, adj=adj1, clr=clrA)
         addP.qq(q, Kqq, adj=adj2, colA=clrNQ, colB=clrNQ)

         if(withHist == TRUE) hist(rowSums(Kqq), col="darkblue")
       }
  )}






makeKGV = function(modelObject){
  modelObject = makeKbq(modelObject)
  modelObject = makeKqb(modelObject)
  Kbb = with(modelObject, Kqb %*% Kbq)
  modelObject$Kbb = Kbb
  Kqq = with(modelObject, Kbq %*% Kqb)
  modelObject$Kqq = Kqq
  modelObject = computeG(modelObject)
  modelObject = computeGG(modelObject)
  modelObject = computeV(modelObject)
  modelObject = computeVC(modelObject)
  return(modelObject)
}



