
PsiProfile_BQ = function(b,q,
                         Psi_bb, Psi_qb,
                         Psi_bq, Psi_qq){
  par (mfrow = c(2,2), mar = c(2,2,2,2))
  ## b -> b
  frame_bq(b,q, mtl=expression(Psi[bb]))
  arrowsX2X(b, Psi_bb)
  addP.bb(b, Psi_bb, adj=3)

  ## q -> b
  frame_bq(b,q, mtl = expression(Psi[qb]))
  arrowsX2Y(q, b, Psi_qb)
  addP.q(q, adj=1)
  addP.b(b, rowSums(Psi_qb), adj=3)

  ## b -> q
  frame_bq(b,q, mtl = expression(Psi[bq]))
  arrowsX2Y(b, q, Psi_bq, clr="blue")
  addP.q(q, rowSums(Psi_bq), adj=3)
  addP.b(b, adj=1)

  ##q -> q
  frame_bq(b,q, mtl = expression(Psi[qq]))
  arrowsX2X(q,Psi_qq, clrA = "blue")
  addP.qq(q, Psi_qq, adj=3)
}

PsiProfile = function(simObject){
  UseMethod("PsiProfile", simObject)
}

PsiProfile.BQ = function(simObject){with(simObject,{
  PsiProfile_BQ(b,q,Psi_bb, Psi_qb, Psi_bq, Psi_qq)
})}

PsiProfile_BQS = function(b,q,s,
                          Psi_bb, Psi_qb, Psi_sb,
                          Psi_bq, Psi_qq, Psi_sq,
                          Psi_bs, Psi_qs, Psi_ss){
  par (mfrow = c(3,3), mar = c(2,2,2,2))
  ## b -> b
  frame_bqs(b,q,s, mtl = expression(Psi[b%->%b]))
  symX2X(b, Psi_bb)
  flowX2X(b, Psi_bb)
  addP.bb(b, Psi_bb, adj=3)

  ## q -> b
  frame_bqs(b,q,s, mtl = expression(Psi[q%->%b]))
  arrowsX2Y(q, b, Psi_qb)
  addP.q(q, adj=1)
  addP.b(b, rowSums(Psi_qb), adj=3)

  ## s -> b
  frame_bqs(b,q,s, mtl = expression(Psi[s%->%b]))
  arrowsX2Y(s, b, Psi_sb, clr="red")
  addP.s(s, adj=1)
  addP.b(b, rowSums(Psi_sb), adj=3)

  ## b -> q
  frame_bqs(b,q,s, mtl = expression(Psi[b%->%q]))
  arrowsX2Y(b, q, Psi_bq, clr="blue")
  addP.q(q, rowSums(Psi_bq), adj=3)
  addP.b(b, adj=1)

  ##q -> q
  frame_bqs(b,q,s, mtl = expression(Psi[q%->%q]))
  symX2X(q,Psi_qq)
  flowX2X(q, Psi_qq, clr = "blue")
  addP.qq(q, Psi_qq, adj=3)

  ##s -> q
  frame_bqs(b,q,s, mtl = expression(Psi[s%->%q]))
  arrowsX2Y(s, q, Psi_sq, clr="blue")
  addP.q(q, rowSums(Psi_sq), adj=3)
  addP.s(s, adj=1)

  ## b -> s
  frame_bqs(b,q,s, mtl = expression(Psi[b%->%s]))
  arrowsX2Y(b, s, Psi_bs, clr="orange")
  addP.b(b, adj=1)
  addP.s(s, rowSums(Psi_bs), adj=3)

  ##q -> s
  plotPoints_bqs(b,q,s, mtl = expression(Psi[q%->%s]))
  arrowsX2Y(q, s, Psi_qs, clr="orange")
  addP.q(q, adj=1)
  addP.s(s, rowSums(Psi_qs), adj=3)

  ##s -> s
  frame_bqs(b,q,s, mtl = expression(Psi[s%->%s]))
  symX2X(s, Psi_ss)
  flowX2X(s, Psi_ss, clr="orange")
  addP.ss(s, Psi_ss, adj=3)
}

PsiProfile.BQS = function(simObject){with(simObject,{
  PsiProfile_BQS(b,q,s, Psi_bb, Psi_qb, Psi_sb,
                 Psi_bq, Psi_qq, Psi_sq,
                 Psi_bs, Psi_qs, Psi_ss)
})}
