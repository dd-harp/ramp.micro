

#' Setup a square lattice model with one point set offset and nested within the other
#'
#' @param N the size of the big lattice
#' @param kFb a blood feeding dispersal kernel
#' @param kFq an egg laying dispersal kernel
#' @param q_outside a [logical] switch: if true, habitats define the big lattice
#'
#' @return a model as a compound [list]
#' @export
make_model_squareLattice=function(N, kFb, kFq, q_outside=TRUE){
  x0 = (N-1)/2
  x1 = (N-2)/2
  if(q_outside == TRUE){
    q0 = lattice(N, -x0, x0)
    b0 = lattice(N-1, -x1, x1)
  } else {
    b0 = lattice(N, -x0, x0)
    q0 = lattice(N-1, -x1, x1)
  }
  lattice_mod = makeSimObj_BQ(b0, q0, kFb, kFq)
  lattice_mod = steadyState(lattice_mod)
  lattice_mod = steadyState(lattice_mod)
  lattice_mod$tot = with(lattice_mod$steadyState,c(sum(B), sum(Q)))
  lattice_mod$means = with(lattice_mod$steadyState,c(mean(B), mean(Q)))
  lattice_mod$cv = with(lattice_mod$steadyState,c(cv(B), cv(Q)))
  return(lattice_mod)
}

#' Setup a model on a random uniform set of points
#'
#' @param N the size of the big lattice
#' @param kFb a blood feeding dispersal kernel
#' @param kFq an egg laying dispersal kernel
#' @param q_outside a [logical] switch: if true, habitats define the big lattice
#'
#' @return a model as a compound [list]
#' @export
make_model_unif=function(N, kFb, kFq, q_outside=TRUE){
  x0 = (N-1)/2
  x1 = (N-2)/2
  if(q_outside == TRUE){
    b0 = cbind(runif((N-1)^2, -x1, x1), runif((N-1)^2, -x1, x1))
    q0 = cbind(runif(N^2, -x0, x0), runif(N^2, -x0, x0))
  } else {
    b0 = cbind(runif(N^2, -x0, x0), runif(N^2, -x0, x0))
    q0 = cbind(runif((N-1)^2, -x1, x1), runif((N-1)^2, -x1, x1))
  }
  unif_mod = makeSimObj_BQ(b0, q0, kFb, kFq)
  unif_mod = steadyState(unif_mod)
  unif_mod = steadyState(unif_mod)
  unif_mod$tot = with(unif_mod$steadyState,c(sum(B), sum(Q)))
  unif_mod$means = with(unif_mod$steadyState,c(mean(B), mean(Q)))
  unif_mod$cv = with(unif_mod$steadyState,c(cv(B), cv(Q)))
  return(unif_mod)
}
