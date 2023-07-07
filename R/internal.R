segment <- function(A, B, n) {
  t_ <- seq(0, 1, length.out = n)
  t(vapply(t_, function(t) {A + t*(B-A)}, numeric(2L)))
}

ray <- function(O, A, n, tmax, OtoA) {
  if(OtoA) {
    t_ <- seq(0, tmax, length.out = n)
  } else {
    t_ <- seq(-tmax, 0, length.out = n)
  }
  t(vapply(t_, function(t) {O + t*(A-O)}, numeric(2L)))
}
