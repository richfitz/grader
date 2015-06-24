## TODO: Option to do this with forward differences, where we'll reuse
## the first element.
##
## TODO: Option to use f(x) wherever possible.
##
## TODO: function requirements are quite strict; the matrix thing is
## going to be annoying.  Might be possible to relax that for a list
## of vectors.

##' Compute gradient of a function
##' @title Compute gradient of a function
##' @param f Function
##' @param x Point to compute gradient at
##' @param eps See \code{numDeriv::grad}
##' @param d see \code{numDeriv::grad}
##' @param r see \code{numDeriv::grad}
##' @param zero_tol see \code{numDeriv::grad}
##' @export
gradient <- function(f, x, eps=1e-4, d=0.0001, r=4,
                     zero_tol=sqrt(.Machine$double.eps/7e-7)) {
  pts <- gradient_points(x, eps, d, r, zero_tol)
  y <- gradient_eval(f, pts)
  gradient_extrapolate(y, pts)
}

gradient_points <- function(x,
                            eps=1e-4, d=0.0001, r=4,
                            zero_tol=sqrt(.Machine$double.eps/7e-7)) {
  v <- 2
  ## This all corresponds to numDeriv's 'case 2'
  n <- length(x) #  number of parameters (theta)

  ## Initial offset:
  h <- abs(d * x) + eps * (abs(x) < zero_tol)
  pts <- vector("list", r * n)
  dim(pts) <- c(r, n)
  dx <- matrix(NA, r, n)

  j <- seq_len(n)
  for (k in seq_len(r)) {
    for (i in seq_len(n)) {
      dx_i <- h * (i == j)
      pts[[k, i]] <- rbind(x + dx_i, x - dx_i)
    }
    dx[k, ] <- 2 * h
    h <- h / v
  }

  ## Instead, let's save this as a 3d array?
  ret <- do.call("rbind", pts)
  attr(ret, "dim_y") <- c(2L, r * n)
  attr(ret, "dx") <- dx
  attr(ret, "n") <- n
  attr(ret, "r") <- r
  class(ret) <- "gradient_points"
  ret
}

gradient_eval <- function(f, pts) {
  y <- f(pts)
  dim(y) <- attr(pts, "dim_y")
  y
}

gradient_extrapolate <- function(y, pts) {
  dx <- attr(pts, "dx")
  r  <- attr(pts, "r")
  n  <- attr(pts, "n")
  a <- (y[1,] - y[2,]) / dx # series of estimates; central fd.
  for (m in seq_len(r - 1L)) {
    four_m <- 4.0 ^ m
    a_next <- matrix(NA, r - m, n)
    for (i in seq_len(nrow(a_next))) {
      a_next[i,] <- (a[i + 1L,] * four_m - a[i,]) / (four_m - 1.0)
    }
    a <- a_next
  }
  drop(a)
}
