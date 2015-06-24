##' Compute Hessian of a function
##' @title Compute Hessian of a function
##' @param f Function
##' @param x Point to compute Hessian at
##' @param ... Additional arguments passed to \code{\link{bates_watts_D}}.
##' @export
Hessian <- function(f, x, ...) {
  hessian_D(bates_watts_D(f, x, ...))
}

gradient_D <- function(obj) {
  obj$D[seq_along(obj$x)]
}

hessian_D <- function(obj) {
  D <- obj$D
  n <- obj$p
  H <- matrix(NA, n, n)
  H[upper.tri(H, TRUE)] <- D[-seq_len(n)]
  H[lower.tri(H)] <- t(H)[lower.tri(H)]
  H
}

##' Compute information about slopes
##' @title Compute information about slopes
##' @param f Function
##' @param x Point
##' @param ... Additional arguments to \code{\link{bates_watts_D}}.
##' @seealso \code{\link{taylor2}}
##' @export
slope_info <- function(f, x, ...) {
  obj <- bates_watts_D(f, x, ...)
  list(x=x, fx=obj$f0, gr=gradient_D(obj), H=hessian_D(obj))
}

##' Create a 2nd order Taylor series (Laplace's method) approximation
##' to a function.
##' @title Approximate a function
##' @param obj Output from \code{\link{slope_info}}
##' @export
taylor2 <- function(obj) {
  if (!all(c("x", "fx", "gr", "H") %in% names(obj))) {
    stop("Invalid input")
  }
  x  <- obj$x
  fx <- obj$fx
  gr <- obj$gr
  H  <- obj$H
  len <- length(x)
  stopifnot(length(fx) == 1)
  stopifnot(length(gr) == len)
  stopifnot(is.matrix(H) && ncol(H) == len && nrow(H) == len)
  function(a) {
    if (length(a) != len) {
      stop("Expected input of length ", len)
    }
    v <- a - x
    drop(fx + gr %*% v + (v %*% H %*% v)/2)
  }
}
