##' @title Computes Bates and Watts D matrix
##' @param f Fucntion
##' @param x Point
##' @param eps See \code{numderiv::genD}
##' @param d See \code{numderiv::genD}
##' @param r See \code{numderiv::genD}
##' @param v See \code{numderiv::genD}
##' @param zero_tol See \code{numderiv::genD}
##' @export
bates_watts_D <- function(f, x, eps=1e-4, d=0.0001, r=4, v=2,
                          zero_tol=sqrt(.Machine$double.eps/7e-7)) {
  pts <- bates_watts_D_points(x, eps, d, r, zero_tol)
  y <- bates_watts_D_eval(f, pts)
  bates_watts_D_extrapolate(x, y, pts)
}

bates_watts_D_points <- function(x, eps=1e-4, d=0.0001, r=4, v=2,
                                 zero_tol=sqrt(.Machine$double.eps/7e-7)) {
  v <- 2 # required
  n <- length(x) # number of parameters (theta)

  h <- abs(d * x) + eps * (abs(x) < zero_tol)

  ## Points for computing first derivatives (same as gradient_points)
  pts <- vector("list", r * n)
  dx <- numeric(r * n)
  dim(pts) <- dim(dx) <- c(r, n)

  ## For second derivatives:
  pts2 <- dx2 <- vector("list", r * n * n)
  dim(pts2) <- dim(dx2) <- c(r, n, n)

  ## Build matrices of points:
  ## Note that the dx2 matrix is lower triangle (diagonal exclusive)
  ## only.
  idx <- seq_len(n)
  for (k in seq_len(r)) {
    for (i in seq_len(n)) {
      dx_i <- h * (i == idx)
      pts[[k, i]] <- rbind(x + dx_i, x - dx_i)
      for (j in seq_len(i - 1L)) {
        dx_j <- h * (j == idx)
        pts2[[k, j, i]] <- rbind(x + dx_i + dx_j,
                                 x - dx_i - dx_j)
        dx2[[k, j, i]] <- c(h[i], h[j])
      }
    }
    dx[k, ] <- h
    h <- h / v
  }

  ret1 <- x
  ret2 <- do.call("rbind", pts)
  ret3 <- do.call("rbind", pts2)

  ret <- rbind(ret1, ret2, ret3, deparse.level=0)

  attr(ret, "n") <- n
  attr(ret, "r") <- r
  attr(ret, "dx") <- dx
  attr(ret, "dx2") <- dx2
  attr(ret, "idx") <- rep(1:3, c(1L, nrow(ret2), nrow(ret3)))
  attr(ret, "dim_y") <- c(2L, dim(pts))
  class(ret) <- "bates_watts_D_points"
  ret
}

bates_watts_D_eval <- function(f, pts) {
  y <- f(pts)
  idx <- attr(pts, "idx")
  list(y[[which(idx == 1L)]],
       array(y[idx == 2L], attr(pts, "dim_y")),
       y[idx == 3L])
}

bates_watts_D_extrapolate <- function(x, y, pts) {
  len_f <- 1L # assume scalar output

  n <- attr(pts, "n")
  r <- attr(pts, "r")
  dx <- attr(pts, "dx")
  dx2 <- attr(pts, "dx2")

  f0 <- y[[1L]]
  f1 <- y[[2L]]
  f2 <- y[[3L]]

  ## Given that we're not allowing multivalue functions here, does
  ## this make sense?
  D <- matrix(0, len_f, (n * (n + 3)) / 2)
  Daprox <- matrix(0, len_f, r)
  Hdiag  <- matrix(0, len_f, n)
  Haprox <- matrix(0, len_f, r)

  for (i in seq_len(n)) { # over parameters
    for (k in seq_len(r)) {
      Daprox[, k] <- (f1[1, k, i] -          f1[2, k, i]) / (2 * dx[k, i])
      Haprox[, k] <- (f1[1, k, i] - 2 * f0 + f1[2, k, i]) / dx[k, i]^2
    }
    for (m in seq_len(r - 1L)) {
      for (k in seq_len(r - m)) {
        Daprox[, k] <- (Daprox[,  k + 1L] * (4^m) - Daprox[, k]) / (4^m - 1)
        Haprox[, k] <- (Haprox[,  k + 1L] * (4^m) - Haprox[, k]) / (4^m - 1)
      }
    }
    D[, i]     <- Daprox[, 1]
    Hdiag[, i] <- Haprox[, 1]
  }

  u <- n
  idx <- 1L
  for (i in seq_len(n)) {
    for (j in seq_len(i)) {
      u <- u + 1L
      if (i == j) {
        D[, u] <- Hdiag[, i]
      } else {
        for (k in seq_len(r)) {
          tmp_f <- f2[c(idx, idx + 1L)]
          tmp_h <- dx2[[k, j, i]]
          Daprox[, k] <-
            (tmp_f[[1]] - 2 * f0 + tmp_f[[2]] -
             Hdiag[, i] * tmp_h[1L]^2 -
             Hdiag[, j] * tmp_h[2L]^2) /
               (2 * tmp_h[1L] * tmp_h[2L])
          idx <- idx + 2L
        }
        for (m in seq_len(r - 1L)) {
          for (k in seq_len(r - m)) {
            Daprox[, k] <- (Daprox[, k+1] * (4^m) - Daprox[, k]) / (4^m - 1)
          }
        }
        D[, u] <- Daprox[, 1]
      }
    }
  }

  list(D=D, p=n, f0=f0, x=x)
}
