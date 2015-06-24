make_f1 <- function(n) {
  force(n)
  function(x) {
    if (is.matrix(x)) {
      if (ncol(x) != n) {
        stop("Wrong number of columns")
      }
    } else {
      if (length(x) != n) {
        stop("Wrong number of elements")
      }
      x <- rbind(x, deparse.level=0)
    }
    rowSums(rep(seq_len(n), each=nrow(x)) * (exp(x) - x)) / n
  }
}

make_f2 <- function(n) {
  vcv <- diag(seq_len(n))
  set.seed(1)
  vcv[upper.tri(vcv)] <- runif(n * (n - 1) / 2, -.5, .5)
  vcv[lower.tri(vcv)] <- t(vcv)[lower.tri(vcv)]
  function(x) {
    if (!(length(x) == n || ncol(x) == n)) {
      stop("Invalid size x")
    }
    mvtnorm::dmvnorm(x, sigma=vcv, log=TRUE)
  }
}
