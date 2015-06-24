context("grader")

test_that("bootstrap", {
  sc2_f <- function(x){
    n <- length(x)
    sum((1:n) * (exp(x) - x)) / n
  }
  n <- 5
  f <- make_f1(n)
  set.seed(1)
  x <- runif(n)
  expect_that(f(x), equals(sc2_f(x)))

  xx <- matrix(runif(n * 7), ncol=n)
  yy1 <- apply(xx, 1, sc2_f)
  yy2 <- apply(xx, 1, f)
  yy3 <- f(xx)
  expect_that(yy2, equals(yy1))
  expect_that(yy3, is_identical_to(yy2))
})

test_that("basic", {
  f <- make_f1(2L)
  set.seed(1)
  x <- runif(2L)

  pts <- gradient_points(x)
  expect_that(pts, is_a("gradient_points"))
  expect_that(dim(attr(pts, "dx")), equals(c(4L, 2L)))
  expect_that(attr(pts, "dim_y"), equals(c(2L, 8L)))
  expect_that(attr(pts, "n"), equals(2L))
  expect_that(attr(pts, "r"), equals(4L))

  ## These should be *really* similar:
  expect_that(gradient(f, x),
              equals(numDeriv::grad(f, x), tolerance=1e-16))

  ## general derivatives (bates watts D)
  D <- bates_watts_D_points(x)
  ans <- bates_watts_D(f, x)
  expect_that(bates_watts_D(f, x),
              equals(old_bates_watts_D(f, x), tolerance=1e-16))

  cmp <- numDeriv::genD(f, x)
  expect_that(ans$D, equals(ans$D, tolerance=1e-16))

  expect_that(ans, equals(cmp[names(ans)], tolerance=1e-16))
})

test_that("surface", {
  set.seed(10)
  f <- make_f2(4L)
  x <- runif(4L)

  ## These should be *really* similar:
  expect_that(gradient(f, x),
              equals(numDeriv::grad(f, x), tolerance=1e-16))

  ## general derivatives (bates watts D)
  D <- bates_watts_D_points(x)
  ans <- bates_watts_D(f, x)
  expect_that(bates_watts_D(f, x),
              equals(old_bates_watts_D(f, x), tolerance=1e-16))

  cmp <- numDeriv::genD(f, x)
  expect_that(ans$D, equals(ans$D, tolerance=1e-16))

  expect_that(ans, equals(cmp[names(ans)], tolerance=1e-16))

  obj <- slope_info(f, x)
  fa <- taylor2(obj)

  ## This should be identical:
  expect_that(fa(x), is_identical_to(f(x)))

  ## Honestly surprised at this, but whatever:
  for (i in 1:10) {
    xx <- runif(length(x), -1, 1)
    expect_that(fa(xx), equals(f(xx), tolerance=1e-4))
  }
})
