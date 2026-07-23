for (pkg in c("dqrng", "Rcpp")) {
  if (!require(pkg, quietly = TRUE, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

library(dqrng)

Rcpp::sourceCpp("DLA.cpp")

DLA <- function(nr, nc, x, y, n_particles, seed = NULL) {

  if (!is.vector(x, "numeric") ||
      !is.vector(y, "numeric") ||
      length(x) == 0L ||
      length(x) != length(y)) {
    stop("'x' and 'y' must be integer vectors with the same length > 0.")
  }

  if (!is.vector(nr, "numeric") ||
      !is.vector(nc, "numeric") ||
      length(nr) != 1L ||
      length(nc) != 1L ||
      nr < 2L ||
      nc < 2L) {
    stop("'nr' and 'nc' must be length-1 integer vectors >= 2.")
  }

  if (!is.vector(n_particles, "numeric") ||
      length(n_particles) != 1L ||
      n_particles < 1) {
    stop("'n_particles' must be a length-1 integer vector >= 1.")
  }

  if (!is.null(seed)) {
    if (!is.vector(seed, "numeric") || length(seed) != 1L) {
      stop("'seed' must be NULL or a length-1 integer vector.")
    }

    seed <- as.integer(seed)
  }

  nr <- as.integer(nr)
  nc <- as.integer(nc)
  n_particles <- as.integer(n_particles)

  storage.mode(x) <- storage.mode(y) <- "integer"

  x_oob <- which(x < 1L | x > nc)

  if (length(x_oob)) {
    stop("All elements of 'x' must be between 1 and 'nc'.")
  }

  y_oob <- which(y < 1L | y > nr)

  if (length(y_oob)) {
    stop("All elements of 'y' must be between 1 and 'nr'.")
  }

  mat <- .Cpp_DLA(nr, nc, x, y, n_particles, seed)

  return(mat)
}
