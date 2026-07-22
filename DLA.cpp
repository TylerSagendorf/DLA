// [[Rcpp::depends(dqrng)]]
#include <Rcpp.h>
#include <dqrng.h> // fast random sampling

#define MIN(x, y) (x) < (y) ? (x) : (y)
#define MAX(x, y) (x) > (y) ? (x) : (y)

// New particles will start from a position that is at most this many units from
// the particle that is currently closest to the starting side
#define LOCAL_REGION 4

// Helpers for matrix indexing
#define TOP(pout, NR, x_coord, y_coord) pout[(x_coord) * NR + y_coord - 1] != NA_INTEGER
#define BOTTOM(pout, NR, x_coord, y_coord) pout[(x_coord) * NR + y_coord + 1] != NA_INTEGER
#define LEFT(pout, NR, x_coord, y_coord) pout[((x_coord) - 1) * NR + y_coord] != NA_INTEGER
#define RIGHT(pout, NR, x_coord, y_coord) pout[(x_coord + 1) * NR + y_coord] != NA_INTEGER
#define TOPRIGHT(pout, NR, x_coord, y_coord) pout[(x_coord + 1) * NR + y_coord - 1] != NA_INTEGER
#define BOTTOMRIGHT(pout, NR, x_coord, y_coord) pout[(x_coord + 1) * NR + y_coord + 1] != NA_INTEGER
#define TOPLEFT(pout, NR, x_coord, y_coord) pout[(x_coord - 1) * NR + y_coord - 1] != NA_INTEGER
#define BOTTOMLEFT(pout, NR, x_coord, y_coord) pout[(x_coord - 1) * NR + y_coord + 1] != NA_INTEGER

// Random sampling
#define RAND_2 dqrng::dqsample_int(2, 1, true)[0] // sample 0, 1
#define RAND_3 dqrng::dqsample_int(3, 1, true, R_NilValue, -1)[0] // sample -1, 0, 1


// [[Rcpp::export(.Cpp_DLA)]]
SEXP DLA(SEXP Rnr,
         SEXP Rnc,
         SEXP x,
         SEXP y,
         SEXP Rn_particles,
         Rcpp::Nullable<Rcpp::IntegerVector> seed) {

  const int NR = INTEGER(Rnr)[0];
  const int NC = INTEGER(Rnc)[0];
  const int n_particles = INTEGER(Rn_particles)[0];

  const int *px = INTEGER(x);
  const int *py = INTEGER(y);

  const int x_len = Rf_length(x);

  SEXP out = PROTECT(Rf_allocMatrix(INTSXP, NR, NC));
  int *pout = INTEGER(out);

  for (int i = 0; i < NR * NC; ++i) {
    pout[i] = NA_INTEGER;
  }

  int min_x = NC - 1;
  int max_x = 0;
  int min_y = NR - 1;
  int max_y = 0;

  int particle_order = 0;

  for (int i = 0; i < x_len; ++i) {
    const int xi = px[i] - 1;
    const int yi = py[i] - 1;

    min_x = MIN(min_x, xi);
    max_x = MAX(max_x, xi);
    min_y = MIN(min_y, yi);
    max_y = MAX(max_y, yi);

    pout[xi * NR + yi] = particle_order;
  }

  ++particle_order;

  dqrng::dqset_seed(seed);

  for (int p = 0; p < n_particles; ++p) {

    // TODO start particles from any side - not just the top. Use a switch
    // statement for the sides as integers?

    // TODO if a side is saturated, stop adding particles to that side. Stop
    // adding particles if all sides are saturated (check that the number of
    // saturated sides is 4).

    int x_coord = dqrng::dqsample_int(NC, 1, true)[0];
    const int y_start = MAX(0, min_y - LOCAL_REGION);
    int y_coord = y_start;

    bool freeze_particle = false;

    while (true) {

      // Reset particles that get too far away from the starting y coordinate
      y_coord = y_coord < y_start - LOCAL_REGION ? y_start : y_coord;

      if (x_coord == 0) {

        if (y_coord == 0) {

          /* Upper-left corner */

          freeze_particle =
            BOTTOM(pout, NR, x_coord, y_coord) ||
            RIGHT(pout, NR, x_coord, y_coord) ||
            BOTTOMRIGHT(pout, NR, x_coord, y_coord);

          if (freeze_particle) break;

          y_coord += RAND_2;

        } else if (y_coord == NR - 1) {

          /* Bottom-left corner */

          freeze_particle =
            TOP(pout, NR, x_coord, y_coord) ||
            TOPRIGHT(pout, NR, x_coord, y_coord) ||
            RIGHT(pout, NR, x_coord, y_coord);

          if (freeze_particle) break;

          y_coord -= RAND_2;

        } else {

          /* Left side */

          freeze_particle =
            TOP(pout, NR, x_coord, y_coord) ||
            BOTTOM(pout, NR, x_coord, y_coord) ||
            TOPRIGHT(pout, NR, x_coord, y_coord) ||
            RIGHT(pout, NR, x_coord, y_coord) ||
            BOTTOMRIGHT(pout, NR, x_coord, y_coord);

          if (freeze_particle) break;

          y_coord += RAND_3;

        }

        x_coord += RAND_2;

      } else if (x_coord == NC - 1) {

        if (y_coord == 0) {

          /* Upper-right corner */

          freeze_particle =
            LEFT(pout, NR, x_coord, y_coord) ||
            BOTTOMLEFT(pout, NR, x_coord, y_coord) ||
            BOTTOM(pout, NR, x_coord, y_coord);

          if (freeze_particle) break;

          y_coord += RAND_2;

        } else if (y_coord == NR - 1) {

          /* Bottom-right corner */

          freeze_particle =
            TOPLEFT(pout, NR, x_coord, y_coord) ||
            LEFT(pout, NR, x_coord, y_coord) ||
            TOP(pout, NR, x_coord, y_coord);

          if (freeze_particle) break;

          y_coord -= RAND_2;

        } else {

          /* Right side */

          freeze_particle =
            TOPLEFT(pout, NR, x_coord, y_coord) ||
            LEFT(pout, NR, x_coord, y_coord) ||
            BOTTOMLEFT(pout, NR, x_coord, y_coord) ||
            TOP(pout, NR, x_coord, y_coord) ||
            BOTTOM(pout, NR, x_coord, y_coord);

          if (freeze_particle) break;

          y_coord += RAND_3;

        }

        x_coord -= RAND_2;

      } else {

        if (y_coord == 0) {

          /* Top side */

          freeze_particle =
            LEFT(pout, NR, x_coord, y_coord) ||
            BOTTOMLEFT(pout, NR, x_coord, y_coord) ||
            BOTTOM(pout, NR, x_coord, y_coord) ||
            RIGHT(pout, NR, x_coord, y_coord) ||
            BOTTOMRIGHT(pout, NR, x_coord, y_coord);

          if (freeze_particle) break;

          y_coord += RAND_2;

        } else if (y_coord == NR - 1) {

          /* Bottom side */

          freeze_particle =
            TOPLEFT(pout, NR, x_coord, y_coord) ||
            LEFT(pout, NR, x_coord, y_coord) ||
            TOP(pout, NR, x_coord, y_coord) ||
            TOPRIGHT(pout, NR, x_coord, y_coord) ||
            RIGHT(pout, NR, x_coord, y_coord);

          if (freeze_particle) break;

          y_coord -= RAND_2;

        } else {

          /* Center */

          freeze_particle =
            TOPLEFT(pout, NR, x_coord, y_coord) ||
            LEFT(pout, NR, x_coord, y_coord) ||
            BOTTOMLEFT(pout, NR, x_coord, y_coord) ||
            TOP(pout, NR, x_coord, y_coord) ||
            BOTTOM(pout, NR, x_coord, y_coord) ||
            TOPRIGHT(pout, NR, x_coord, y_coord) ||
            RIGHT(pout, NR, x_coord, y_coord) ||
            BOTTOMRIGHT(pout, NR, x_coord, y_coord);

          if (freeze_particle) break;

          y_coord += RAND_3;

        }

        x_coord += RAND_3;

      }

    }

    min_y = MIN(min_y, y_coord);

    pout[x_coord * NR + y_coord] = particle_order++;
  }

  UNPROTECT(1);

  return out;

}


/*** R
library(dqrng)
library(ComplexHeatmap)

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

nr <- 2^8
nc <- 2^8

# The bottom row of the output matrix will be filled with particles to start
x <- seq_len(nc)
y <- rep.int(nr, length(x))
n_particles <- 2^12

# bench::mark(
#   DLA(nr = nr, nc = nc, y = y, x = x, n_particles = n_particles, seed = 0L),
#   iterations = 30L
# )

res <- DLA(nr = nr, nc = nc, y = y, x = x, n_particles = n_particles, seed = 999L)

# Heatmap ----
grid_size <- unit(2, "pt")

Heatmap(
  matrix = res,
  col = circlize::colorRamp2(
    breaks = range(res, na.rm = TRUE),
    colors = c("black", "white")
  ),
  na_col = "lightblue",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  height = grid_size * nr,
  width = grid_size * nc,
  heatmap_legend_param = list(
    title = "Arrival order"
  )
)
*/
