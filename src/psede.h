
#ifndef PSEDE_H
#define PSEDE_H

#include "psede_fct.h"
#include "psede_diff.h"

/**
 * Compute a 1D extrema grid (Chebyshev-Gauss-Lobatto nodes)
 */
void
psede_nodes(double *x, int n, int stride);

#endif /* PSEDE_H */
