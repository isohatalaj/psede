
#ifndef PSEDE_H
#define PSEDE_H

#include "psede_fct.h"
#include "psede_diff.h"

/**
 * Compute a 1D extrema grid (Chebyshev-Gauss-Lobatto
 * nodes). Calculates `n` nodes with output in array `x` and
 * sequential elements separated by `stride` elements.  Parameter
 * `howmany` selects number of arrays to fill and `dist` selects the
 * distance between such arrays.
 */
void
psede_Tx_nodes(double *x, int n, int stride, int howmany, int dist);

/**
 * Same as `psede_Tx_nodes` but computes a single node array with unit
 * stride.
 */
void
psede_Tx_nodes_0(double *x, int n);


/**
 * Compute a multi-dimensinal extrema grid (Chebyshev-Gauss-Lobatto
 * nodes). Number of dimensions is `dims`, with desired number of grid
 * points for each dimension in `n`. The output array `x`, minimum
 * size `dims*n[0]*n[1]*...*n[dims-1]`, will contain the grid points
 * in dense "row-major" order, ie. the last dimension varies fastest
 * in memory.
 */
void
psede_Tx_nodes_multi_0(double *x, int dims, const int *n);


#endif /* PSEDE_H */
