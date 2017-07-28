
#ifndef PSEDE_FCT_H
#define PSEDE_FCT_H

#include "psede_multi.h"

#include <fftw3.h>

/* Transform function signatures should conform to psede_transf_call_t,
 * for easy extension to multidimensional case.
 */

/**
 * Compute a 1D extrema grid (Chebyshev-Gauss-Lobatto
 * nodes). Calculates `n` nodes with output in array `x` and
 * sequential elements separated by `stride` elements.  Parameter
 * `howmany` selects number of arrays to fill and `dist` selects the
 * distance between such arrays. `params` is ignored.
 */
int
psede_Tx_nodes(double *x, int n, int stride, int howmany, int dist,
	       void *params);

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
 * in memory with no gaps between consequtive rows.
 */
void
psede_Tx_nodes_multi_0(double *x, int dims, const int *n);

/**
 * Fast Chebyshev transform workspace
 */
typedef struct {
  int size;
  int stride;
  int howmany;
  int dist;
  int align;
  fftw_plan plan;
} psede_fct_t;

/**
 * Allocate a Fast Chebyshev Transform object
 */
psede_fct_t*
psede_fct_alloc();

/**
 * Free a Fast Chebyshev Transform object
 */ 
void
psede_fct_free(psede_fct_t *self);

/**
 * Allocate memory best suited for a Fast Chebyshev Transform, in
 * units of sizeof(double). Presently this is just a wrapper for
 * FFTW's aligned allocation functions, given to present a dependency
 * independent API (ie. to accommodate different FFT backends).
 * Arrays allocated via this function must be freed using
 * psede_fct_free_array.
 */

double*
psede_fct_alloc_array(int size);

/**
 * Free an array alloceted using psede_fct_alloc_array
 */
void
psede_fct_free_array(double *array);

/**
 * Apply a Fast Chebyshev Transform to a vector or vectors; converts
 * point values to mode coefficients. Return zero for success.
 */
int
psede_Tx_fct_apply(double *x, int size, int stride,
		   int howmany, int dist, psede_fct_t *self);

/**
 * Apply an Inverse Fast Chebyshev Transform to a vector or vectors;
 * converts mode coefficients to point values. Return zero for
 * success.
 */
int
psede_Tx_fct_apply_inv(double *x, int size, int stride,
		       int howmany, int dist, psede_fct_t *self);


/* ******************************************************************************** */
/* Transformer interface to FCTs. */

/**
 * Initialize a transf to the Chebyshev extrema nodes
 * transform (i.e. the transform that sets an input array
 * to the grid points). Parameter `params` is ignored.
 */
int
psede_init_Tx_nodes(psede_transf_t *, void *params);


/**
 * Initialize a transf to the Chebyshev extrema points-to-modes
 * transform. Parameter `params` is ignored.
 */
int
psede_init_Tx_fct(psede_transf_t *, void *params);

/**
 * Initialize a transf to the Chebyshev extrama modes-to-points
 * transform. Parameter `params` is ignored.
 */
int
psede_init_Tx_ifct(psede_transf_t *, void *params);


#endif /* PSEDE_FCT_H */
