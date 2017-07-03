
#ifndef PSEDE_DIFF_H
#define PSEDE_DIFF_H

#include "psede_fct.h"

/**
 * Apply the differentiation operator on a vector or vectors of
 * Chebyshev mode coefficients; this is an O(size) operation per
 * vector
 */
void
psede_Tx_diff_mode_apply(double *x, int size, int stride, 
			 int howmany, int dist);

/**
 * Apply the differentiation operator on a vector or vectors of
 * Chebyshev extrema grid values; this is an O(size*log(size))
 * operation per vector
 */
int
psede_Tx_diff_point_apply(double *x, int size, int stride,
			  int howmany, int dist, psede_fct_t *fct);

/**
 * Apply the differentiation operator on a vector of multi-dimensional
 * Chebyshev extrema grid values. Function point data is in array `x`
 * with `dims` dimensions, row-major or last dimension fastest varying
 * layout, with `sizes` points per each dimension. Differentation
 * computed along `diff_dim` dimension.
 */
int
psede_Tx_diff_point_apply_multi_0(double *x, int diff_dim,
				  int dims, const int *sizes, psede_fct_t *fct);

/**
 * Initialize given array to mode space differentiation matrix
 */
void
psede_Tx_diff_mode_matrix(double *diff_mode, int size, int stride,
			  int howmany, int dist);

/**
 * Initialize given array to mode space differentiation matrix
 */
int
psede_Tx_diff_point_matrix(double *diff_point, int size, int stride,
			   int howmany, int dist,
			   psede_fct_t *fct);


#endif /* PSEDE_DIFF_H */
