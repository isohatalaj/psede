
#ifndef PSEDE_DIFF_H
#define PSEDE_DIFF_H

#include "psede.h"

/**
 * Apply the differentiation operator on a vector or vectors of
 * Chebyshev mode coefficients; this is an O(size) operation per
 * vector
 */
void
psede_Tx_diff_mode_apply(double *x, int size, int stride, 
			 int howmany, int dist);

/**
 * Apply the integration operator on a vector or vectors of
 * Chebyshev mode coefficients; this is an O(size) operation per
 * vector
 */
void
psede_Tx_integ_mode_apply(double *x, int size, int stride, 
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
 * Apply the integration operator on a vector or vectors of
 * Chebyshev extrema grid values; this is an O(size*log(size))
 * operation per vector
 */
int
psede_Tx_integ_point_apply(double *x, int size, int stride,
			   int howmany, int dist, psede_fct_t *fct);

/**
 * Apply the differentiation operator on a vector of multi-dimensional
 * Chebyshev extrema grid values. Function point data is in array `x`
 * with `dims` dimensions, row-major ie. last dimension fastest varying
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
 * Apply differentiation to a matrix (`init == 0`), or generate the a
 * differentiation matrix (`init != 0`).
 */
int
psede_Tx_diff_point_matrix(double *diff_point, int size, int stride,
			   int howmany, int dist,
			   int init,
			   psede_fct_t *fct);

/**
 * Apply a multidimensional differentiation to a matrix (`init == 0`),
 * or generate the a differentiation matrix (`init != 0`).
 */
int
psede_Tx_diff_point_matrix_multi_0(double *diff_point,
				   int diff_dim,
				   int dims, const int *sizes,
				   int init,
				   psede_fct_t *fct);

/**
 * Apply a multidimensional integration to a matrix (`init == 0`), or
 * generate the a differentiation matrix (`init != 0`).
 */
int
psede_Tx_integ_point_matrix_multi_0(double *integ_point,
				    int integ_dim,
				    int dims, const int *sizes,
				    int init,
				    psede_fct_t *fct);


#endif /* PSEDE_DIFF_H */
