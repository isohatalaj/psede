
#ifndef PSEDE_DIFF_H
#define PSEDE_DIFF_H

#include "psede.h"

/* Transform function signatures should conform to psede_transform_t,
 * inorder for them to be easily generalized to applying to
 * multidimensional grids. */

/**
 * Apply the differentiation operator on a vector or vectors of
 * Chebyshev mode coefficients; this is an O(size) operation per
 * vector
 */
int
psede_Tx_diff_mode_apply(double *x, int size, int stride, 
			 int howmany, int dist, void *params_is_null);

/**
 * Apply the integration operator on a vector or vectors of
 * Chebyshev mode coefficients; this is an O(size) operation per
 * vector
 */
int
psede_Tx_integ_mode_apply(double *x, int size, int stride, 
			  int howmany, int dist, void *params_is_null);

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
 * Apply a multidimensional differentiation to a matrix (`init == 0`),
 * or generate the a differentiation matrix (`init != 0`).
 */
/* int */
/* psede_Tx_diff_point_matrix_multi_0(double *diff_point, */
/* 				   int diff_dim, */
/* 				   int dims, const int *sizes, */
/* 				   int init, */
/* 				   psede_fct_t *fct); */

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
