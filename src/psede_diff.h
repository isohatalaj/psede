
#ifndef PSEDE_DIFF_H
#define PSEDE_DIFF_H

#include "psede_fct.h"

/* Transform function signatures should conform to psede_transf_call_t,
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

/* ******************************************************************************** */
/* Transformer API functions. The above routines can be freely used to
 * do low-level transform operations, however, using them via
 * transf objects can greatly simplify code that requires
 * composition and linear combination of multiple transformations.
 */

/**
 * Initialize a given transform to be differentiation operator for
 * Chebyshev extrema grid (point space). The parameter `order`
 * specifies the order of the differentiation; negative orders signify
 * integration. Implicitly uses the standard FCT backend to transfer
 * between mode and point representations. `params` is ignored.
 */
int
psede_init_Tx_diff(psede_transf_t *transf, int order, void *params);


#endif /* PSEDE_DIFF_H */
