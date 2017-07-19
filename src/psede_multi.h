
#ifndef PSEDE_MULTI_H
#define PSEDE_MULTI_H

#include "psede.h"

/**
 * One-dimensional transform type. Here, `target` is on input the
 * array to be transformed, on output it will be the transformed data.
 * Parameter `size` is the number of elements in each array to be
 * transformed, with consequtive elements separated by `stride`
 * units. The transform is repeatedly applied to `howmany` vectors,
 * which are in turn separated by `dist` elements (to be very clear,
 * the distance between the initial elements of two consequtive arrays
 * is `dist`).
 */
typedef int (psede_transform_t)(double *target, int size, int stride,
				int howmany, int dist, void *params);

/**
 * Convenience function for applying a 1D transform on rows or columns
 * of a tight-packed square matrix. The matrix is assumed to be
 * row-major `size`-by-`size`, lead-dimension `size`, with the
 * transform applied column-wise if `trans` is 0. If `init` is
 * non-zero, the target array is initialized to identity matrix.
 */
int
psede_apply_matrix_0(double *target, int size, int trans, int init,
		     psede_transform_t *transform, void *params);

/**
 * Apply a one-dimensional transformation to a multi-dimensional
 * row-major array along a selected axis.
 */
int
psede_apply_multi(double *multi_target, int n_dimensions, const int *sizes,
		  psede_transform_t *transform, int target_dimension, void *params);


int
psede_apply_multi_matrix_0(double *multi_target, int n_dimensions, const int *sizes,
			   int transp, int init,
			   psede_transform_t *transform,
			   int target_dimension, void *params);

#endif /* PSEDE_MULTI_H */
