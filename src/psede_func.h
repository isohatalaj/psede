
#ifndef PSEDE_FUNC_H
#define PSEDE_FUNC_H

#include "psede.h"

/**
 * Apply a multi-dimensional function multiplication onto a matrix.
 * Here, `mat` is the target matrix (if `init` is true, this is
 * initialized to the identity matrix), `func` is the function to
 * apply where `x` is a vector of argument values and `params` are
 * parameters to that function, `nodes` is an array of collocation
 * point arrays, `dims` is the number of dimensions with `sizes`
 * giving the number of points along each dimension. Returns
 * zero on success.
 *
 */
int
psede_function_multiply_multi_0(double *mat,
				double (*func)(const double *x, void *params),
				void *params,
				double **nodes,
				int dims, const int *sizes,
				int init);



#endif /* PSEDE_FUNC_H */
