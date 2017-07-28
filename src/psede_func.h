
#ifndef PSEDE_FUNC_H
#define PSEDE_FUNC_H

#include "psede_multi.h"
#include "psede_collocation.h"

typedef int (psede_func_t)(double x, double *y, void *params);

/**
 * Initialize a transform to be an in-place function application.
 */
int
psede_init_func(const psede_colloc_t *colloc,
		psede_transf_t *transform,
		psede_func_t *func,
		void *func_params);



/* ******************************************************************************** */
/* DEPRECIATED --> */

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
 * TODO: This needs to go!
 */
int
psede_function_multiply_multi_0(double *mat,
				double (*func)(const double *x, void *params),
				void *params,
				double **nodes,
				int dims, const int *sizes,
				int init);


/* ******************************************************************************** */

typedef struct {
  double (*func)(const double *x, void *func_params);
  void *func_params;
  int *index;
  double *x;
  const double **nodes;  
} psede_func_params_t;

psede_func_params_t *
psede_func_params_alloc(const double **nodes);

#endif /* PSEDE_FUNC_H */
