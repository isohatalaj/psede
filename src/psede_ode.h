
#ifndef PSEDE_ODE_H
#define PSEDE_ODE_H

#include "psede.h"

/**
 * Linear ODE function type. 
 */
typedef int (*psede_ode_linear_fun_t)(double x, 
				      double *A, double *Q, double *C, 
				      void *params);

/**
 * ODE solver work structure
 */
typedef struct {
  int size;
  int dim;

  double *x;
  double *y;
  double *D;

  psede_fct_t *fct;

  int work_size;
  double *work;
  psede_linsolve_work_t *linsolve_work;
} psede_ode_t;

/**
 * Allocate ODE solver object and initialize it for solving `dim`
 * dimensional ordinary differential equation system, with `size`
 * point approximation, ie `size - 1` order polynomial approximation
 */ 
psede_ode_t*
psede_ode_alloc(int size, int dim);

/**
 * Free ODE solver object
 */
void
psede_ode_free(psede_ode_t *self);

/**
 * Solve a generic linear ordinary differential equation. The equation
 * is assumed to be written in the form
 *
 *     Q(x)Y'(x) = A(x)Y(x) + C(x),
 *
 * where Y(x) is the function to be solved at x, a `dim` dimensional
 * vector; A(x) and Q(x) are `dim`-by-`dim` matrices, and C(x) is a
 * `dim` dimensional vector.
 *
 * Number of boundary conditions is `nbcs` (which may be greater or
 * less than `dim` -- the user is assumed to know what they're
 * doing). The conditions are given in the form B.(Y(x0), Y(x1)) =
 * gamma where B is a `nbcs`-by-2`dim` matrix and `gamma` is
 * `nbcs`-dimensional vector.
 */
int
psede_ode_linear_solve(psede_ode_t *self, 
		       double x0, double x1,
		       psede_ode_linear_fun_t fun, void *params,
		       int nbcs,
		       const double *B, 
		       const double *gamma,
		       const int *ls);


#endif
