/**
 * @file cps_ode.h Ordinary differential equation solver routines.
 */

#ifndef CPS_ODE_H
#define CPS_ODE_H

#include "cps.h"

/**
 * Linear ODE function type. 
 */
typedef int (*cps_linode_fun_t)(double x, 
				double *A, double *Q, double *C, 
				void *params);

/**
 *  CGL pseudospectral ODE solver object. 
 */
typedef struct {
  int n_size;    /**< Number of points, also order of the method */
  int n_dim;     /**< Dimension of the ODE system */
  
  double *x;     /**< Node x-positions, size n_size */
  double *y;     /**< For <= n_dim functions f at xs, size n_size*n_dim */
  double *M;     /**< ODE operator matrix, size (n_size*n_dim)^2 */
  double *D;     /**< CGL derivative matrix, size (n_size*n_size) */

  double *F;     /**< Derivatives at given point, size n_dim */
  double *A;     /**< Deriv. matrix or Jacobian, size n_dim*n_dim */
  double *Q;     /**< Blah */


  fftw_plan ply; /**< FFTW plan for fast C'shev transform at vector y,
		   size n_size */
  
  int *ipiv;    /**< Work array for LAPACK (pivot indices) */
} cps_ode_t;


/**
 * Allocate ODE solver structure. 
 */
cps_ode_t *
cps_ode_make(int n_size, int n_dim);

/**
 * Free ODE solver structure. 
 */
void
cps_ode_free(cps_ode_t *self);

/**
 * Solve a linear ODE boundary or initial value problem.
 * The equations are taken to have the form 
 * \f[ Q(x)Y'(x) = A(x)Y(x) + C(x), \f]
 * 
 *
 * @param x0 Solution range lower bound.
 * @param x1 Range upper bound.
 * @param fun Function that evaluates the ODE system.

 * @param nbcs Number of boundary conditions. Note that although
 * typically this equals the dimension (the differential order) of the
 * problem, this may not always be the case.

 * @param bcs Boundary conditions, left-hand side. If the dimension of
 * the problem is ndim, then this is a nbcs-by-2*ndim matrix. Each row
 * of the matrix gives one boundary condition, with each element
 * corresponding to a coefficient of
 *
 * @param gammas Boundary conditions, right-hand side.
 * @param ls Positions for the boundary conditions in matrix form
 * equation.
 */
int
cps_linode_solve(cps_ode_t *self, 
		 double x0, double x1,
		 cps_linode_fun_t fun, void *params,
		 int nbcs,
		 const double *bcs, 
		 const double *gammas,
		 const int *ls);

#endif /* CPS_ODE_H */
