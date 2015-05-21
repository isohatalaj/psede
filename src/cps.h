/**
 * @file cps.h Chebyshev pseudspectral differential equation solver
 * routines. These routines are intended as general purpose, simple
 * ODE/PDE boundary/initial value problem solvers. 
 *
 * The solvers we give are following types of problems:
 * 1. Linear systems of the form
 * \f[  
 *        Q(x)Y'(x) = A(x)Y(x) + C(x),   
 * \f]
 * where \f$Y(x), Y'(x), C(x)\f$ are \f$R^n\f$ vectors and
 * \f$Q(x)\f$, \f$A(x)\f$ are \f$R^{n \times n}\f$ matrices
 * over a finite interval \f$I = [x_0, x_1]\f$, with Robin boundary
 * conditions.
 * 2. Non-linear systems of the form
 * \f[  
 *        F(x, Y(x), Y'(x)) = 0,     
 * \f]
 * where \f$F\f$ is some function \f$[x_0,x_1] \otimes R^n \otimes R^n
 * \to R^n\f$ that is at least once differentiable in all of its
 * arguments.
 *
 */

#ifndef CPS_H
#define CPS_H

#include <stdlib.h>
#include <fftw3.h>

/**
 * Fill a vector with the nth order CGL nodes. 
 */
void
cps_nodes(double *x,  /**< Output array to contain the nodes */
	  int n,      /**< Length of array `x` */
	  int stride  /**< Array `x` stride */
	  );

/**
 * Set a vector of length `n` to the \f$i\f$th unit vector, \f$u_j =
 * \delta_{ij}\f$, that is `u[j] = i == j ? 1.0 : 0.0`.
 */
void
cps_unit(int i,      /**< Index of the non-zero element in output */
	 double *u,  /**< Output array, possibly strided  */
	 int n,      /**< Length of the array `u` */
	 int stride  /**< Array `u` stride  */
	 );

/**
 * Generate an FFTW plan needed for fast Chebyshev transforms, the
 * generated plan will be associated with the array `x`.
 *
 * @return The newly constructed FFTW plan, or `NULL` on failure.
 */
fftw_plan
cps_plan_T(double *x, /**< Generate the plan for this array */
	   int n,     /**< Length of array `x` */
	   int stride /**< Array `x` stride */
	   );

/**
 * Apply the differentiation operator to vector of Chebyshev mode
 * coefficients.
 *
 * @param x Array containing Chebyshev mode coefficients.
 * @param n Length of array `x`.
 * @param stride Array `x` stride.
 */
void
cps_apply_D_mod(double *x, int n, int stride);

/**
 * Transform a vector of CGL node values to mode coefficients. Uses a
 * fast Chebyshev transform that is performed by executing a
 * previously prepared FFTW plan. The array x must be associated with
 * that plan.
 *
 * @param x Input array containing function values at CGL nodes. Must
 * be associated with the given FFTW plan.
 * @param n Length of array `x`.
 * @param stride Array `x` stride.
 * @param fftw_plan FFTW associated with the array `x`. Should
 * be constructed using the function `cps_plan_T`.
 *
 * @see cps_plan_T
 */
void
cps_apply_T(double *x, int n, int stride, const fftw_plan plan);

/**
 *  Transform a vector of mode coefficients to CGL node values. 
 */
void
cps_apply_T_inv(double *x, int n, int stride, const fftw_plan plan);

#endif /* CPS_H */
