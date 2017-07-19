/**
 * @file lintest.c An example linear differential equation problem.
 * The sample problem we consider is \f[ y''(x) + p y'(x) + q y(x) =
 * 0, \f] which translated to our canonical form \f$QY' = AY\f$, gives
 * \f[ 
 *    Y = \left(\begin{array}{c} y \\ y' \end{array}\right),
 *    \qquad
 *    Q = \left(\begin{array}{cc} 1 & 0 \\ 0 & 1 \end{array}\right),
 *    \qquad
 *    A = \left(\begin{array}{cc} 0 & 1 \\ -q & -p \end{array}\right).
 * \f]
 * At boundaries, we use some generic Robin/periodic conditions
 * \f[
 *    \left(\begin{array}{cccc}
 *      b_{y_0}^1 & b_{y'_0}^1 & b_{y_1}^1 & b_{y'_1}^1 \\
 *      b_{y_0}^2 & b_{y'_0}^2 & b_{y_1}^2 & b_{y'_1}^2 
 *    \end{array}\right)
 *    \left(\begin{array}{c} 
 *      y(x_0) \\ y'(x_0) \\ y(x_1) \\ y'(x_1) 
 *    \end{array}\right)
 *    =
 *    \left(\begin{array}{c} 
 *      \gamma_0^1 \\ \gamma_1^1 \\ \gamma^1_0 \\ \gamma^1_1
 *    \end{array}\right).
 * \f]
 * This specification allows us to impose arbitrary initial value /
 * boundary value / periodic conditions.  See the code for #main for
 * actual parameter choices.
 *  
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "psede.h"
#include "psede_ode.h"

/** Test linear differential equation function.  Evaluates the
 *  differential equation in the form required by #cps_linode_solve,
 *  ie has the type #cps_linode_fun_t. 
 *
 *  @see cps_linode_solve
 *  @see cps_linode_fun_t
 */
int
testfun(double x, double *A, double *Q, double *C, void *vpars)
{
  double *par = vpars;
  double p = par[0], q = par[1];

  A[0] = 0;   A[1] = 1;
  A[2] = -q;  A[3] = -p;

  Q[0] = Q[3] = 1.0;
  Q[1] = Q[2] = 0.0;

  C[0] = 0;
  C[1] = 0;

  return 0;
}

int
main()
{
  int status = 0;

  psede_ode_t *odesolver;

  int n_dim = 2;	/*< Dimension of our ODE proble, here two */
  int n_size = 16;      /*< Order of the solver to use */
  double x0 = 0.0;      /*< Solution interval start */
  double x1 = 1.0;      /*< Solution interval end */

  /* ODE system parameters, here vector of (p, q) */
  double pars[2] = {1.0, -1.0}; 

  /* Boundary conditions */
  int nbcs = 2;         /*< Number of boundary conditions */
  double B[8] = {1.0, 0.0, 0.0, 0.0,
		 0.0, 1.0, 0.0, 0.0};
  double gamma[2] = {0.0, 1.0};
  int ls[2] = {1, 2};

  /* Analytic solution is C1*exp(b1) + C2*exp(b2) where: */
  double b1 = 0.5*(-pars[0] + sqrt(pars[0]*pars[0]-4*pars[1]));
  double b2 = 0.5*(-pars[0] - sqrt(pars[0]*pars[0]-4*pars[1]));
  double C2 = (gamma[1] - b1*gamma[0])/(b2 - b1);
  double C1 = gamma[0] - C2;

  printf("# Solving differential equation\n"
	 "#\n"
	 "#    y''(x) + p y'(x) + q y(x) = 0\n"
	 "#\n"
	 "# over x in [x0, x1] with boundary/initial value conditions\n"
	 "#\n"
	 "#    b00 y(x0) + b01 y'(x0) + b02 y(x1) + b03 y'(x1) = gamma0 \n"
	 "#    b10 y(x0) + b11 y'(x0) + b12 y(x1) + b13 y'(x1) = gamma1 \n"
	 "#\n"
	 "# where\n"
	 "#\n"
	 "#    p      = %21.15lf, q      = %21.15lf\n"
	 "#\n"
	 "#    b00    = %21.15lf, b01    = %21.15lf, b02    = %21.15lf, b03    = %21.15lf\n"
	 "#    b10    = %21.15lf, b11    = %21.15lf, b12    = %21.15lf, b13    = %21.15lf\n"
	 "#\n"
	 "#    gamma0 = %21.15lf, gamma1 = %21.15lf\n"
	 "#\n"
	 "# Analytic solution is\n"
	 "#\n"
	 "#    y = C1 exp(beta1) + C2 exp(beta2)\n"
	 "#\n"
	 "# where\n"
	 "#\n"
	 "#    C1    = %21.15lf, C2    = %21.15lf\n"
	 "#    beta1 = %21.15lf, beta2 = %21.15lf\n"
	 "#\n"
	 , pars[0], pars[1]
	 , B[0], B[1], B[2], B[3]
	 , B[4], B[5], B[6], B[7]
	 , gamma[0], gamma[1]
	 , C1, C2
	 , b1, b2
	 );

  odesolver = psede_ode_alloc(n_size, n_dim);
  if (odesolver == NULL)
    {
      fprintf(stderr, "Failed creating solver object.\n");
      return 1;
    }

  status = psede_ode_linear_solve(odesolver,
  				  x0, x1,
  				  testfun, pars,
  				  nbcs, B, gamma, ls);
  if (status)
    {
      fprintf(stderr, "Linear ODE solver failed.\n");
      psede_ode_free(odesolver);
      return 2;
    }


  int i;
  double *xs = odesolver->x;
  double *ys = odesolver->y;

  for (i = 0; i < n_size; ++i)
    {
      double x = 0.5*(x1-x0)*xs[i] + 0.5*(x0+x1);
      double ya = C1*exp(b1*x) + C2*exp(b2*x);
      double dya = C1*b1*exp(b1*x) + C2*b2*exp(b2*x);

      double yn = ys[n_dim*i + 0];
      double dyn = ys[n_dim*i + 1];
	
      double yerr = log10(fabs(yn - ya));
      double dyerr = log10(fabs(dyn - dya));

      if (i % 20 == 0)
	{
	  printf("#%21s %21s %21s %21s %21s %21s %21s\n"
		 , "x"
		 , "y(x) (numeric)", "y'(x) (numeric)"
		 , "y(x) (analytic)", "y'(x) (analytic)"
		 , "y(x) logerror", "y'(x) logerror"
		 );
	}
      printf(" %21.15lf %21.15lf %21.15lf %21.15lf %21.15lf"
	     " %21.15lf %21.15lf\n"
	     , x
	     , yn, dyn
	     , ya, dya
	     , yerr, dyerr
	     );
    }
			    
  psede_ode_free(odesolver);
  return 0;
}
