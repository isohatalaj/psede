
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "psede.h"

int
myfunc(double x, double *y, double *dy, double a)
{
  const double omega = 2*M_PI*a;
  *y = sin(omega*x);
  *dy = omega*cos(omega*x);
  return 0;
}

int
myfunc_wrap(double x, double *y, double *a)
{
  int status;
  double dy;
  status = myfunc(x, y, &dy, *a);
  return status;
}

int
main()
{
  const int n = 32;
  const psede_colloc_t *colloc = &psede_Tx;

  psede_transf_t
    nodes = psede_tnil,
    diff = psede_tnil,
    func = psede_tnil;

  double *xs = NULL, *fs = NULL, *ds = NULL;
  double a = 2.0;

  /*
   * Allocate arrays for the x coordinates (nodes), the function
   * values (our function f evaluated at the nodes), and the
   * derivative values (f' at nodes).
   */
  xs = psede_fct_alloc_array(n);  
  fs = psede_fct_alloc_array(n);
  ds = psede_fct_alloc_array(n); 

  /*
   * Construct the needed transformations: nodes sets an array to the
   * node coordinates, diff affects a differentiation on input of
   * function values, and func is a function multiplication (takes an
   * input vector and multiplies each element with the value of the
   * function evaluated at the corresponding node).
   */
  psede_init_nodes(colloc, &nodes);
  psede_init_diff(colloc, &diff, 1); 
  psede_init_func(colloc, &func, (psede_func_t*) myfunc_wrap, &a);

  /*
   * Apply transformations. First set the array to the node
   * coordinates, then prepare the function values vector by setting
   * it to ones (using psede_transf_ones transform). The function
   * values at nodes are obtained by using the function multiplication
   * transform on the vector of ones. We make a separate copy of the
   * function values into the array ds, upon which we apply the
   * differentiation operator to get the approximate function
   * derivatives.
   */
  psede_transf_apply_0(&nodes, xs, n);
  psede_transf_apply_0(&psede_transf_ones, fs, n);
  psede_transf_apply_0(&func, fs, n);
  psede_copy(ds, fs, n, 1, 1, 0);
  psede_transf_apply_0(&diff, ds, n);

  /*
   * Print results and compute the max norm error.
   */
  double max_err = 0.0;  

  printf("#%25s %25s %25s %25s\n",
	 "x", "f(x)", "f'(x) num", "f'(x) exact");
  int i;
  for (i = 0; i < n; ++i)
    {
      double y, dy, err;
      myfunc(xs[i], &y, &dy, a);
      err = fabs(ds[i] - dy);
      if (err > max_err) max_err = err;
      
      printf(" %25.15lf %25.15lf %25.15lf %25.15lf\n",
	     xs[i], fs[i], ds[i], dy);
    }

  printf("\n# max. error = %le\n", max_err);

  /*
   * Clean up. For each transform we have initialized, there should be
   * a matching destroy call (note that some routines, not used here,
   * have the option of transferring the cleanup responsibility to the
   * called routined). Also all arrays should be freed.
   */

  psede_transf_destroy(&func);
  psede_transf_destroy(&nodes);
  psede_transf_destroy(&diff);

  psede_fct_free_array(ds);
  psede_fct_free_array(fs);
  psede_fct_free_array(xs);
  
  return 0;
}




