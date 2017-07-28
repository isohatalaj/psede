
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "psede.h"

#define DIM 3

int
myfunc(double *x, double *y, double *dy, double *a)
{
  const double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  const double fx = sin(a[0]*x[0]);
  const double fy = cos(a[1]*x[1]);
  const double fz = cos(a[2]*x[2]);

  *y = fx*fy*fz*exp(-r*r);

  dy[0] = a[0]*cos(a[0]*x[0])*fy*fz*exp(-r*r) - 2*x[0]*(*y);
  dy[1] = fx*(-1)*a[1]*sin(a[1]*x[1])*fz*exp(-r*r) - 2*x[1]*(*y);
  dy[2] = fx*fy*(-1)*a[2]*sin(a[2]*x[2])*exp(-r*r) - 2*x[2]*(*y);
    
  return 0;
}

int
myfunc_wrap(double *x, double *y, double *a)
{
  double dy[3];
  return myfunc(x, y, dy, a);
}

double
l2norm(int n, const double *x, const double *y)
{
  int i;
  double s = 0.0;
  for (i = 0; i < n; ++i)
    {
      const double t = x[i] - y[i];
      s += t*t;
    }

  return sqrt(s);
}

#define EXIT_WITH(msg) do { fprintf(stderr, msg); status = 1; goto exit; } while (0)

int
main()
{
  int i;
  int status = 0;
  const int ns[DIM] = {24, 25, 31};
  const int m = ns[0]*ns[1]*ns[2];
  const psede_colloc_t *collocs[DIM] = {&psede_Tx, &psede_Tx, &psede_Tx};

  double a[DIM] = {1.0, 2.0, 3.0};

  psede_multicolloc_t multicolloc;
  psede_multitransf_t nodes[DIM], diffs[DIM], func;
  double *xs[DIM], *fs, *ds[DIM];

  /* ******************************************************************************** */
  /* Pre-initialization */
  
  /* Recommended practice is to set all objects to their corresponding
   * nil values before initialization. */
  multicolloc = psede_multicolloc_nil;  
  for (i = 0; i < DIM; ++i)
    {
      xs[i] = ds[i] = NULL;
      nodes[i] = psede_multitransf_nil;
      diffs[i] = psede_multitransf_nil;
    }
  fs = NULL;
  func = psede_multitransf_nil;

  /* ******************************************************************************** */
  /* Initialization */ 

  for (i = 0; i < DIM; ++i)
    {
      xs[i] = psede_fct_alloc_array(m);
      ds[i] = psede_fct_alloc_array(m);
      
      if (xs[i] == NULL || ds[i] == NULL) EXIT_WITH("Out of memory.\n");
    }

  fs = psede_fct_alloc_array(m);
  if (fs == NULL) EXIT_WITH("Out of memory.\n");

  status = psede_init_multicolloc(&multicolloc, DIM, collocs,
				  PSEDE_RETAIN_OWNERSHIP);
  if (status) EXIT_WITH("Multicollocation object init failed.\n");


  status = psede_init_multifunc(&multicolloc, &func,
				(psede_multifunc_t*) myfunc_wrap, a);
  if (status) EXIT_WITH("Function tranform init failed.\n");


  for (i = 0; i < DIM; ++i)
    {
      status = psede_init_multinodes(&multicolloc, &nodes[i], i);
      if (status) EXIT_WITH("Failed initializing nodes transform.\n");
    }

  for (i = 0; i < DIM; ++i)
    {
      int orders[DIM] = {0};
      orders[i] = 1;
      
      status = psede_init_multidiff(&multicolloc, &diffs[i], orders);
      if (status) EXIT_WITH("Failed initializing diffs transform.\n");
    }

  /* ******************************************************************************** */
  /* Run transforms */

  /* Set input arrays to ones. Function application and nodes
   * operators multiply the inputs by their own output values to
   * produce the final result. */
  psede_multitransf_apply(&psede_multitransf_ones, xs[0], DIM, ns, 1, 1, 0);
  psede_multitransf_apply(&psede_multitransf_ones, xs[1], DIM, ns, 1, 1, 0);
  psede_multitransf_apply(&psede_multitransf_ones, xs[2], DIM, ns, 1, 1, 0);
  psede_multitransf_apply(&psede_multitransf_ones, fs, DIM, ns, 1, 1, 0);
  
  for (i = 0; i < DIM; ++i)
    {
      status = psede_multitransf_apply(&nodes[i], xs[i], DIM, ns, 1, 1, 0);
      if (status) EXIT_WITH("Failed applying nodes transform.\n");
    }

  status = psede_multitransf_apply(&func, fs, DIM, ns, 1, 1, 0);
  if (status) EXIT_WITH("Failed applying nodes transform.\n");

  for (i = 0; i < DIM; ++i) psede_copy(ds[i], fs, m, 1, 1, 0);

  for (i = 0; i < DIM; ++i)
    {
      status = psede_multitransf_apply(&diffs[i], ds[i], DIM, ns, 1, 1, 0);
      if (status) EXIT_WITH("Failed applying diff transform.\n");
    }

  /* ******************************************************************************** */
  /* Output */

  double max_err = 0.0;
  for (i = 0; i < m; ++i)
    {
      double y_x, dy_x[3];
      double x_x[3];
      double err;

      x_x[0] = xs[0][i];
      x_x[1] = xs[1][i];
      x_x[2] = xs[2][i];

      myfunc(x_x, &y_x, dy_x, a);

      if (y_x != fs[i]) EXIT_WITH("Function in transform array does not match "
				  "the function value\n");

      double dy_n[3];
      dy_n[0] = ds[0][i];
      dy_n[1] = ds[1][i];
      dy_n[2] = ds[2][i];
      
      err = l2norm(DIM, dy_x, dy_n);
      if (err > max_err) max_err = err;

      printf(" %25.15lf %25.15lf %25.15lf %25.15lf %25.15lf"
	     " %25.15lf %25.15lf %25.15lf %25.15lf %25.15lf\n",
	     xs[0][i], xs[1][i], xs[2][i], fs[i],
	     ds[0][i], ds[1][i], ds[2][i],
	     dy_x[0], dy_x[1], dy_x[2]);
    }

  fprintf(stderr, "# log10 ||err|| = %lf\n", log10(max_err));

  
  /* ******************************************************************************** */
  /* Cleanup and exit */
  
 exit:

  if (fs) psede_fct_free_array(fs);  
  for (i = 1; i < DIM; ++i) if (xs[i]) psede_fct_free_array(xs[i]);
  for (i = 1; i < DIM; ++i) if (ds[i]) psede_fct_free_array(ds[i]);
  
  for (i = 1; i < DIM; ++i) psede_multitransf_destroy(&nodes[i]);
  for (i = 1; i < DIM; ++i) psede_multitransf_destroy(&diffs[i]);

  psede_multitransf_destroy(&func);
  psede_multicolloc_destroy(&multicolloc);

  return status;
}
