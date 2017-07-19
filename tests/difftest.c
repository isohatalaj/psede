/**
 * @file difftest.c Example program demonstrating differentiation
 * methods supplied in psede.
 *  
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "psede.h"


void
fun_1d(double x, double *y, double *dy)
{
  const double w = 2*M_PI;
  const double a = 1.0;
  const double f = 0.0;

  const double c = cos(w*x + f);
  const double s = sin(w*x + f);
  const double u = 1 / (1 + a*x*x);
  
  *y = s * u;
  *dy = (w*c - 2*a*x*s*u) * u;
}

void
fun_3d(double x[3], double *fx, double dfx[3])
{
  const double s = 1.0;
  const double r2 = s*(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  const double wx = 2*M_PI;
  const double wy = 2*M_PI;
  const double wz = 2*M_PI;
  const double a = wx*x[0] + wy*x[1] + wz*x[2];

  *fx = sin(a) / (1 + r2);

  dfx[0] = cos(a)*wx / (1 + r2) - sin(a)*2*s*x[0] / ((1 + r2)*(1 + r2));
  dfx[1] = cos(a)*wy / (1 + r2) - sin(a)*2*s*x[1] / ((1 + r2)*(1 + r2));
  dfx[2] = cos(a)*wz / (1 + r2) - sin(a)*2*s*x[2] / ((1 + r2)*(1 + r2));
}

int
test_3d(int n[3], int out, double *max_err, double *elapsed)
{
  const int dim = 3;
  const int m = n[0]*n[1]*n[2];
  int status = 0;
  int i[3];

  double *x = NULL;
  double *y = NULL;
  double *dy_exact = NULL;
  double *dy_numer = NULL;

  psede_fct_t *fct = NULL;

  x = malloc(dim*m*sizeof(*x));
  y = malloc(m*sizeof(*y));
  dy_exact = malloc(dim*m*sizeof(*dy_exact));
  dy_numer = psede_fct_alloc_array(dim*m);

  if (x == NULL || y == NULL || dy_exact == NULL || dy_numer == NULL) 
    {
      status = 1;
      printf("# Out of memory, exiting...\n");
      goto exit;
    }

  fct = psede_fct_alloc();
  if (fct == NULL)
    {
      status = 1;
      printf("# Failed allocating Fast Chebyshev Transform work object\n");
      goto exit;
    }

  psede_Tx_nodes_multi_0(x, dim, n);

  for (i[0] = 0; i[0] < n[0]; ++i[0])
    {
      for (i[1] = 0; i[1] < n[1]; ++i[1])
  	{
  	  for (i[2] = 0; i[2] < n[2]; ++i[2])
  	    {
  	      const int ix = i[2] + n[2]*(i[1] + n[1]*i[0]);
	      
  	      fun_3d(&x[dim*ix], &y[ix], &dy_exact[dim*ix]);
  	      dy_numer[ix + 0*m] = y[ix];
  	      dy_numer[ix + 1*m] = y[ix];
  	      dy_numer[ix + 2*m] = y[ix];

	      /* printf(" %25.15le %25.15le %25.15le %25.15le\n", */
	      /* 	     x[dim*ix+0], x[dim*ix+1], x[dim*ix+2], y[ix]); */
  	    }
  	}
    }

  clock_t start = clock();

  int j;
  for (j = 0; j < dim; ++j)
    {
      status = psede_apply_multi(dy_numer + j*m, dim, n,
				 (psede_transform_t *) psede_Tx_diff_point_apply,
				 j,
				 (void *) fct);
      if (status)
	{
	  printf("# Failed differentiating\n");
	  goto exit;
	}
    }

  int n_fftelems = n[0]*n[1] + n[0]*n[2] + n[1]*n[2];

  clock_t end = clock();
  *elapsed = 1000.0 * (end - start) / ((double) CLOCKS_PER_SEC);
  *elapsed /= n_fftelems;

  int first = 1;

  for (i[0] = 0; i[0] < n[0]; ++i[0])
    {
      for (i[1] = 0; i[1] < n[1]; ++i[1])
  	{
  	  for (i[2] = 0; i[2] < n[2]; ++i[2])
  	    {
	      for (j = 0; j < dim; ++j)
		{
		  const int ix = i[2] + n[2]*(i[1] + n[1]*i[0]);
		  const double dy_ix_numer = dy_numer[ix + j*m];
		  const double dy_ix_exact = dy_exact[dim*ix + j];
		  const double diff = fabs(dy_ix_numer - dy_ix_exact);

		  if (first) { *max_err = diff; first = 0; }
		  else if (diff > *max_err) *max_err = diff;
		}
	    }
	}
    }
  
 exit:
  if (fct) psede_fct_free(fct);
  if (x) free(x);
  if (y) free(y);
  if (dy_exact) free(dy_exact);
  if (dy_numer) psede_fct_free_array(dy_numer);

  return status;
}



int
test_1d(int n, int out, double *max_err, double *elapsed)
{
  int status = 0;
  int i;

  double *x = NULL;
  double *y = NULL;
  double *dy_exact = NULL;
  double *dy_numer = NULL;

  psede_fct_t *fct = NULL;

  x = malloc(n*sizeof(*x));
  y = malloc(n*sizeof(*y));
  dy_exact = malloc(n*sizeof(*dy_exact));
  dy_numer = malloc(n*sizeof(*dy_numer));

  if (x == NULL || y == NULL || dy_exact == NULL || dy_numer == NULL) 
    {
      status = 1;
      printf("# Out of memory, exiting...\n");
      goto exit;
    }

  fct = psede_fct_alloc();
  if (fct == NULL)
    {
      status = 1;
      printf("# Failed allocating Fast Chebyshev Transform work object\n");
      goto exit;
    }

  psede_Tx_nodes_0(x, n);

  for (i = 0; i < n; ++i) 
    {
      fun_1d(x[i], &y[i], &dy_exact[i]);
      dy_numer[i] = y[i];
    }

  clock_t start = clock();
  status = psede_Tx_diff_point_apply(dy_numer, n, 1, 1, n, fct);
  clock_t end = clock();
  *elapsed = 1000.0 * (end - start) / ((double) CLOCKS_PER_SEC);
  *elapsed /= n;

  if (status) 
    {
      printf("# Failed differentiating\n");
      goto exit;
    }

  if (out)
    {
      printf("#%21s %21s %21s %21s %25s\n",
	     "x", "f(x)", "f'(x) exact", "f'(x) numer", "diff");
      for (i = 0; i < n; ++i)
	{
	  printf(" %21.15lf %21.15lf %21.15lf %21.15lf %25.15le\n",
		 x[i], y[i], dy_exact[i], dy_numer[i], dy_exact[i] - dy_numer[i]);
	}
    }

  *max_err = fabs(dy_exact[0] - dy_numer[0]);
  for (i = 1; i < n; ++i)
    {
      const double diff = fabs(dy_exact[i] - dy_numer[i]);
      if (diff > *max_err) *max_err = diff;
    }
  
 exit:
  if (fct) psede_fct_free(fct);
  if (x) free(x);
  if (y) free(y);
  if (dy_exact) free(dy_exact);
  if (dy_numer) free(dy_numer);

  return status;
}

int
main()
{
  int status = 0;
  int n;

  printf("# Testing 1D case...\n");

  printf("#%10s %25s %25s\n", "n", "log10 |err|", "time [ms]");
  for (n = 16; n <= 1 << 20; n = 3*n/2)
    {
      double max_err, elapsed;

      status = test_1d(n, 0, &max_err, &elapsed);
      
      printf(" %10d %25.15lf %25.15lf\n", n, log10(max_err), elapsed);
      if (status) goto exit;
    }

  printf("\n\n# Testing 3D case...\n");


  printf("#%5s %25s %25s\n", "n", "log10 |err|", "time [ms]");
  for (n = 8; n <= 128; n += 2)
    {
      double max_err, elapsed;
      int nn[3];
      nn[0] = n;
      nn[1] = n;
      nn[2] = n;
      /* int m = nn[0]*nn[1]*nn[2]; */

      status = test_3d(nn, 0, &max_err, &elapsed);
      
      printf(" %5d %25.15lf %25.15lf\n", n, log10(max_err), elapsed);
      if (status) goto exit;
    }
  

 exit:

  if (status) printf("# Failed!\n");
  else printf("# Success.\n");

  return status;
}
