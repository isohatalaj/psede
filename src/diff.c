
#include "psede_diff.h"

int
psede_Tx_diff_mode_apply(double *x, int n, int stride, 
			 int howmany, int dist, void *params)
{
  int i, k;
  double xi, xii;

  for (k = 0; k < howmany; ++k)
    {
      xi = x[(n-1)*stride];
      x[(n-1)*stride] = 0.0;
      
      xii = x[(n-2)*stride];
      x[(n-2)*stride] = 2.0*(n-1)*xi;
      xi = xii;

      for (i = n - 3; i > 0; --i)
	{
	  xii = x[i*stride];
	  x[i*stride] = x[(i+2)*stride] + 2.0*(i+1)*xi;
	  xi = xii;
	}

      x[0] = 0.5*x[2*stride] + xi;

      x += dist;
    }

  return 0;
}

int
psede_Tx_integ_mode_apply(double *x, int n, int stride, 
			  int howmany, int dist, void *params)
{
  int i, k;
  double xprev, xtemp;

  for (k = 0; k < howmany; ++k)
    {
      double s = 0.0;
      double p = n % 2 == 0 ? 1.0 : -1.0;

      xtemp = x[(n - 1)*stride];
      x[(n - 1)*stride] = 0.5*x[(n - 2)*stride] / (n - 1);
      s += p * xtemp / ((n - 1)*(n - 1) - 1);
      p = -p;
      xprev = xtemp;

      for (i = n - 2; i >= 2; --i)
	{
	  xtemp = x[i*stride];
	  x[i*stride] = 0.5*(x[(i-1)*stride] - xprev) / i;
	  s += p * xtemp / (i*i - 1);
	  p = -p;
	  xprev = xtemp;
	}

      xtemp = x[1*stride];
      x[1*stride] = x[0*stride] - 0.5*xprev;
      xprev = xtemp;
      
      x[0*stride] = x[0*stride] - 0.25*xprev + s;
      
      x += dist;
    }

  return 0;
}

int
psede_Tx_diff_point_apply(double *x, int size, int stride,
			  int howmany, int dist, psede_fct_t *fct)
{
  int status;

  status = psede_Tx_fct_apply(x, size, stride, howmany, dist, fct);
  if (status) return status;

  status = psede_Tx_diff_mode_apply(x, size, stride, howmany, dist, NULL);
  if (status) return status;
  
  status = psede_Tx_fct_apply_inv(x, size, stride, howmany, dist, fct);
  if (status) return status;

  return 0;
}

int
psede_Tx_integ_point_apply(double *x, int size, int stride,
			   int howmany, int dist, psede_fct_t *fct)
{
  int status;

  status = psede_Tx_fct_apply(x, size, stride, howmany, dist, fct);
  if (status) return status;

  status = psede_Tx_integ_mode_apply(x, size, stride, howmany, dist, NULL);
  if (status) return status;

  status = psede_Tx_fct_apply_inv(x, size, stride, howmany, dist, fct);
  if (status) return status;

  return 0;
}

/* int */
/* psede_Tx_diff_point_apply_multi_0(double *x, int diff_dim, */
/* 				  int dims, const int *sizes, psede_fct_t *fct) */
/* { */
/*   int status; */
/*   int n_below = 1, n_above = 1; */

/*   int i; */
/*   for (i = 0; i < dims; ++i) */
/*     { */
/*       if (i < diff_dim) n_above *= sizes[i]; */
/*       else if (i > diff_dim) n_below *= sizes[i]; */
/*     } */

/*   if (n_below == 1) */
/*     { */
/*       for (i = 0; i < n_below; ++i) */
/* 	{ */
/* 	  status = psede_Tx_diff_point_apply(x + i, sizes[diff_dim], n_below, */
/* 					     n_above, sizes[diff_dim]*n_below, fct); */
/* 	  if (status) return status; */
/* 	} */
/*     } */
/*   else */
/*     { */
/*       for (i = 0; i < n_above; ++i) */
/*   	{ */
/*   	  status = psede_Tx_diff_point_apply(x + i*sizes[diff_dim]*n_below, */
/*   					     sizes[diff_dim], n_below, n_below, 1, fct); */
/*   	  if (status) return status; */
/*   	} */
/*     } */

/*   return status; */
/* } */

int
psede_Tx_integ_point_apply_multi_0(double *x, int diff_dim,
				   int dims, const int *sizes, psede_fct_t *fct)
{
  int status;
  int n_below = 1, n_above = 1;

  int i;
  for (i = 0; i < dims; ++i)
    {
      if (i < diff_dim) n_above *= sizes[i];
      else if (i > diff_dim) n_below *= sizes[i];
    }

  if (n_below == 1)
    {
      for (i = 0; i < n_below; ++i)
	{
	  status = psede_Tx_integ_point_apply(x + i, sizes[diff_dim], n_below,
					      n_above, sizes[diff_dim]*n_below, fct);
	  if (status) return status;
	}
    }
  else
    {
      for (i = 0; i < n_above; ++i)
  	{
  	  status = psede_Tx_integ_point_apply(x + i*sizes[diff_dim]*n_below,
					      sizes[diff_dim], n_below, n_below, 1, fct);
  	  if (status) return status;
  	}
    }

  return status;
}

/* int */
/* psede_Tx_diff_point_matrix_multi_0(double *diff_point, */
/* 				   int diff_dim, */
/* 				   int dims, const int *sizes, */
/* 				   int init, */
/* 				   psede_fct_t *fct) */
/* { */
/*   int status; */
/*   int i, m = 1; */
/*   for (i = 0; i < dims; ++i) m *= sizes[i]; */

/*   if (init) psede_identity(diff_point, m, 1, m, m); */

/*   /\* TODO: Make these transpositions unnecessary *\/ */
/*   psede_transp_0(m, diff_point); */
  
/*   for (i = 0; i < m; ++i) */
/*     { */
/*       status = psede_apply_multi(diff_point + i*m, dims, sizes, */
/* 				 (psede_transform_t *) psede_Tx_diff_point_apply, */
/* 				 diff_dim, */
/* 				 (void *) fct); */
/*       if (status) return status; */
/*     } */

/*   psede_transp_0(m, diff_point); */
  
/*   return 0; */
/* } */

int
psede_Tx_integ_point_matrix_multi_0(double *diff_point,
				    int diff_dim,
				    int dims, const int *sizes,
				    int init,
				    psede_fct_t *fct)
{
  int status;
  int i, m = 1;
  for (i = 0; i < dims; ++i) m *= sizes[i];

  if (init) psede_identity(diff_point, m, 1, m, m);

  /* TODO: Make these transpositions unnecessary */
  psede_transp_0(m, diff_point);

  for (i = 0; i < m; ++i)
    {
      status = psede_Tx_integ_point_apply_multi_0(diff_point + i*m,
						  diff_dim, dims, sizes, fct);
      if (status) return status;
    }

  psede_transp_0(m, diff_point);

  return 0;
}
