
#include <stdlib.h>

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

/* ******************************************************************************** */

/**
 * Transform parameters for one-dimensional differentiation, Chebyshev
 * extrema in point space.
 */
typedef struct {
  int order;
  psede_transf_t fct;
  psede_transf_t ifct;
} psede_Tx_diff_params_t;

psede_Tx_diff_params_t *
psede_Tx_diff_params_alloc();

void
psede_Tx_diff_params_free(psede_Tx_diff_params_t*);

int
psede_Tx_diff_apply(double *target, int size, int stride,
		    int howmany, int dist,
		    psede_Tx_diff_params_t *params);



psede_Tx_diff_params_t *
psede_Tx_diff_params_alloc(int order)
{
  int status;
  psede_Tx_diff_params_t *params;

  params = malloc(sizeof(*params));
  if (params == NULL) return NULL;
  
  status = psede_init_Tx_fct(&params->fct, NULL);
  if (status) goto fail;

  status = psede_init_Tx_ifct(&params->ifct, NULL);
  if (status) goto fail;

  params->order = order;

  return params;

 fail:
  psede_Tx_diff_params_free(params);

  return NULL;
}

void
psede_Tx_diff_params_free(psede_Tx_diff_params_t *params)
{
  if (params)
    {
      psede_transf_destroy(&params->fct);
      psede_transf_destroy(&params->ifct);

      free(params);
    }
}

int
psede_Tx_diff_apply(double *target, int size, int stride,
		    int howmany, int dist,
		    psede_Tx_diff_params_t *params)
{
  int status;
  int i;
  
  if (params->order == 0) return 0;

  status = psede_transf_apply(&params->fct,
				   target, size,
				   stride, howmany, dist);
  if (status) return status;

  if (params->order > 0)
    {
      for (i = 0; i < params->order; ++i)
	{
	  status = psede_Tx_diff_mode_apply(target, size,
					    stride, howmany, dist,
					    NULL);
	  if (status) return status;
	}
    }
  else
    {
      for (i = 0; i < -params->order; ++i)
	{
	  status = psede_Tx_integ_mode_apply(target, size,
					     stride, howmany, dist,
					     NULL);
	  if (status) return status;
	}
    }

  status = psede_transf_apply(&params->ifct,
				   target, size,
				   stride, howmany, dist);
  if (status) return status;

  return 0;
}

int
psede_init_Tx_diff(psede_transf_t *transf, int order, void *null)
{
  psede_Tx_diff_params_t *params;

  params = psede_Tx_diff_params_alloc(order);
  if (params == NULL) return 1;

  transf->transform = (psede_transf_call_t*) psede_Tx_diff_apply;
  transf->finalize = (psede_transf_finalize_t*) psede_Tx_diff_params_free;
  transf->params = params;

  return 0;
}
