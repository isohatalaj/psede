
#include <lapacke.h>

#include "cps_ode.h"

cps_ode_t *
cps_ode_make(int n_size,
	     int n_dim)
{
  int n_sys_size = n_size * n_dim;
  int n_size_sqr = n_size * n_size;
  int n_dim_sqr = n_dim * n_dim;

  cps_ode_t *self = malloc(sizeof(*self));
  if (self == NULL) return NULL;

  self->n_size = n_size;
  self->n_dim = n_dim;

  self->x = NULL;
  self->y = NULL;
  self->M = NULL;
  self->D = NULL;
  self->F = NULL;
  self->A = NULL;
  self->Q = NULL;
  self->ipiv = NULL;
  self->ply = NULL;

  self->x = malloc(n_size*sizeof(*self->x));
  if (!self->x) goto fail;

  self->y = fftw_malloc(n_sys_size*sizeof(*self->y));
  if (!self->y) goto fail;

  self->M = malloc(n_sys_size*n_sys_size*sizeof(*self->M));
  if (!self->M) goto fail;

  self->D = malloc(n_size_sqr*sizeof(*self->D));
  if (!self->D) goto fail;

  self->F = malloc(n_dim*sizeof(*self->F));
  if (!self->F) goto fail;

  self->A = malloc(n_dim_sqr*sizeof(*self->A));
  if (!self->A) goto fail;

  self->Q = malloc(n_dim_sqr*sizeof(*self->Q));
  if (!self->Q) goto fail;

  self->ipiv = malloc(n_sys_size*sizeof(*self->ipiv));
  if (!self->ipiv) goto fail;

  self->ply = cps_plan_T(self->y, n_size, 1);
  if (!self->ply) goto fail;

  cps_nodes(self->x, n_size, 1);

  int i, j;
  for (i = 0; i < n_size; ++i)
    {
      cps_unit(i, self->y, n_size, 1);
      cps_apply_T(self->y, n_size, 1, self->ply);
      cps_apply_D_mod(self->y, n_size, 1);
      cps_apply_T_inv(self->y, n_size, 1, self->ply);

      for (j = 0; j < n_size; ++j)
	{
	  self->D[n_size*j+i] = self->y[j];
	}
    }

  return self;


 fail:

  cps_ode_free(self);
  return NULL;
}

void
cps_ode_free(cps_ode_t *self)
{
  if (self->ply) fftw_destroy_plan(self->ply);
  if (self->ipiv) free(self->ipiv);
  if (self->F) free(self->F);
  if (self->A) free(self->A);
  if (self->Q) free(self->Q);
  if (self->D) free(self->D);
  if (self->M) free(self->M);
  if (self->y) fftw_free(self->y);
  if (self->x) free(self->x);
  free(self);
}

int
cps_linode_solve(cps_ode_t *self, 
		 double x0, double x1,
		 cps_linode_fun_t fun, void *params,
		 int nbcs,
		 const double *bcs, 
		 const double *gammas,
		 const int *ls)
{
  int status;
  int i, j;
  int n = self->n_size;
  int k = self->n_dim;
  int nk = n * k;
  int kk = k * k;

  double s1 = 0.5*(x1-x0);
  double s0 = 0.5*(x0+x1);

  /* Setup finite dimensional representation of the ODE using the CGL
   * pseudospectral differentiation operator.  */

  for (j = 0; j < n; ++j)
    {
      status = fun(s1*self->x[j]+s0, self->A, self->Q, self->F, params);
      if (status) return status;

      for (i = 0; i < n; ++i)
	{
	  double Dij = self->D[n*j+i]/s1;
	  
	  if (i == j)
	    {
	      cps_mscladd(k,
			  self->A, k, -1.0,
			  self->Q, k, Dij,
			  &self->M[i*k + n*kk*j], nk);
	    }
	  else
	    {
	      cps_msclset(k,
			  self->Q, k, Dij,
			  &self->M[i*k + n*kk*j], nk);
	    }
	}

      for (i = 0; i < k; ++i)
	{
	  self->y[k*j + i] = self->F[i];
	}
    }

  for (j = 0; j < nbcs; ++j)
    {
      /* Determine row number l onto which this boundary condition is
       * injected. */
      int l;

      if (ls[j] > 0) l = nk - (k - ls[j]) - 1;
      else l = -ls[j] - 1;

      /* Update the M matrix: */
      for (i = 0; i < k; ++i)
      	{
      	  self->M[l*nk+i] = bcs[2*k*j+i+k];
      	}
      for (; i < nk - k; ++i) self->M[l*nk+i] = 0.0;
      for (i = 0; i < k; ++i)
      	{
      	  self->M[l*nk + nk - k + i] = bcs[2*k*j+i];
      	}
      
      /* Update rhs vector: */
      self->y[l] = gammas[j];
    }

  int info;
  
  info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 
			nk, nk, self->M, nk, self->ipiv);
  if (info)
    {
      fprintf(stderr, "# error: dgetrf yielded bad info %d\n", info);
      return info;
    }
  
  info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', 
			nk, 1, self->M, nk, self->ipiv, self->y, 1);
  if (info)
    {
      fprintf(stderr, "# error: dgetrs yielded bad info %d\n", info);
      return info;
    }

  return 0;
}
