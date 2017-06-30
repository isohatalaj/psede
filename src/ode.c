
#include <stdlib.h>

#include "psede_ode.h"

extern void
dgetrf_(int*, int*, double*, int*, int*, int*);

extern void
dgetrs_(char*, int*, int*, double*, int*, int*, double*, int*, int*);


static inline int
get_work_size(int size, int dim)
{
  int cap = size*dim;
  return cap*cap + 2*dim*dim + dim;
}

static inline int
get_iwork_size(int size, int dim)
{
  int cap = size*dim;
  return cap;
}

static void
zero_self(psede_ode_t *self)
{
  self->size = 0;
  self->dim = 0;
  self->x = NULL;
  self->y = NULL;
  self->D = NULL;

  self->fct = NULL;
  
  self->work_size = 0;
  self->iwork_size = 0;
  self->work = NULL;
  self->iwork = NULL;
}

psede_ode_t*
psede_ode_alloc(int size, int dim)
{
  int status;
  psede_ode_t *self;

  self = malloc(sizeof(*self));
  if (self == NULL) return NULL;

  zero_self(self);

  int cap = size*dim;

  self->size = size;
  self->dim = dim;

  self->x = psede_fct_alloc_array(size);
  if (self->x == NULL) goto fail;

  self->y = psede_fct_alloc_array(cap);
  if (self->y == NULL) goto fail;

  self->D = psede_fct_alloc_array(size*size);
  if (self->D == NULL) goto fail;

  self->work_size = get_work_size(size, dim);
  self->work = psede_fct_alloc_array(self->work_size);
  if (self->work == NULL) goto fail;

  self->iwork_size = get_iwork_size(size, dim);
  self->iwork = malloc(self->iwork_size*sizeof(*self->iwork));
  if (self->iwork == NULL) goto fail;

  self->fct = psede_fct_alloc();
  if (self->fct == NULL) goto fail;

  psede_nodes(self->x, size, 1);

  status = psede_diff_point_matrix(self->D, size, 1, size, size,
				   self->fct);
  if (status) goto fail;

  return self;

 fail:
  psede_ode_free(self);

  return NULL;
}

void
psede_ode_free(psede_ode_t *self)
{
  if (self->x) psede_fct_free_array(self->x);
  if (self->y) psede_fct_free_array(self->y);
  if (self->D) psede_fct_free_array(self->D);
  if (self->fct) psede_fct_free(self->fct);
  if (self->work) psede_fct_free_array(self->work);
  if (self->iwork) free(self->iwork);

  free(self);
}

static void
msclset(int n,
	const double *A, int lda, double alpha,
	double *C, int ldc)
{
  int i, j;
  for (i = 0; i < n; ++i)
    {
      for (j = 0; j < n; ++j)
	{
	  C[i*ldc+j] = alpha*A[i*lda+j];
	}
    }
}

static void
mscladd(int n,
	const double *A, int lda, double alpha,
	const double *B, int ldb, double beta,
	double *C, int ldc)
{
  int i, j;
  for (i = 0; i < n; ++i)
    {
      for (j = 0; j < n; ++j)
	{
	  C[i*ldc+j] = alpha*A[i*lda+j] + beta*B[i*ldb+j];
	}
    }
}

int
psede_ode_linear_solve(psede_ode_t *self, 
		       double x0, double x1,
		       psede_ode_linear_fun_t fun, void *params,
		       int nbcs,
		       const double *B, 
		       const double *gamma,
		       const int *ls)
{
  int status;
  int i, j;
  int n = self->size;
  int k = self->dim;
  int nk = n * k;
  int kk = k * k;

  double s1 = 0.5*(x1-x0);
  double s0 = 0.5*(x0+x1);

  /* Setup finite dimensional representation of the ODE using the CGL
   * pseudospectral differentiation operator.  */

  double *M = self->work;
  double *A = M + nk*nk;
  double *Q = A + kk;
  double *C = Q + kk;

  for (j = 0; j < n; ++j)
    {
      double x = s1*self->x[j]+s0;

      status = fun(x, A, Q, C, params);
      if (status) return status;

      for (i = 0; i < n; ++i)
	{
	  double Dij = self->D[n*j+i]/s1;
	  
	  if (i == j)
	    {
	      mscladd(k, 
		      A, k, -1.0,
		      Q, k, Dij,
		      M + i*k + n*kk*j, nk);
	    }
	  else
	    {
	      msclset(k,
		      Q, k, Dij,
		      M + i*k + n*kk*j, nk);
	    }
	}

      for (i = 0; i < k; ++i)
	{
	  self->y[k*j + i] = C[i];
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
      	  M[l*nk+i] = B[2*k*j+i+k];
      	}
      for (; i < nk - k; ++i) M[l*nk+i] = 0.0;
      for (i = 0; i < k; ++i)
      	{
      	  M[l*nk + nk - k + i] = B[2*k*j+i];
      	}
      
      /* Update rhs vector: */
      self->y[l] = gamma[j];
    }

  int info;
  
  dgetrf_(&nk, &nk, M, &nk, self->iwork, &info);
  if (info)
    {
      fprintf(stderr, "# error: dgetrf yielded bad info %d\n", info);
      return info;
    }

  int one = 1;

  dgetrs_("T", &nk, &one, M, &nk, self->iwork, self->y, &nk, &info);
  if (info)
    {
      fprintf(stderr, "# error: dgetrs yielded bad info %d\n", info);
      return info;
    }

  return 0;
}
