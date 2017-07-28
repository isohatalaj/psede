
#include <math.h>
#include <stdlib.h>

#include "psede_fct.h"

int
psede_Tx_nodes(double *x, int n, int stride, int howmany, int dist, void *null)
{
  int i, j;
  const double z = M_PI/(n-1);

  for (i = 0; i < n; ++i)
    {
      x[i*stride] = cos(i*z);
    }

  for (j = 1; j < howmany; ++j)
    {
      for (i = 0; i < n; ++i)
	{
	  x[i*stride + j*dist] = x[i*stride];
	}
    }

  return 0;
}


void
psede_Tx_nodes_0(double *x, int n)
{
  psede_Tx_nodes(x, n, 1, 1, n, NULL);
}


void
psede_Tx_nodes_multi_0(double *x, int dims, const int *n)
{
  int dim;
  int stride = 1; /* stride in x,y,z,.. tuple units */

  int m = 1; /* size of x in x,y,z,... tuples */
  for (dim = 0; dim < dims; ++dim) m *= n[dim];

  for (dim = dims - 1; dim >= 0; --dim)
    {
      /* Do one set of points... */
      psede_Tx_nodes(x + dim, n[dim], dims*stride, 1, n[dim], NULL);

      /* ..and then duplicate it for the rest of the array. */
      int i, j = 0, k = 0;

      for (i = 0; i < m; ++i)
	{
	  x[dims*i + dim] = x[dims*stride*k + dim];

	  j++;
	  if (j == stride)
	    {
	      j = 0;
	      k++;
	      if (k == n[dim])
		{
		  k = 0;
		}
	    }
	}

      stride *= n[dim];
    }  
}


void
psede_identity(double *id, int size, int stride, int howmany, int dist)
{
  int i, j;
  
  for (i = 0; i < size; ++i)
    {
      for (j = 0; j < i; ++j)
	{
	  id[j*dist + i*stride] = 0.0;
	}

      id[j*dist + i*stride] = 1.0;

      for (j = i + 1; j < size; ++j)
	{
	  id[j*dist + i*stride] = 0.0;
	}
    }
}

static fftw_plan 
make_plan(double *x, 
	  int size, int stride,
	  int howmany, int dist)
{
  fftw_plan plan;
  int ns[1];
  fftw_r2r_kind kinds[1];
  int rank;
  unsigned flags;

  rank = 1;
  ns[0] = size;
  kinds[0] = FFTW_REDFT00;
  flags = FFTW_ESTIMATE;

  plan = fftw_plan_many_r2r(rank, ns, howmany,
			    x, NULL, stride, dist,
			    x, NULL, stride, dist,
			    kinds, flags);
  return plan;
}

static int
update_plan(psede_fct_t *self,
	    double *x, 
	    int size, int stride,
	    int howmany, int dist)
{
  int align = fftw_alignment_of(x);

  if (self->plan == NULL
      || self->size != size
      || self->stride != stride
      || self->howmany != howmany
      || self->dist != dist
      || self->align != align)
    {
      if (self->plan) fftw_destroy_plan(self->plan);
      self->plan = make_plan(x, size, stride, howmany, dist);
      if (self->plan == NULL) return 1;

      self->size = size;
      self->stride = stride;
      self->howmany = howmany;
      self->dist = dist;
      self->align = align;
    }

  return 0;
}

psede_fct_t *
psede_fct_alloc()
{
  psede_fct_t *self;

  self = malloc(sizeof(*self));
  if (self == NULL) return NULL;

  self->plan = NULL;

  return self;
}

void
psede_fct_free(psede_fct_t *self)
{
  if (self->plan) fftw_destroy_plan(self->plan);

  free(self);
}

double*
psede_fct_alloc_array(int size)
{
  return fftw_malloc(size*sizeof(double));
}

void
psede_fct_free_array(double *array)
{
  fftw_free(array);
}


int
psede_Tx_fct_apply(double *x,
		   int n,
		   int stride,
		   int howmany,
		   int dist,
		   psede_fct_t *self)
{
  int status;

  status = update_plan(self, x, n, stride, howmany, dist);
  if (status) return status;

  fftw_execute_r2r(self->plan, x, x);

  const double norm = 1.0/(n - 1);

  int i, j;
  for (j = 0; j < howmany; ++j)
    {
      x[0] *= 0.5*norm;
      for (i = 1; i < n - 1; ++i)
	{
	  x[i*stride] *= norm;
	} 
      x[(n-1)*stride] *= 0.5*norm;

      x += dist;
    }

  return 0;
}

int
psede_Tx_fct_apply_inv(double *x, 
		       int size, int stride,
		       int howmany, int dist,
		       psede_fct_t *self)
{
  int status;

  status = update_plan(self, x, size, stride, howmany, dist);
  if (status) return status;

  int i, j;
  for (j = 0; j < howmany; ++j)
    {
      double *xj = x + j*dist;

      for (i = 1; i < size - 1; ++i)
	{
	  xj[i*stride] *= 0.5;
	} 
    }

  fftw_execute_r2r(self->plan, x, x);

  return 0;
}

static int
psede_init_any_fct(psede_transf_t *tformer,
		   psede_transf_call_t *tform)
{
  psede_fct_t *fct = NULL;

  fct = psede_fct_alloc();
  if (fct == NULL) goto fail;

  tformer->transform = tform;
  tformer->finalize = (psede_transf_finalize_t*) psede_fct_free;
  tformer->params = fct;

  return 0;

 fail:

  if (fct) psede_fct_free(fct);

  return 1;
}

int
psede_init_Tx_fct(psede_transf_t *t, void *null)
{
  return psede_init_any_fct(t, (psede_transf_call_t*) psede_Tx_fct_apply);
}

int
psede_init_Tx_ifct(psede_transf_t *t, void *null)
{
  return psede_init_any_fct(t, (psede_transf_call_t*) psede_Tx_fct_apply_inv);
}

int
psede_init_Tx_nodes(psede_transf_t *t, void *null)
{
  t->transform = (psede_transf_call_t*) psede_Tx_nodes;
  t->finalize = NULL;
  t->params = NULL;

  return 0;
}

