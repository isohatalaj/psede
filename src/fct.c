
#include <stdlib.h>

#include "psede_fct.h"

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
psede_fct_apply(psede_fct_t *self, 
		double *x,
		int n,
		int stride,
		int howmany,
		int dist)
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
psede_fct_apply_inv(psede_fct_t *self, 
		    double *x, 
		    int size, int stride,
		    int howmany, int dist)
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

