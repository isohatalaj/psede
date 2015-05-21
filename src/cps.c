
#include <stdio.h>
#include <math.h>

#include <lapacke.h>

#include "cps.h"
#include "cps_ode.h"


/** Set C = alpha A. Small matrices assumed. */
void
cps_msclset(int n,
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

/** Set C = alpha A + beta B. Small matrices assumed. */
void
cps_mscladd(int n,
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


void
cps_nodes(double *x, int n, int stride)
{
  int i;
  const double z = M_PI/(n-1);
  for (i = 0; i < n; ++i)
    {
      x[i*stride] = cos(i*z);
    }
}

void
cps_unit(int i, double *u, int n, int stride)
{
  int j;

  for (j = 0; j < i; ++j)
    {
      u[j*stride] = 0.0;
    }

  u[i*stride] = 1.0;

  for (j = i + 1; j < n; ++j)
    {
      u[j*stride] = 0.0;
    }
  
}

void
cps_id(double *id, int n, int lda)
{
  int i, j;

  for (i = 0; i < n; ++i)
    {
      for (j = 0; j < i; ++j)
	{
	  id[i*lda+j] = 0.0;
	}
      
      id[i*lda+i] = 1.0;

      for (++j; j < n; ++j)
	{
	  id[i*lda+j] = 0.0;
	}
    }
}

void
cps_apply_D_mod(double *x, int n, int stride)
{
  int i, j;
  double s;

  for (i = 0; i < n; ++i)
    {
      s = 0.0;

      if (i)
	{
	  for (j = i + 1; j < n - 1; j += 2)
	    {
	      s += 2*j*x[j*stride];
	    }
	  if (j == n - 1) s += 2*j*x[j*stride];
	}
      else
	{
	  for (j = i + 1; j < n - 1; j += 2)
	    {
	      s += j*x[j*stride];
	    }
	  if (j == n - 1) s += j*x[j*stride];
	}
      
      x[i*stride] = s;
    }
}

fftw_plan 
cps_plan_T(double *x, int n, int stride)
{
  fftw_plan plan;
  int ns[1];
  fftw_r2r_kind kinds[1];
  int rank, howmany;
  int *inembed, *onembed;
  int idist, odist;
  unsigned flags;

  howmany = 1;
  rank = 1;
  ns[0] = n;
  kinds[0] = FFTW_REDFT00;
  inembed = onembed = NULL;
  idist = odist = 0;
  flags = FFTW_ESTIMATE;

  plan = fftw_plan_many_r2r(rank, ns, howmany,
			    x, inembed, stride, idist,
			    x, onembed, stride, odist,
			    kinds, flags);

  return plan;
}

void
cps_apply_T(double *x, int n, int stride, const fftw_plan plan)
{
  int i;
  const double norm = 1.0/(n - 1);

  fftw_execute(plan);

  x[0] *= 0.5*norm;
  for (i = 1; i < n - 1; ++i)
    {
      x[i*stride] *= norm;
    }
  x[(n-1)*stride] *= 0.5*norm;
}


void
cps_apply_T_inv(double *x, int n, int stride, const fftw_plan plan)
{
  int i;

  for (i = 1; i < n - 1; ++i)
    {
      x[i*stride] *= 0.5;
    }

  fftw_execute(plan);
}
