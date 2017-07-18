
#include <stdlib.h>

#include "psede_util.h"

extern void
dgetrf_(int*, int*, double*, int*, int*, int*);

extern void
dgetrs_(char*, int*, int*, double*, int*, int*, double*, int*, int*);

extern void
dgemv_(char *trans,
       int *m, int *n,
       double *alpha, double *a, int *lda,
       double *x, int *incx,
       double *beta, double *y, int *incy);

extern void
dgemm_(char *transa, char *transb,
       int *m, int *n, int *k,
       double *alpha, double *a, int *lda, double *b, int *ldb,
       double *beta, double *c, int *ldc);

extern void
daxpy_(int *n, double *alpha, double *x, int *incx, double *y, int *incy);

extern void
dgels_(char *trans, int *n, int *m, int *nrhs,
       double *a, int *lda,
       double *b, int *ldb,
       double *work, int *lwork, int *info);


psede_linsolve_work_t *
psede_linsolve_work_alloc(int n)
{
  int status = 0;
  psede_linsolve_work_t *work = NULL;

  work = malloc(sizeof(*work));
  if (work == NULL) { status = 1; goto exit; }

  work->iwork = NULL;
  work->dgels_work = NULL;

  work->iwork_size = n;
  work->iwork = malloc(n*sizeof(*work->iwork));
  if (work->iwork == NULL) { status = 1; goto exit; }

  /* Query dgels for optimal work size. */
  /* TODO */

 exit:
  if (status)
    {
      psede_linsolve_work_free(work);
      return NULL;
    }

  return work;
}

void
psede_linsolve_work_free(psede_linsolve_work_t *work)
{
  if (work)
    {
      if (work->iwork) free(work->iwork);
      if (work->dgels_work) free(work->dgels_work);
      free(work);
    }
}

void
psede_transp_0(int n, double *m)
{
  int i, j;
  for (i = 0; i < n; ++i)
    {
      for (j = i + 1; j < n; ++j)
	{
	  const double temp = m[i + n*j];
	  m[i + n*j] = m[j + n*i];
	  m[j + n*i] = temp;
	}
    }
}


int
psede_linsolve_solve_0(int n, double *m, double *b,
		       psede_linsolve_work_t *work)
{
  int info;

  if (work->iwork_size < n) return 1;
  
  dgetrf_(&n, &n, m, &n, work->iwork, &info);
  if (info) return 1;

  int one = 1;

  dgetrs_("T", &n, &one, m, &n, work->iwork, b, &n, &info);
  if (info) return 1;

  return 0;
}

void
psede_gemv_0(int n, double alpha, double *a, double *x, double beta, double *y)
{
  int one = 1;
  dgemv_("T", &n, &n, &alpha, a, &n, x, &one, &beta, y, &one);
}

void
psede_gemm_0(int n, double alpha, double *a, double *b, double beta, double *c)
{
  dgemm_("T", "T", &n, &n, &n, &alpha, a, &n, b, &n, &beta, c, &n);
}

void
psede_daxpy_0(int n, double alpha, double *a, double *b)
{
  int one = 1;
  daxpy_(&n, &alpha, a, &one, b, &one);
}
