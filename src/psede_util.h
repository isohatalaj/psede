
#ifndef PSEDE_UTIL_H
#define PSEDE_UTIL_H

typedef struct {
  int iwork_size;
  int *iwork;

  int dgels_lwork;
  double *dgels_work;
} psede_linsolve_work_t;

psede_linsolve_work_t *
psede_linsolve_work_alloc(int n);

void
psede_linsolve_work_free(psede_linsolve_work_t *work);

/**
 * In-place transpose for `n`-by-`n` matrix.
 */
void
psede_transp_0(int n, double *m);

/**
 * Solve linear equation m.x = b, where `m` is `n`-by-`n` matrix and
 * `b` is a vector of `n` elements. A pre-allocated workspace is to be
 * provided in `work`. The input matrix is destroyed in the process.
 */
int
psede_linsolve_solve_0(int n, double *m, double *b,
		       psede_linsolve_work_t *work);


/**
 * Matrix-vector produxt, computes `y <- alpha*A*x + beta*y. 
 * Simply BLAS dgemv fpr dense-packed (size == lead dimension)
 * square-row major matrices.
 */
void
psede_gemv_0(int n, double alpha, double *a, double *x, double beta, double *y);

/**
 * Matrix-matrix product, computes `c <- alpha*a.b + beta*c`. Simply
 * BLAS dgemm for dense-packed (size == lead dimension) square
 * row-major matrices.
 */
void
psede_gemm_0(int n, double alpha, double *a, double *b, double beta, double *c);

/**
 * Array scaled add, computes `b <- alpha a + b`. Simply BLAS daxpy
 * with stride one arrays.
 */
void
psede_daxpy_0(int n, double alpha, double *a, double *b);

/**
 * Initialize a set of vectors to sequence of unit vectors, ie.
 * construct a given number of columns/rows of an identity matrix.
 */
void
psede_identity(double *matrix, int size, int stride, int howmany, int dist);



#endif
