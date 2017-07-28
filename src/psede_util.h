
#ifndef PSEDE_UTIL_H
#define PSEDE_UTIL_H

#include <stdio.h>

#define Psede_debug(str, ...) do { fprintf(stderr, "# debug %s, %s (%d): " str, __FUNCTION__, __FILE__, __LINE__, ##__VA_ARGS__); fflush(stderr); } while(0)

/**
 * Enum used when passing pointers to other objects. The value
 * `PSEDE_PASS_OWNERSHIP` indicates that the responsibility for
 * finalising a pointer argument is transferred to the routine being
 * called. The opposite case where the caller retains its
 * responsibility for the pointer is flagged by
 * `PSEDE_RETAIN_OWNERSHIP`.
 */
typedef enum {
  PSEDE_PASS_OWNERSHIP,
  PSEDE_RETAIN_OWNERSHIP
} psede_ownership_t;

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

typedef void (psede_finalize_t)(void *);


/**
 * Linked list for storing pointers that need to be finalized.  Used
 * by objects that allocate dynamically varying numbers of other
 * objects.
 */
/* typedef struct psede_cleanup_list_t; */
/* typedef struct { */
/*   psede_finalize_t *finalize; */
/*   void *params; */
/*   psede_cleanup_list_t *next; */
/* } psede_cleanup_list_t; */

/* int */
/* psede_cleanup_list_push(psede_cleanup_list_t **list, */
/* 			psede_finalize_t *finalize, void *params); */

/* void */
/* psede_cleanup_list_finalize(psede_cleanup_list_t *list); */


/**
 * Copy multiple strided arrays.
 */
void
psede_copy(double *dest, const double *source,
	   int n, int stride, int howmany, int dist);

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
