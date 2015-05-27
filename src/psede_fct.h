
#ifndef PSEDE_FCT_H
#define PSEDE_FCT_H

#include <fftw3.h>

/**
 * Fast Chebyshev transform workspace
 */
typedef struct {
  int size;
  int stride;
  int howmany;
  int dist;
  int align;
  fftw_plan plan;
} psede_fct_t;

/**
 * Allocate a Fast Chebyshev Transform object
 */
psede_fct_t*
psede_fct_alloc();

/**
 * Free a Fast Chebyshev Transform object
 */ 
void
psede_fct_free(psede_fct_t *self);

/**
 * Allocate memory best suited for a Fast Chebyshev Transform, in
 * units of sizeof(double). Presently this is just a wrapper for
 * FFTW's aligned allocation functions, given to present a dependency
 * independent API (ie. to accommodate different FFT backends).
 * Arrays allocated via this function must be freed using
 * psede_fct_free_array.
 */
double*
psede_fct_alloc_array(int size);

/**
 * Free an array alloceted using psede_fct_alloc_array
 */
void
psede_fct_free_array(double *array);

/**
 * Apply a Fast Chebyshev Transform to a vector or vectors; converts
 * point values to mode coefficients. Return zero for success.
 */
int
psede_fct_apply(psede_fct_t *self,
		double *x, int size, int stride,
		int howmany, int dist);

/**
 * Apply an Inverse Fast Chebyshev Transform to a vector or vectors;
 * converts mode coefficients to point values. Return zero for
 * success.
 */
int
psede_fct_apply_inv(psede_fct_t *self,
		    double *x, int size, int stride,
		    int howmany, int dist);

#endif /* PSEDE_FCT_H */
