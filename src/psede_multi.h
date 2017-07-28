
#ifndef PSEDE_MULTI_H
#define PSEDE_MULTI_H

#include "psede_util.h"


/* ******************************************************************************** */

/* To provide flexible/expandable API, the operations such as
 * pseudospectral differentiation, integration, multiplication by a
 * function, etc. are abstracted into "transforms". These can
 * subsequently manipulated, combined, etc. to create representations
 * for more general operators, that e.g. map to the left-hand sides of
 * equations we wish to solve. The transforms can be subsequently
 * plugged into e.g. matrix free solvers, or used directly to
 * construct an equivalent matrix.
 *
 * In the future, the transforms provided here can be improved to
 * include e.g. parallelized versions of the routines, and due to the
 * generic API, these new transforms can be used as drop-in
 * replacements for the old ones.
 *
 * We distinguish two types of transforms: One-dimensional, or just
 * plain, transforms, and multi-dimensional transforms.  This
 * distinction is only for the sake of convenience/practicality and
 * bears no loss of generality: It is easier to think in terms of 1D
 * transforms, and to actually write 1D transforms. When a
 * multidimensional equivalent is needed, the 1D transforms are easily
 * embedded into a multidimensional domain via the provided routines.
 *
 * All transforms designed operate in-place, and transform objects
 * consume a minimal amount of memory, so as to facilitate even large
 * problems where RAM quickly becomes a scarce resource.
 *
 */

/* ******************************************************************************** */

/**
 * One-dimensional transform type. Here, `target` is on input the
 * array to be transformed, on output it will be the transformed data.
 * Parameter `size` is the number of elements in each array to be
 * transformed, with consequtive elements separated by `stride`
 * units. The transform is repeatedly applied to `howmany` vectors,
 * which are in turn separated by `dist` elements (to be very clear,
 * the distance between the initial elements of two consequtive arrays
 * is `dist`).
 */
typedef int (psede_transf_call_t)(double *target, int size, int stride,
				  int howmany, int dist, void *params);

/**
 * Function for finalizing a 1D transf object. This deletes the
 * `params` object in the transf struct.
 */
typedef void (psede_transf_finalize_t)(void*);

/**
 * Wrapper for the 1D transform and associated parameters.  
 */
typedef struct {
  psede_transf_call_t *transform;
  psede_transf_finalize_t *finalize;
  void *params;
} psede_transf_t;

/**
 * Uninitialized transf. All transfs are advised to be set
 * to `psede_tnil` before initialization.
 */
extern const
psede_transf_t psede_tnil;

/**
 * Ones transform. Sets all elements in input to one.
 */
extern const
psede_transf_t psede_transf_ones;

/**
 * Make a (shallow) copy of a transformation. 
 */
void
psede_transf_copy(psede_transf_t *dest,
		  const psede_transf_t *source);

/**
 * Apply a transf. Simply calls the transform function with the
 * parameters supplied in the transf.
 */
int
psede_transf_apply(const psede_transf_t *transf,
		   double *target,
		   int size, int stride, int howmany, int dist);

/**
 * Apply a transf. Short-hand for `psede_transf_apply` but with unit
 * stride, single target.
 */
int
psede_transf_apply_0(const psede_transf_t *transf,
		     double *target,
		     int size);

/**
 * Destroys a transf object. Finalizes the object by calling the
 * finalize function, freeing any objects that initialization
 * allocated. Allocation and deallocation of the transf object is
 * up to the user; as these are rather small, one can safely create
 * them on the stack.
 */
void
psede_transf_destroy(psede_transf_t *transf);


/* ******************************************************************************** */

/**
 * Multidimensional transform type. Similar to the 1D case, but the
 * number of dimensions, `dimension`, and numbers of grid points along
 * each dimension, `sizes`, are need to be given.
 */
typedef int (psede_multitransf_call_t)(double *target,
				       int dimension, const int *sizes,
				       int stride, int howmany,
				       int dist, void *params);

/**
 * Function for finalizing a multidim transf objects. This
 * deletes the `params` object, but does not free the object itself.
 */
typedef void (psede_multitransf_finalize_t)(void*);

/**
 * Wrapper for the multi-dimensional transform and parameters.
 */
typedef struct {
  psede_multitransf_call_t *multitransform;
  psede_multitransf_finalize_t *finalize;
  void *params;
} psede_multitransf_t;

extern const
psede_multitransf_t psede_multitransf_nil;

extern const
psede_multitransf_t psede_multitransf_ones;


/**
 * Make a (shallow) copy of a transformation. 
 */
void
psede_multitransf_copy(psede_multitransf_t *dest,
		       const psede_multitransf_t *source);

/**
 * Apply a multitransf object. Wraps a call to the transform
 * function with parameters supplied in the transf struct.
 */
int
psede_multitransf_apply(const psede_multitransf_t *multitransf,
			double *target,
			int dimension, const int *sizes,
			int stride, int howmany, int dist);

/**
 * Destroy a multitransf object. Finalizes the object by calling
 * the finalize function.
 */
void
psede_multitransf_destroy(psede_multitransf_t *mtf);

/**
 * Convenience function for applying a 1D transform on rows or columns
 * of a tight-packed square matrix. The matrix is assumed to be
 * row-major `size`-by-`size`, lead-dimension `size`, with the
 * transform applied column-wise if `trans` is 0. If `init` is
 * non-zero, the target array is initialized to identity matrix.
 */
int
psede_apply_matrix_0(double *target, int size, int trans, int init,
		     psede_transf_call_t *transform, void *params);

/**
 * Apply a one-dimensional transformation to a multi-dimensional
 * row-major array along a selected axis.
 */
int
psede_apply_multi(double *multi_target, int n_dimensions, const int *sizes,
		  int stride, int howmany, int dist,
		  psede_transf_call_t *transform,
		  int target_dimension, void *params);

/**
 * Convenience function for applying a multi-dimensional transform on
 * rows or columns of a tight-packed square matrix. Extension of
 * `psede_apply_matrix_0`.
 */
int
psede_apply_multi_matrix_0(double *multi_target,
			   int n_dimensions, const int *sizes,
			   int transp, int init,
			   psede_transf_call_t *transform,
			   int target_dimension, void *params);

/* ******************************************************************************** */

/**
 * Initialize a multitransf as a 1D transf `tfmer` applied along
 * selected dimension `dimension`. Deleting the object initialized by
 * this function will not destroy the transf object `tfmer` passed in
 * as argument. Parameter `owner` indicates whether caller wishes to
 * retain ownership of the pointer `tfmer`, or allow the called
 * routine to finalise
 */
int
psede_init_embedding(psede_multitransf_t *mtfmer,
		     psede_transf_t *tfmer, int dimension,
		     psede_ownership_t ownership);


/**
 * Initialize a multitransf to a linear combination of `n_ts`
 * individual transforms, with combination coefficients `ws`.
 *
 * @todo This needs fixing. Currently not conforming to the ownership
 * passing scheme used in similar functions.
 */
int
psede_init_lincomb(psede_multitransf_t *t,
		   int n_ts,
		   psede_multitransf_t **ts,
		   const double *ws);


/**
 * Composition multitransform. Takes n multitransform and sequences
 * them together to form one composite multitransform. 
 */
int
psede_init_composite(psede_multitransf_t *mt,
		     int n_ts,
		     psede_multitransf_t **ts,
		     psede_ownership_t ownership);


#endif /* PSEDE_MULTI_H */
