
#include <stdlib.h>
#include <string.h>

#include "psede_multi.h"
#include "psede_util.h"

const psede_transf_t psede_tnil = {NULL, NULL, NULL};
const psede_multitransf_t psede_multitransf_nil = {NULL, NULL, NULL};

int
psede_transf_ones_apply(double *target, int size, int stride, int howmany, int dist,
			void *params)
{
  int i, j;

  /* todo should let blas handle low level stuff */
  for (i = 0; i < howmany; ++i)
    {
      for (j = 0; j < size; ++j)
	{
	  target[i*dist + j*stride] = 1.0;
	}
    }

  return 0;
}

const psede_transf_t psede_transf_ones = {
  psede_transf_ones_apply, NULL, NULL
};

int
psede_multitransf_ones_apply(double *target, int dim, const int *sizes,
			     int stride, int howmany, int dist, void *params)
{
  int m = 1;
  int i;
  for (i = 0; i < dim; ++i) m *= sizes[i];
  return psede_transf_ones_apply(target, m, stride, howmany, dist, params);
}

const psede_multitransf_t psede_multitransf_ones = {
  psede_multitransf_ones_apply, NULL, NULL
};

void
psede_transf_copy(psede_transf_t *dest,
		  const psede_transf_t *source)

{  
  memcpy(dest, source, sizeof *dest);
}

void
psede_multitransf_copy(psede_multitransf_t *dest,
		       const psede_multitransf_t *source)

{  
  memcpy(dest, source, sizeof *dest);
}

void
psede_transf_destroy(psede_transf_t *transf)
{
  if (transf->finalize)
    transf->finalize(transf->params);
}

void
psede_multitransf_destroy(psede_multitransf_t *mtfmer)
{
  if (mtfmer->finalize)
      mtfmer->finalize(mtfmer->params);
}

int
psede_transf_apply(const psede_transf_t *transf,
		   double *target,
		   int size, int stride, int howmany, int dist)
{
  return transf->transform(target, size, stride, howmany, dist,
			   transf->params);
}

int
psede_transf_apply_0(const psede_transf_t *transf,
		     double *target,
		     int size)
{
  return psede_transf_apply(transf, target, size, 1, 1, 0);
}


int
psede_multitransf_apply(const psede_multitransf_t *multitransf,
			double *target,
			int dimension, const int *sizes,
			int stride, int howmany, int dist)
{
  return multitransf->multitransform(target, dimension, sizes,
				     stride, howmany, dist,
				     multitransf->params);
}


int
psede_apply_matrix_0(double *target, int size, int trans, int init,
		     psede_transf_call_t *transform, void *params)
{
  int status;
  
  if (init) psede_identity(target, size, 1, size, size);

  if (trans)
    {
      status = transform(target, size, 1, size, size, params);      
    }
  else
    {
      status = transform(target, size, size, size, 1, params);
    }

  return status;
}

int
psede_apply_multi_matrix_0(double *multi_target,
			   int n_dimensions, const int *sizes,
			   int transp, int init,
			   psede_transf_call_t *transform,
			   int target_dimension, void *params)
{
  int status;
  int i, m = 1;
  for (i = 0; i < n_dimensions; ++i) m *= sizes[i];

  if (init) psede_identity(multi_target, m, 1, m, m);

  if (transp)
    {
      status = psede_apply_multi(multi_target, n_dimensions, sizes,
				 1, m, m,
				 transform, target_dimension, params);
    }
  else
    {
      status = psede_apply_multi(multi_target, n_dimensions, sizes,
				 m, m, 1,
				 transform, target_dimension, params);
    }

  return status;
}

int
psede_apply_multi(double *multi_target, int n_dimensions, const int *sizes,
		  int stride, int howmany, int dist,
		  psede_transf_call_t *transform, int target_dimension, void *params)
{
  int status;
  int n_below = 1, n_above = 1;

  /* Calculate some necessary structure parameters. These could
   * (should?) be precomputed/put into a convenient struct describing
   * the multidim grid, but the performance/convenience gain is likely
   * to be completely negligible, hence I'm KISSing this. */
  int i;
  for (i = 0; i < n_dimensions; ++i)
    {
      if (i < target_dimension) n_above *= sizes[i];
      else if (i > target_dimension) n_below *= sizes[i];
    }

  int j;
  for (j = 0; j < howmany; ++j)
    {
      double *data = multi_target + j*dist;
      
      /* There is a choice on how to execute the transforms on the
       * multigrid. You can either loop over the number of dimensions
       * "above" (n_above) or "below" (n_below) the target dimension, and
       * use the stride and howmany parameters to get everything
       * right. Here, for the time being, I've opted to make the outermost
       * loop the shortest possible, trusting the transform itself to do
       * what it does best. Proper benchmarks should be done to determine
       * optimal choice; some initial testing suggested that it matters
       * little which way you do this. */
      if (n_below < n_above)
	{
	  for (i = 0; i < n_below; ++i)
	    {
	      status = transform(data + stride*i,
				 sizes[target_dimension],
				 stride*n_below, n_above,
				 stride*sizes[target_dimension]*n_below, params);
	      if (status) return status;
	    }
	}
      else
	{
	  for (i = 0; i < n_above; ++i)
	    {
	      status = transform(data + stride*i*sizes[target_dimension]*n_below,
				 sizes[target_dimension],
				 stride*n_below, n_below, stride, params);
	      if (status) return status;
	    }
	}
    }

  return status;
}

/* ******************************************************************************** */

/**
 * Embedding transf. Takes a one-dimensional transform
 * and constructs a multi-dimensional version, where the original
 * transform applies selectively on a given axis (dimension).
 */
typedef struct {
  int dimension;
  psede_transf_t transf;
  psede_ownership_t transf_caller_policy;
} psede_embed_params_t;

/**
 * Multitransform that applies a 1D transform along selected
 * dimension.
 */
int
psede_apply_embed(double *multi_target, int n_dimensions, const int *sizes,
		  int stride, int howmany, int dist,
		  psede_embed_params_t *params)
{
  return psede_apply_multi(multi_target, n_dimensions, sizes,
			   stride, howmany, dist,
			   params->transf.transform,
			   params->dimension,
			   params->transf.params);
}

psede_embed_params_t *
psede_embed_params_alloc(psede_transf_t *tfmer, int dimension,
			 psede_ownership_t ownership)
{
  psede_embed_params_t *params = NULL;
  
  params = malloc(sizeof(*params));
  if (params == NULL) return NULL;

  params->dimension = dimension;
  psede_transf_copy(&params->transf, tfmer);
  params->transf_caller_policy = ownership;

  return params;
}

void
psede_embed_params_free(psede_embed_params_t *params)
{
  if (params)
    {
      if (params->transf_caller_policy == PSEDE_PASS_OWNERSHIP)
	  psede_transf_destroy(&params->transf);
      
      free(params);
    }
}


int
psede_init_embedding(psede_multitransf_t *mt,
		     psede_transf_t *tfmer, int dimension,
		     psede_ownership_t ownership)
{
  psede_embed_params_t *params = NULL;

  params = psede_embed_params_alloc(tfmer, dimension, ownership);
  if (params == NULL) return 1;

  mt->multitransform = (psede_multitransf_call_t*) psede_apply_embed;
  mt->finalize = (psede_multitransf_finalize_t*) psede_embed_params_free;
  mt->params = params;

  return 0;
}

/* ******************************************************************************** */

/**
 * Linear combination transf. Takes n multitransforms and weights
 * and produces an operator that applies each of the transforms to the
 * input and then takes the weighted sum of the resulting outputs.
 */
typedef struct {
  int n_transfs;
  psede_multitransf_t **transfs;
  double *weights;
  int n_work;
  double *work_temp;
  double *work_accum;
} psede_lincomb_params_t;

/**
 * Multitransform that does linear combinations. 
 */
int
psede_apply_lincomb(double *multi_target, int n_dimensions, const int *sizes,
		    int stride, int howmany, int dist,
		    psede_lincomb_params_t *params);

void
psede_finalize_lincomb(psede_lincomb_params_t *params);

static void
psede_lincomb_params_free(psede_lincomb_params_t *params)
{
  if (params)
    {
      if (params->work_temp) free(params->work_temp);
      if (params->work_accum) free(params->work_accum);
      if (params->weights) free(params->weights);
      
      free(params);
    }
}

static psede_lincomb_params_t *
psede_lincomb_params_alloc(int n_ts,
			   psede_multitransf_t **ts,
			   const double *ws)
{
  psede_lincomb_params_t *params = NULL;

  params->weights = NULL;
  params->transfs = NULL;

  params->weights = malloc(n_ts*sizeof(double));
  if (params->weights == NULL) goto fail;

  params->transfs = malloc(n_ts*sizeof(*ts));
  if (params->transfs == NULL) goto fail;

  params->n_transfs = n_ts;
  memcpy(params->transfs, ts, n_ts*sizeof(*ts));
  memcpy(params->weights, ws, n_ts*sizeof(double));

  return params;

 fail:

  psede_lincomb_params_free(params);
  return NULL;
}


int
psede_init_lincomb(psede_multitransf_t *t,
		   int n_ts,
		   psede_multitransf_t **ts,
		   const double *ws)
{
  psede_lincomb_params_t *params = NULL;

  params = psede_lincomb_params_alloc(n_ts, ts, ws);
  if (params == NULL) return 1;

  t->multitransform = (psede_multitransf_call_t*) psede_apply_lincomb;
  t->finalize = (psede_multitransf_finalize_t*) psede_lincomb_params_free;
  t->params = params;

  return 0;
}


int
psede_apply_lincomb(double *multi_target, int n_dimensions, const int *sizes,
		    int stride, int howmany, int dist,
		    psede_lincomb_params_t *params)
{
  int status;
  int m = 1;
  int i, j, k;

  /* There are a few ways of doing this. In all cases we need employ
   * some temporary storage to carry out the sum over the output of
   * the individual transforms. Since each transform needs to use the
   * same input data, we will make a copy of the input for each of the
   * transforms. We do not need to make a copy of all of the input,
   * however, since we can make the howmany transforms individually.
   * This might not be optimal in terms of performance, but saves
   * memory. Alternative would be to make howmany copies of the input
   * arrays and run the transforms on those, or one could even
   * make such full copies for each component transform individually.
   * 
   * Here, we choose to minimize the memory usage, and allocate a
   * minimal amount of work storage, and then looping the transforms
   * howmany times to build the final output. A more performance
   * oriented, parallelized version could perhaps use a different
   * scheme. For transforms running in single thread, the performance
   * boost from more RAM usage would proabably be tiny.
   * 
   */

  for (i = 0; i < n_dimensions; ++i) m *= sizes[i];

  if (params->n_work < m)
    {
      /* TODO: These allocations should be done using an aligned
	 malloc. */

      free(params->work_temp);
      free(params->work_accum);

      params->work_temp = malloc(m*sizeof(double));
      if (params->work_temp == NULL) return 1;

      params->work_accum = malloc(m*sizeof(double));
      if (params->work_accum == NULL) return 1;

      params->n_work = m;
    }

  double *accum = params->work_accum;
  double *temp = params->work_temp;

  for (k = 0; k < howmany; ++k)
    {
      double *data = multi_target + k * dist;

      for (j = 0; j < m; ++j)
	{
	  accum[j] = 0.0;
	  temp[j] = data[j*stride];
	}

      for (j = 0; j < params->n_transfs; ++j)
	{
	  status = psede_multitransf_apply(params->transfs[j],
						temp, n_dimensions, sizes, 1, 1, 0);
	  if (status) return status;

	  for (i = 0; i < m; ++i) accum[i] += params->weights[j]*temp[i];
	}

      for (j = 0; j < m; ++j)
	{
	  data[j*stride] = accum[j];
	}
    }

  return 0;
}

/* ******************************************************************************** */

typedef struct {
  int n_transfs;
  psede_multitransf_t *transfs;
  psede_ownership_t transfs_caller_policy;
} psede_composite_params_t;

void
psede_composite_params_free(psede_composite_params_t *params)
{
  if (params == NULL) return;

  if (params->transfs)
    {

      if (params->transfs_caller_policy == PSEDE_PASS_OWNERSHIP)
	{

	  int i;
	  for (i = 0; i < params->n_transfs; ++i)
	    {
	      psede_multitransf_destroy(&params->transfs[i]);
	    }
	}

      free(params->transfs);
    }

  free(params);
}

psede_composite_params_t *
psede_composite_params_alloc(int n_transfs,
			     psede_multitransf_t **transfs,
			     psede_ownership_t ownership)
{
  psede_composite_params_t *params;

  params = malloc(sizeof(*params));
  if (params == NULL) return NULL;

  params->transfs = malloc(n_transfs*sizeof(*params->transfs));
  if (params->transfs == NULL) goto fail;

  params->n_transfs = n_transfs;

  int i;
  for (i = 0; i < n_transfs; ++i)
    {
      psede_multitransf_copy(&params->transfs[i], transfs[i]);
    }
  
  params->transfs_caller_policy = ownership;
  
  return params;

 fail:

  psede_composite_params_free(params);

  return NULL;
}

int
psede_composite_transform(double *target,
			  int dimension, const int *sizes,
			  int stride, int howmany,
			  int dist, psede_composite_params_t *params)
{
  int status;
  int i;

  for (i = 0; i < params->n_transfs; ++i)
    {
      status = psede_multitransf_apply(&params->transfs[i],
				       target, dimension, sizes,
				       stride, howmany, dist);

      if (status) return status;
    }

  return 0;
}


int
psede_init_composite(psede_multitransf_t *mt,
		     int n_ts,
		     psede_multitransf_t **ts,
		     psede_ownership_t ownership)
{
  psede_composite_params_t *params;

  params = psede_composite_params_alloc(n_ts, ts, ownership);
  if (params == NULL) return 1;

  mt->multitransform = (psede_multitransf_call_t*) psede_composite_transform;
  mt->finalize = (psede_multitransf_finalize_t*) psede_composite_params_free;
  mt->params = params;

  return 0;  
}
