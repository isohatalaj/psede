#include <stdlib.h>
#include <string.h>

#include "psede_collocation.h"

const
psede_colloc_t psede_colloc_nil = { NULL, NULL, NULL, NULL, NULL, NULL };

const
psede_multicolloc_t psede_multicolloc_nil = { NULL, NULL, NULL, NULL, NULL, NULL };

const psede_colloc_t psede_Tx = {
  psede_init_Tx_diff,
  psede_init_Tx_fct,
  psede_init_Tx_ifct,
  psede_init_Tx_nodes,
  NULL,
  NULL
};

void
psede_colloc_copy(psede_colloc_t *dest,
		  const psede_colloc_t *source)
{
  memcpy(dest, source, sizeof(*dest));
}

int
psede_init_diff(const psede_colloc_t *colloc,
		psede_transf_t *transform,
		int order)
{
  return colloc->init_diff(transform, order, colloc->params);
}


int
psede_init_P2M(const psede_colloc_t *colloc,
	       psede_transf_t *transform)
{
  return colloc->init_P2M(transform, colloc->params);
}

int
psede_init_M2P(const psede_colloc_t *colloc,
	       psede_transf_t *transform)
{
  return colloc->init_M2P(transform, colloc->params);
}

int
psede_init_nodes(const psede_colloc_t *colloc,
		 psede_transf_t *transform)
{
  return colloc->init_nodes(transform, colloc->params);
}

void
psede_finalize(psede_colloc_t *colloc)
{
  colloc->finalize(colloc->params);
}

int
psede_init_multidiff(const psede_multicolloc_t *colloc,
		     psede_multitransf_t *transform,
		     int *orders)
{
  return colloc->init_multidiff(transform, orders, colloc->params);
}


int
psede_init_multiP2M(const psede_multicolloc_t *colloc,
		    psede_multitransf_t *transform)
{
  return colloc->init_multiP2M(transform, colloc->params);
}

int
psede_init_multiM2P(const psede_multicolloc_t *colloc,
		    psede_multitransf_t *transform)
{
  return colloc->init_multiM2P(transform, colloc->params);
}

int
psede_init_multifunc(const psede_multicolloc_t *colloc,
		     psede_multitransf_t *transform,
		     psede_multifunc_t *func,
		     void *func_params)
{
  return colloc->init_multifunc(transform, func, func_params,
				colloc->params);
}


int
psede_init_multinodes(const psede_multicolloc_t *colloc,
		      psede_multitransf_t *transform, int axis)
{
  return colloc->init_multinodes(transform, axis, colloc->params);
}

void
psede_multicolloc_destroy(psede_multicolloc_t *colloc)
{
  colloc->finalize(colloc->params);
}

/* The canonical multicollocation object: A direct product of 1D
 * collocations.
 */
typedef struct {
  int n_dims;
  psede_colloc_t *cs;
  psede_ownership_t cs_caller_policy;
} psede_multicolloc_canon_params_t;


void
psede_multicolloc_canon_params_free(psede_multicolloc_canon_params_t *params)
{
  if (params == NULL) return;

  if (params->cs_caller_policy == PSEDE_PASS_OWNERSHIP)
    {
      int i;
      for (i = 0; i < params->n_dims; ++i)
	{
	  psede_finalize(&params->cs[i]);
	}
    }
  
  free(params->cs);
  free(params);
}

psede_multicolloc_canon_params_t *
psede_multicolloc_canon_params_alloc(int n_dims,
				     const psede_colloc_t **cs,
				     psede_ownership_t ownership)
{
  psede_multicolloc_canon_params_t *params;

  params = malloc(sizeof(*params));
  if (params == NULL) return NULL;

  params->cs = NULL;

  params->n_dims = n_dims;
  params->cs = malloc(n_dims*sizeof(*params->cs));
  if (params->cs == NULL) goto fail;

  int i;
  for (i = 0; i < n_dims; ++i)
    {
      psede_colloc_copy(&params->cs[i], cs[i]);
    }

  params->cs_caller_policy = ownership;

  return params;

 fail:
  psede_multicolloc_canon_params_free(params);

  return NULL;
}

/**
 * Construct the canonical multidimensional differentiation
 * transform.
 *
 * @todo The order of the differentiations is fixed here.  Also,
 * psede_init_canon_multiany, is nearly identical to this -- rather
 * disgusting repetition of complex code. Should be fused somehow.
 */
int
psede_init_canon_multidiff(psede_multitransf_t *mdiff, int *orders,
			   psede_multicolloc_canon_params_t *params)
{
  int status;
  int i;
  int n_composite = 0;

  for (i = 0; i < params->n_dims; ++i)
    {
      if (orders[i] > 0) n_composite++;
    }

  psede_multitransf_t **mt = malloc(n_composite*sizeof(*mt));
  if (mt == NULL) return 0;

  for (i = 0; i < n_composite; ++i) mt[i] = NULL;

  int j = 0;
  for (i = 0; i < params->n_dims; ++i)
    {
      if (orders[i] == 0) continue;

      mt[j] = malloc(sizeof(*mt[i]));
      if (mt[j] == NULL) goto fail;

      psede_transf_t t;
      status = psede_init_diff(&params->cs[i], &t, orders[i]);
      if (status) goto fail;

      status = psede_init_embedding(mt[j], &t, i,
				    PSEDE_PASS_OWNERSHIP);
      if (status)
	{
	  psede_transf_destroy(&t);
	  goto fail;
	}
      
      j++;
    }

  status = psede_init_composite(mdiff, n_composite, mt,
				PSEDE_PASS_OWNERSHIP);
  if (status) goto fail;

  for (i = 0; i < n_composite; ++i) free(mt[i]);
  free(mt);

  return 0;

 fail:

  for (i = 0; i < n_composite; ++i)
    {
      if (mt[i])
	{
	  psede_multitransf_destroy(mt[i]);
	  free(mt[i]);
	}
    }  
  free(mt);

  return 1;
}

int
psede_init_canon_multiany(psede_multitransf_t *many,
			  int (initializer)(const psede_colloc_t *, psede_transf_t*),
			  psede_multicolloc_canon_params_t *params)
{
  int status;
  int i;

  psede_multitransf_t **mt = malloc(params->n_dims*sizeof(*mt));
  if (mt == NULL) return 0;

  for (i = 0; i < params->n_dims; ++i) mt[i] = NULL;


  for (i = 0; i < params->n_dims; ++i)
    {
      mt[i] = malloc(sizeof(*mt[i]));
      if (mt[i] == NULL) goto fail;

      psede_transf_t t;
      status = initializer(&params->cs[i], &t);
      if (status) goto fail;

      status = psede_init_embedding(mt[i], &t, i,
				    PSEDE_PASS_OWNERSHIP);
      if (status)
	{
	  psede_transf_destroy(&t);
	  goto fail;
	}
    }

  status = psede_init_composite(many, params->n_dims, mt,
				PSEDE_PASS_OWNERSHIP);
  if (status) goto fail;

  for (i = 0; i < params->n_dims; ++i) free(mt[i]);
  free(mt);

  return 0;

 fail:

  for (i = 0; i < params->n_dims; ++i)
    {
      if (mt[i])
	{
	  psede_multitransf_destroy(mt[i]);
	  free(mt[i]);
	}
    }
  free(mt);

  return 1;
}

typedef struct {
  psede_multifunc_t *func;
  void *params;

  int dim;
  psede_transf_t *nodes;

  int *sizes;
  double **xs;
  double *x;
  int *ix;
} psede_mftransf_params_t;

void
psede_mftransf_params_free(psede_mftransf_params_t *params)
{
  if (params == NULL) return;

  if (params->nodes)
    {
      int i;
      for (i = 0; i < params->dim; ++i)
	{
	  psede_transf_destroy(&params->nodes[i]);
	  free(params->xs[i]);
	}
      free(params->nodes);
    }
  
  if (params->xs) free(params->xs);
  if (params->x) free(params->x);
  if (params->ix) free(params->ix);
  if (params->sizes) free(params->sizes);

  free(params);
}

int
psede_mftransf_apply(double *target, int dimension, const int *sizes,
		     int stride, int howmany, int dist,
		     psede_mftransf_params_t *params)
{
  if (dimension != params->dim) return 1;

  int status;
  int i, m = 1;
  for (i = 0; i < dimension; ++i)
    {      
      if (params->sizes[i] != sizes[i])
	{
	  double *ptr;
	  ptr = realloc(params->xs[i], sizes[i]*sizeof(*ptr));
	  if (ptr == NULL) return 1;

	  params->xs[i] = ptr;
	  params->sizes[i] = sizes[i];

	  status = psede_transf_apply(&params->nodes[i], params->xs[i],
				      params->sizes[i], 1, 1, 0);
	  if (status) return status;
	}

      params->ix[i] = 0;
      params->x[i] = params->xs[i][0];
      m *= sizes[i];
    }

  int j;
  for (j = 0; j < howmany; ++j)
    {
      double *ptr = target + j*dist;

      for (i = 0; i < m; ++i)
	{
	  double y;

	  status = params->func(params->x, &y,
				params->params);
	  if (status) return status;

	  ptr[i*stride] *= y;

	  int k;
	  for (k = params->dim - 1; k >= 0; --k)
	    {
	      if (params->ix[k] == params->sizes[k] - 1)
		params->ix[k] = 0;
	      else
		params->ix[k]++;

	      params->x[k] = params->xs[k][params->ix[k]];

	      if (params->ix[k] != 0) break;
	    }
	}
    }

  return 0;
}

psede_mftransf_params_t *
psede_mftransf_params_alloc(psede_multifunc_t *func, void *func_params,
			    int dim, const psede_colloc_t **cs)
{
  int status;
  psede_mftransf_params_t *params;

  params = malloc(sizeof(*params));
  if (params == NULL) return NULL;

  params->nodes = NULL;
  params->sizes = NULL;
  params->xs = NULL;
  params->x = NULL;
  params->ix = NULL;

  params->func = func;
  params->params = func_params;
  params->dim = dim;

  params->nodes = malloc(dim*sizeof(*params->nodes));
  if (params->nodes == NULL) goto fail;

  params->sizes = malloc(dim*sizeof(*params->sizes));
  if (params->sizes == NULL) goto fail;

  params->xs = malloc(dim*sizeof(*params->sizes));
  if (params->xs == NULL) goto fail;

  params->x = malloc(dim*sizeof(*params->x));
  if (params->x == NULL) goto fail;

  params->ix = malloc(dim*sizeof(*params->ix));
  if (params->ix == NULL) goto fail;

  int i;
  for (i = 0; i < dim; ++i)
    {
      params->nodes[i] = psede_tnil;
      params->sizes[i] = 0;
      params->xs[i] = NULL;
    }

  for (i = 0; i < dim; ++i)
    {
      status = psede_init_nodes(cs[i], &params->nodes[i]);
      if (status) goto fail;
    }

  return params;

 fail:
  psede_mftransf_params_free(params);

  return NULL;
}

int
psede_init_canon_multiP2M(psede_multitransf_t *mt, 
			  psede_multicolloc_canon_params_t *params)
{
  return psede_init_canon_multiany(mt,
				   psede_init_P2M,
				   params);
}

int
psede_init_canon_multiM2P(psede_multitransf_t *mt, 
			  psede_multicolloc_canon_params_t *params)
{
  return psede_init_canon_multiany(mt, psede_init_M2P, params);
}


int psede_init_canon_multifunc(psede_multitransf_t *mt,
			       psede_multifunc_t *func,
			       void *func_params,
			       psede_multicolloc_canon_params_t *params)
{
  psede_colloc_t **cs;

  cs = malloc(params->n_dims*sizeof(*cs));
  if (cs == NULL) return 1;

  int i;
  for (i = 0; i < params->n_dims; ++i) cs[i] = &params->cs[i];

  psede_mftransf_params_t *transf_params;
  transf_params = psede_mftransf_params_alloc(func, func_params,
					      params->n_dims,
					      (const psede_colloc_t **) cs);
  if (transf_params == NULL) goto fail;

  mt->multitransform = (psede_multitransf_call_t*) psede_mftransf_apply;
  mt->finalize = (psede_finalize_t*) psede_mftransf_params_free;
  mt->params = transf_params;

  free(cs);

  return 0;
  
 fail:

  if (cs) free(cs);

  return 1;
}

/**
 * Function for extracting a coordinate/node value.
 *
 * @todo Make sure that ptr to int cast is kosher.
 */
int
psede_multifunc_pick_x(const double *x, double *y, void *axis)
{
  int index = (int) (long) axis;
  *y = x[index];
  return 0;
}

int
psede_init_canon_multinodes(psede_multitransf_t *mt,
			    int axis,
			    psede_multicolloc_canon_params_t *params)
{
  return psede_init_canon_multifunc(mt, psede_multifunc_pick_x, (void *) (long) axis,
				    params);
}

  
int
psede_init_multicolloc(psede_multicolloc_t *mc,
		       int n_dims,
		       const psede_colloc_t **cs,
		       psede_ownership_t ownership)
{
  psede_multicolloc_canon_params_t *params;

  params = psede_multicolloc_canon_params_alloc(n_dims, cs, ownership);
  if (params == NULL) return 1;

  mc->init_multidiff = (psede_init_multidiff_t*) psede_init_canon_multidiff;
  mc->init_multiP2M = (psede_init_multiP2M_t*) psede_init_canon_multiP2M;
  mc->init_multiM2P = (psede_init_multiM2P_t*) psede_init_canon_multiM2P;
  mc->init_multinodes = (psede_init_multinodes_t*) psede_init_canon_multinodes;
  mc->init_multifunc = (psede_init_multifunc_t*) psede_init_canon_multifunc;
  mc->finalize = (psede_finalize_t*) psede_multicolloc_canon_params_free;
  mc->params = params;

  return 0;
}
