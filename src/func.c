
#include <stdlib.h>

#include "psede_func.h"
#include "psede_util.h"




typedef struct {
  psede_func_t *func;
  void *params;

  int size;
  psede_transf_t nodes;
  double *xs;
} psede_functransf_params_t;

int
psede_func_transform(double *target, int size, int stride,
		     int howmany, int dist, psede_functransf_params_t *params)
{
  int status;
  int i, j;

  if (params->size != size)
    {
      double *ptr;
      ptr = realloc(params->xs, size*sizeof(*params->xs));
      if (ptr == NULL) return 1;

      params->xs = ptr;
      params->size = size;

      status = psede_transf_apply(&params->nodes, params->xs, size, 1, 1, 0);
      if (status) return status;

      for (j = 0; j < size; ++j)
	{
	  double y;
	  status = params->func(params->xs[j], &y, params->params);
	  if (status) return status;
	  params->xs[j] = y;
	}
    }

  /* todo: delegate this to blas. */
  for (i = 0; i < howmany; ++i)
    {
      double *ptr = target + i*dist;
      for (j = 0; j < size; ++j) ptr[j*stride] *= params->xs[j];
    }

  return 0;
}

void
psede_functransf_params_free(psede_functransf_params_t *params)
{
  if (params == NULL) return;

  psede_transf_destroy(&params->nodes);

  free(params);
}

psede_functransf_params_t *
psede_functransf_params_alloc(psede_func_t *func,
			      void *func_params,
			      const psede_colloc_t *colloc)
{
  psede_functransf_params_t *params = NULL;

  params = malloc(sizeof(*params));
  if (params == NULL) return NULL;

  params->func = func;
  params->params = func_params;
  params->size = 0;
  params->nodes = psede_tnil;
  params->xs = NULL;

  int status;
  status = psede_init_nodes(colloc, &params->nodes);
  if (status) goto fail;

  return params;

 fail:
  psede_functransf_params_free(params);

  return NULL;
}

int
psede_init_func(const psede_colloc_t *colloc,
		psede_transf_t *t,
		psede_func_t *func,
		void *func_params)
{
  psede_functransf_params_t *params = NULL;

  params = psede_functransf_params_alloc(func, func_params, colloc);
  if (params == NULL) return 1;

  t->transform = (psede_transf_call_t*) psede_func_transform;
  t->finalize = (psede_transf_finalize_t*) psede_functransf_params_free;
  t->params = params;

  return 0;
}

int
psede_function_multiply_multi_0(double *mat,
				double (*func)(const double *x, void *params),
				void *params,
				double **nodes,
				int dims, const int *sizes,
				int init)
{
  int status = 0;
  int *i = NULL, j, k, m = 1;
  double *x = NULL;

  i = malloc(dims*sizeof(*i));
  if (i == NULL) { status = 1; goto exit; }

  x = malloc(dims*sizeof(*x));
  if (x == NULL) { status = 1; goto exit; }

  for (j = 0; j < dims; ++j) i[j] = 0;
  for (j = 0; j < dims; ++j) x[j] = nodes[j][0];
  for (j = 0; j < dims; ++j) m *= sizes[j];

  if (init) psede_identity(mat, m, 1, m, m);

  k = 0;
  int loop = 1;
  while (loop)
    {
      const double f = func(x, params);

      if (init)
	{
	  mat[k*m + k] = f;
	}
      else
	{
	  for (j = 0; j < m; ++j) mat[k*m + j] *= f;
	}

      /* Increment all indices (last is fastest varying), and if all
	 roll back to zero, exit the main loop */
      loop = 0;
      for (j = dims - 1; j >= 0 && !loop; --j)
	{
	  i[j]++;
	  if (i[j] == sizes[j]) i[j] = 0;
	  else loop = 1;

	  x[j] = nodes[j][i[j]];
	}

      k++;
    }

  /* Internal consistency check: Should have processed all m rows of
     the input matrix. */
  if (k != m) { status = 1; goto exit; }
  
 exit:
  if (i) free(i); 
  if (x) free(x);
  
  return status;
}




  

