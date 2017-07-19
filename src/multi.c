
#include "psede_multi.h"

int
psede_apply_matrix_0(double *target, int size, int trans, int init,
		     psede_transform_t *transform, void *params)
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
psede_apply_multi_matrix_0(double *multi_target, int n_dimensions, const int *sizes,
			   int transp, int init,
			   psede_transform_t *transform, int target_dimension, void *params)
{
  int status;
  int i, m = 1;
  for (i = 0; i < n_dimensions; ++i) m *= sizes[i];

  if (init) psede_identity(multi_target, m, 1, m, m);

  /* TODO: Make these transpositions unnecessary */
  /* Note that we need to transpose if transp is false...  */
  if (!transp) psede_transp_0(m, multi_target);
  
  for (i = 0; i < m; ++i)
    {
      status = psede_apply_multi(multi_target + i*m, n_dimensions, sizes,
				 transform, target_dimension, params);
      if (status) return status;
    }

  if (!transp) psede_transp_0(m, multi_target);
  
  return 0;
}

int
psede_apply_multi(double *multi_target, int n_dimensions, const int *sizes,
		  psede_transform_t *transform, int target_dimension, void *params)
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
	  status = transform(multi_target + i,
			     sizes[target_dimension],
			     n_below, n_above,
			     sizes[target_dimension]*n_below, params);
	  if (status) return status;
	}
    }
  else
    {
      for (i = 0; i < n_above; ++i)
  	{
	  status = transform(multi_target + i*sizes[target_dimension]*n_below,
			     sizes[target_dimension],
			     n_below, n_below, 1, params);
  	  if (status) return status;
  	}
    }

  return status;
}
