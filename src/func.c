
#include <stdlib.h>

#include "psede_func.h"

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
