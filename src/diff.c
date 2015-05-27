
#include "psede_diff.h"

void
psede_diff_mode_apply(double *x, int n, int stride, 
		      int howmany, int dist)
{
  int i, k;
  double xi, xii;

  for (k = 0; k < howmany; ++k)
    {
      xi = x[(n-1)*stride];
      x[(n-1)*stride] = 0.0;
      
      xii = x[(n-2)*stride];
      x[(n-2)*stride] = 2.0*(n-1)*xi;
      xi = xii;

      for (i = n - 3; i > 0; --i)
	{
	  xii = x[i*stride];
	  x[i*stride] = x[(i+2)*stride] + 2.0*(i+1)*xi;
	  xi = xii;
	}

      x[0] = 0.5*x[2*stride] + xi;

      x += dist;
    }
}

int
psede_diff_point_apply(double *x, int size, int stride,
		       int howmany, int dist, psede_fct_t *fct)
{
  int status;

  status = psede_fct_apply(fct, x, size, stride, howmany, dist);
  if (status) return status;

  psede_diff_mode_apply(x, size, stride, howmany, dist);

  status = psede_fct_apply_inv(fct, x, size, stride, howmany, dist);
  if (status) return status;

  return 0;
}

static void
id_matrix(double *id, int size, int stride, int howmany, int dist)
{
  int i, j;
  
  for (i = 0; i < size; ++i)
    {
      for (j = 0; j < i; ++j)
	{
	  id[j*dist + i*stride] = 0.0;
	}

      id[j*dist + i*stride] = 1.0;

      for (j = i + 1; j < size; ++j)
	{
	  id[j*dist + i*stride] = 0.0;
	}
    }
}

/* TODO: The performance of these matrix construction functions is in
 * all likelyhood sub-optimal. Getting these absolutely right is going
 * to need some benchmarking, which would probably more effort than
 * its worth -- the complexity of these routines is O(n^2) and
 * O(n^2*log(n)) which (I think) is as good as it is going to get.
 */
void
psede_diff_mode_matrix(double *diff_mode, int size, int stride,
		       int howmany, int dist)
{
  id_matrix(diff_mode, size, stride, howmany, dist);
  psede_diff_mode_apply(diff_mode, size, dist, howmany, stride);
}

int
psede_diff_point_matrix(double *diff_point, int size, int stride,
			int howmany, int dist,
			psede_fct_t *fct)
{
  int status;
  id_matrix(diff_point, size, stride, howmany, dist);
  status = psede_diff_point_apply(diff_point, size, dist, 
				  howmany, stride,
				  fct);

  return status;
}


