
#include <math.h>

#include "psede.h"

void
psede_Tx_nodes(double *x, int n, int stride, int howmany, int dist)
{
  int i, j;
  const double z = M_PI/(n-1);

  for (i = 0; i < n; ++i)
    {
      x[i*stride] = cos(i*z);
    }

  for (j = 1; j < howmany; ++j)
    {
      for (i = 0; i < n; ++i)
	{
	  x[i*stride + j*dist] = x[i*stride];
	}
    }
}


void
psede_Tx_nodes_0(double *x, int n)
{
  psede_Tx_nodes(x, n, 1, 1, n);
}


void
psede_Tx_nodes_multi_0(double *x, int dims, const int *n)
{
  int dim;
  int stride = 1; /* stride in x,y,z,.. tuple units */

  int m = 1; /* size of x in x,y,z,... tuples */
  for (dim = 0; dim < dims; ++dim) m *= n[dim];

  for (dim = dims - 1; dim >= 0; --dim)
    {
      /* Do one set of points... */
      psede_Tx_nodes(x + dim, n[dim], dims*stride, 1, n[dim]);

      /* ..and then duplicate it for the rest of the array. */
      int i, j = 0, k = 0;

      for (i = 0; i < m; ++i)
	{
	  x[dims*i + dim] = x[dims*stride*k + dim];

	  j++;
	  if (j == stride)
	    {
	      j = 0;
	      k++;
	      if (k == n[dim])
		{
		  k = 0;
		}
	    }
	}

      stride *= n[dim];
    }  
}


