
#include <math.h>

#include "psede.h"

void
psede_nodes(double *x, int n, int stride)
{
  int i;
  const double z = M_PI/(n-1);
  for (i = 0; i < n; ++i)
    {
      x[i*stride] = cos(i*z);
    }
}
