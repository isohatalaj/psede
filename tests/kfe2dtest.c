
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "psede.h"


void
driftdiff(const double r[2], void *params,
	  double mu_out[2], double sigma_out[2*2])
{
  const double sigma = 1.0;
  
  mu_out[0] = -r[0];
  mu_out[1] = r[1];

  sigma_out[0] = sigma; sigma_out[1] = 0.0;
  sigma_out[2] = sigma; sigma_out[3] = sigma;
}



double mu_x(const double r[2], void *params)
{
  double mu_out[2], sigma_out[2*2];
  driftdiff(r, params, mu_out, sigma_out);
  return mu_out[0];
}

double mu_y(const double r[2], void *params)
{
  double mu_out[2], sigma_out[2*2];
  driftdiff(r, params, mu_out, sigma_out);
  return mu_out[1];
}

double Diffs_x(const double r[2], void *params)
{
  double mu_out[2], sigma_out[2*2];
  driftdiff(r, params, mu_out, sigma_out);
  return sigma_out[0]*sigma_out[0] + sigma_out[1]*sigma_out[1];
}

double Diffs_y(const double r[2], void *params)
{
  double mu_out[2], sigma_out[2*2];
  driftdiff(r, params, mu_out, sigma_out);
  return sigma_out[2]*sigma_out[2] + sigma_out[3]*sigma_out[3];
}

double Corrs_xy(const double r[2], void *params)
{
  double mu_out[2], sigma_out[2*2];
  driftdiff(r, params, mu_out, sigma_out);
  return sigma_out[0]*sigma_out[2] + sigma_out[1]*sigma_out[3];
}


int min(int i, int j) { return i < j ? i : j; }


/* For clarity: */
#define KEEP_INPUT 0
#define INPUT_TO_IDENTITY 1
#define VAR_X 0 
#define VAR_Y 1


int
main()
{
  int status = 0;

  const int dim = 2;
  const int n[2] = {32, 32};
  const int m = n[0]*n[1];
  const int m2 = m*m;
  double *nodes[2] = {NULL, NULL};

  int i, j, k, l;

  psede_fct_t *fct = NULL;
  double *LOp = NULL;
  double *Tmp = NULL;
  double *gamma = NULL, *resid = NULL;
  double *JxOp = NULL, *JyOp = NULL;
  double *Jx_sol = NULL, *Jy_sol = NULL;

  void *params = NULL;

  /* ******************************************************************************** */
  /* INITIALIZATION */
  
  fct = psede_fct_alloc();
  if (fct == NULL) { printf("# Out-of-memory.\n"); status = 1; goto exit; }

  for (i = 0; i < dim; ++i)
    {
      nodes[i] = psede_fct_alloc_array(n[i]);
      if (nodes[i] == NULL) { printf("# Out-of-memory.\n"); status = 1; goto exit; }
      
      psede_Tx_nodes_0(nodes[i], n[i]);
    }

  LOp = psede_fct_alloc_array(m2);    /* Linear operator corresponding to the KFE and boundary conditions */
  Tmp = psede_fct_alloc_array(m2);    /* Work memory for temporary objects */
  JxOp  = psede_fct_alloc_array(m2);  /* Prob. current operator, x-direction */
  JyOp  = psede_fct_alloc_array(m2);  /* Prob. current operator, y-direction */

  /* gamma will be the rhs of the linear equation determining the
     solution, LOp f == gamma. It's mostly zero, except for one
     element, that will enforce normalisation. See below. */
  gamma = psede_fct_alloc_array(m);
  resid = psede_fct_alloc_array(m);
  Jx_sol = psede_fct_alloc_array(m);
  Jy_sol = psede_fct_alloc_array(m);

  if (LOp == NULL || Tmp == NULL ||
      JxOp  == NULL || JyOp  == NULL || gamma == NULL || resid == NULL)
    { printf("# Out of memory allocating main matrices.\n"); status = 1; goto exit; }

  /* ******************************************************************************** */
  /* PROBLEM SETUP */

  /* Construct probability current operators JxOp and JyOp */
  /* JxOp f := mu_x f - D_x (1/2 D_x f) - D_y (1/2 C f) */
  /* D_x = sigma_x1^2 + sigma_x2^2,  C = sigma_x1 sigma_y1 + sigma_x2 sigma_y2 */
  psede_function_multiply_multi_0(JxOp, mu_x, params, nodes, dim, n, INPUT_TO_IDENTITY);  /* JxOp <- mu_x */

  psede_function_multiply_multi_0(Tmp, Diffs_x, params, nodes, dim, n, INPUT_TO_IDENTITY); /* Tmp <- Diffs_x */
  psede_Tx_diff_point_matrix_multi_0(Tmp, VAR_X, dim, n, KEEP_INPUT, fct); /* Tmp <- Dx Tmp */
  psede_daxpy_0(m2, -0.5, Tmp, JxOp); /* JxOp <- -0.5*Tmp + JxOp */

  psede_function_multiply_multi_0(Tmp, Corrs_xy, params, nodes, dim, n, INPUT_TO_IDENTITY);
  psede_Tx_diff_point_matrix_multi_0(Tmp, VAR_Y, dim, n, KEEP_INPUT, fct);
  psede_daxpy_0(m2, -0.5, Tmp, JxOp);

  /* JyOp f := mu_y f - D_x (1/2 C f) - D_y (1/2  f) */
  psede_function_multiply_multi_0(JyOp, mu_y, params, nodes, dim, n, INPUT_TO_IDENTITY);

  psede_function_multiply_multi_0(Tmp, Diffs_y, params, nodes, dim, n, INPUT_TO_IDENTITY);
  psede_Tx_diff_point_matrix_multi_0(Tmp, VAR_Y, dim, n, KEEP_INPUT, fct);
  psede_daxpy_0(m2, -0.5, Tmp, JyOp);

  psede_function_multiply_multi_0(Tmp, Corrs_xy, params, nodes, dim, n, INPUT_TO_IDENTITY);
  psede_Tx_diff_point_matrix_multi_0(Tmp, VAR_X, dim, n, KEEP_INPUT, fct);
  psede_daxpy_0(m2, -0.5, Tmp, JyOp);

  /* Construct Kolmogorov Forward Equation operator LOp */
  /* LOp := Dx JxOp + Dy JyOp */
  memcpy(LOp, JxOp, m2*sizeof(double)); /* LOp <- JxOp */
  memcpy(Tmp, JyOp, m2*sizeof(double)); /* Tmp <- JyOp */
  
  psede_Tx_diff_point_matrix_multi_0(LOp, VAR_X, dim, n, KEEP_INPUT, fct); /* LOp <- Dx LOp */ 
  psede_Tx_diff_point_matrix_multi_0(Tmp, VAR_Y, dim, n, KEEP_INPUT, fct); /* Tmp <- Dy Tmp */

  psede_daxpy_0(m2, 1.0, Tmp, LOp); /* LOp <- LOp + Tmp */



  /* Enforce boundary conditions: No current shall flow perpendicular to boundaries of [-1,1]^2 */
  /* (1)  (J f)(x = -1, y).(1, 0) = 0, for all y in [-1,1] */
  /* (2)  (J f)(x = +1, y).(1, 0) = 0, for all y in [-1,1] */
  /* (3)  (J f)(x, y = -1).(0, 1) = 0, for all x in [-1,1] */
  /* (4)  (J f)(x, y = +1).(0, 1) = 0, for all x in [-1,1] */

  /* (1) */
  /* (-1, y) determined by rows (0 + n[1]*0, 1 + n[1]*0, ..., n[1] - 1 + n[1]*0) in LOp.
   * We overwrite those with corresponding rows from JxOp to make the x-current vanish.
   */
  i = 0;
  for (j = 0; j < n[1]; ++j) 
    {
      k = j + n[1]*i;
      for (l = 0; l < m; ++l) LOp[l + m*k] = -JxOp[l + m*k];
    }

  /* (2) */
  /* (+1, y) determined by rows (0 + n[1]*(n[0]-1), 1 + n[1]*(n[0]-1), ..., n[1] - 1 + n[1]*(n[0]-1)) in LOp.
   * Overwrite these as well.
   */
  i = n[0] - 1;
  for (j = 0; j < n[1]; ++j) 
    {
      k = j + n[1]*i;
      for (l = 0; l < m; ++l) LOp[l + m*k] = JxOp[l + m*k];      
    }
  
  /* (3) and (4) */
  /* Similar to above */
  j = 0;
  for (i = 0; i < n[0]; ++i)
    {
      k = j + n[1]*i;
      for (l = 0; l < m; ++l) LOp[l + m*k] = -JyOp[l + m*k];      
    }

  j = n[1] - 1;
  for (i = 0; i < n[0]; ++i)
    {
      k = j + n[1]*i;
      for (l = 0; l < m; ++l) LOp[l + m*k] = JyOp[l + m*k];      
    }

  i = 0;
  j = 0;
  k = j + n[1]*i;
  for (l = 0; l < m; ++l) LOp[l + m*k] = -JxOp[l + m*k] - JyOp[l + m*k];

  i = 0;
  j = n[1] - 1;
  k = j + n[1]*i;
  for (l = 0; l < m; ++l) LOp[l + m*k] = -JxOp[l + m*k] + JyOp[l + m*k];

  i = n[0] - 1;
  j = 0;
  k = j + n[1]*i;
  for (l = 0; l < m; ++l) LOp[l + m*k] = JxOp[l + m*k] - JyOp[l + m*k];

  i = n[0] - 1;
  j = n[1] - 1;
  k = j + n[1]*i;
  for (l = 0; l < m; ++l) LOp[l + m*k] = JxOp[l + m*k] + JyOp[l + m*k];

  for (i = 0; i < m; ++i) gamma[i] = 0;

  /* Enforce normalisation: Replace one row of LOp with the last row 
   * of integration operator. */
  /* Construct integration operator into Tmp. TODO: Quite inefficient,
     as we do not need the full operator, just the last row would
     suffice.*/
  psede_Tx_integ_point_matrix_multi_0(Tmp, VAR_X, dim, n, INPUT_TO_IDENTITY, fct); 
  psede_Tx_integ_point_matrix_multi_0(Tmp, VAR_Y, dim, n, KEEP_INPUT, fct);

  /* Now pick a row in LOp and overwrite it with the last row from Tmp */
  j = n[1]/4; i = n[0]/4;
  k = j + n[1]*i;
  int fix_ind = k;
  for (l = 0; l < m; ++l) LOp[l + m*k] = Tmp[l + m*0];

  /* Unlike all other rows, this should yield a one instead of zero. So: */
  gamma[k] = 1.0;

  /* ******************************************************************************** */
  /* NUMERICAL SOLUTION OF THE LINEAR EQUATIONS  */

  /* Make a copy of LOp, as it will be destroyed upon solving the
     linear equation. */
  memcpy(Tmp, LOp, m2*sizeof(double));
  
  /* We can now solve the equation LOp.f = gamma. */

  psede_linsolve_work_t *work = psede_linsolve_work_alloc(m);
  if (work == NULL)
    {
      printf("# Out of memory allocating solver workspace.\n");
      status = 1;
      goto exit;
    }

  status = psede_linsolve_solve_0(m, LOp, gamma, work);
  if (status) { printf("# Linsolve failed.\n"); goto exit; }

  memcpy(LOp, Tmp, m2*sizeof(double));

  /* ******************************************************************************** */
  /* PRINT RESULTS  */
  
  /* Dump solution to standard output. */
  printf("# Solve succesfull. Result:\n");

  /* Apply the current operators and the KFE itself to the obtained
   * solution to find the probability currents and residual. */
  psede_gemv_0(m, 1.0, JyOp, gamma, 0.0, Jy_sol);
  psede_gemv_0(m, 1.0, JxOp, gamma, 0.0, Jx_sol);
  psede_gemv_0(m, 1.0, LOp, gamma, 0.0, resid);
  resid[fix_ind] -= 1.0;

  double integ_x = 0.0; /* Calculate integral to check normalisation.
			   We could use the pseudospec. integral
			   operator, but for redundancy, lets verify
			   using trapezoidal rule. */
  double prev_if = 0.0;
  double prev_x = 0.0;

  double resid_max = 0.0;
    
  k = 0;
  for (i = 0; i < n[0]; ++i)
    {
      const double x = nodes[0][i];
      
      double integ_y = 0.0;
      double prev_f = 0.0;
      double prev_y = 0.0;
      
      for (j = 0; j < n[1]; ++j)
  	{
	  const double y = nodes[1][j];
	  
  	  if (k % 20 == 0) printf("#%25s %25s %25s %25s %25s %25s\n",
				  "x", "y", "f(x,y)", "log10 |resid|", "Jx", "Jy");
  	  k++;

	  const int ijx = j + i*n[1];
	  const double f = gamma[ijx];

	  if (j > 0) integ_y += 0.5*(prev_f + f) * fabs(prev_y - y);
	  prev_y = y;
	  prev_f = f;

	  const double r = log10(fabs(resid[ijx]));

  	  printf(" %25.15lf %25.15lf %25.15lf %25.15lf %25.15lf %25.15lf\n",
  		 x, y, f, r,
		 Jx_sol[ijx], Jy_sol[ijx]);

	  if ((i == 0 && j == 0) || resid_max < r) resid_max = r;
  	}

      if (i > 0) integ_x += 0.5*(prev_if + integ_y) * fabs(prev_x - x);
      prev_x = x;
      prev_if = integ_y;
    }

  printf("# Normalisation = %lf\n", integ_x);
  printf("# |resid|       = %lf\n", resid_max);

  /* ******************************************************************************** */
  /* CLEANUP AND EXIT */
  
 exit:

  for (i = 0; i < dim; ++i) if (nodes[i]) psede_fct_free_array(nodes[i]);

  if (fct) psede_fct_free(fct);

  if (gamma) psede_fct_free_array(gamma);
  if (resid) psede_fct_free_array(resid);
  if (LOp) psede_fct_free_array(LOp);
  if (Tmp) psede_fct_free_array(Tmp);
  if (JxOp) psede_fct_free_array(JxOp);
  if (JyOp) psede_fct_free_array(JyOp);
  
  if (status) printf("# Test failed.\n");
}


