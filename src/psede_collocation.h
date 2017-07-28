/**
 * Provides an abstraction for different collocation methods.
 */

#ifndef PSEDE_COLLOCATION_H
#define PSEDE_COLLOCATION_H

#include "psede_multi.h"
#include "psede_fct.h"
#include "psede_diff.h"
#include "psede_util.h"

typedef int (psede_init_diff_t)(psede_transf_t *transform, int order, void *params);
typedef int (psede_init_P2M_t)(psede_transf_t *transform, void *params);
typedef int (psede_init_M2P_t)(psede_transf_t *transform, void *params);
typedef int (psede_init_nodes_t)(psede_transf_t *transform, void *params);

/**
 * Definition of a collocation method. Provides constructors for key
 * transforms, namely differentiation/integration and mappings between
 * point (node) and mode representations.
 *
 */
typedef struct {

  psede_init_diff_t *init_diff;
  psede_init_P2M_t *init_P2M;
  psede_init_M2P_t *init_M2P;
  psede_init_nodes_t *init_nodes;

  psede_finalize_t *finalize;  

  void *params;
  
} psede_colloc_t;

extern const psede_colloc_t psede_colloc_nil;

void
psede_colloc_copy(psede_colloc_t *dest,
		  const psede_colloc_t *source);


typedef int (psede_multifunc_t)(const double *x, double *y, void *params);


struct psede_multicolloc;
typedef struct psede_multicolloc psede_multicolloc_t;

int
psede_init_multifunc(const psede_multicolloc_t *colloc,
		     psede_multitransf_t *transform,
		     psede_multifunc_t *func,
		     void *func_params);


typedef int (psede_init_multidiff_t)(psede_multitransf_t *mt, int *orders, void *params);
typedef int (psede_init_multiP2M_t)(psede_multitransf_t *mt, void *params);
typedef int (psede_init_multiM2P_t)(psede_multitransf_t *mt, void *params);
typedef int (psede_init_multinodes_t)(psede_multitransf_t *mt, int axis, void *params);
typedef int (psede_init_multifunc_t)(psede_multitransf_t *mt,
				     psede_multifunc_t *func,
				     void *func_params,
				     void *params);

/**
 * Multicollocation type. Encapsulates functions that construct the
 * multidimensional equivalents of the transforms provided by
 * `psede_colloc_t`. Notable differences is the nodes transform which
 * here is replaced by the function application. For the 1D case, the
 * nodes transform provides means for doing function applications on
 * vectors; this approach becomes trickier in the multidimensional
 * case, hence the provided function transform and the absence of a
 * nodes transform.
 *
 * @todo A proper nodes transform should be included.
 *
 */
struct psede_multicolloc {

  psede_init_multidiff_t *init_multidiff;
  psede_init_multiP2M_t *init_multiP2M;
  psede_init_multiM2P_t *init_multiM2P;
  psede_init_multinodes_t *init_multinodes; 
  psede_init_multifunc_t *init_multifunc;

  psede_finalize_t *finalize;  

  void *params;

};

extern const psede_multicolloc_t psede_multicolloc_nil;

int
psede_init_multicolloc(psede_multicolloc_t *mc,
		       int n_dims,
		       const psede_colloc_t **cs,
		       psede_ownership_t ownersip);


extern const
psede_colloc_t psede_Tx;


/**
 * Construct a collocation
 */
int
psede_colloc_init_scaled_Tx(psede_colloc_t *colloc,
			    double x0, double x1);


int
psede_init_diff(const psede_colloc_t *colloc,
		psede_transf_t *transform, int order);


int
psede_init_P2M(const psede_colloc_t *colloc,
	       psede_transf_t *transform);

int
psede_init_M2P(const psede_colloc_t *colloc,
	       psede_transf_t *transform);

int
psede_init_nodes(const psede_colloc_t *colloc,
		 psede_transf_t *transform);

void
psede_finalize(psede_colloc_t *colloc);

int
psede_init_multidiff(const psede_multicolloc_t *colloc,
		     psede_multitransf_t *transform, int *orders);

int
psede_init_multiP2M(const psede_multicolloc_t *colloc,
		    psede_multitransf_t *transform);

int
psede_init_multiM2P(const psede_multicolloc_t *colloc,
		    psede_multitransf_t *transform);

int
psede_init_multinodes(const psede_multicolloc_t *colloc,
		      psede_multitransf_t *transform, int axis);

int
psede_init_multifunc(const psede_multicolloc_t *colloc,
		     psede_multitransf_t *transform,
		     psede_multifunc_t *func,
		     void *func_params);

void
psede_multicolloc_destroy(psede_multicolloc_t *colloc);

#endif /* PSEDE_COLLOCATION_H */
