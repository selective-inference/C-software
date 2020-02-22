#include <math.h> // for exp, log
#include <stdio.h>
// This provides (as a function of linear predictor), the Cox partial likelihood, 
// its gradient and its Hessian acting on an n \times k matrix on the right.

// This forms the vector of weights
// W_i = \sum_{j: T_j \geq T_{(i)}} \exp(\eta_j)
// so it is ordered according to argsort(T)

void _update_cox_weights(double *linear_pred_ptr, /* Linear term in objective */
			 double *weight_ptr,      /* accumulated weights */
			 double *weight_denom_ptr,/* accumulated denom terms */
			 long *censoring_ptr,     /* censoring indicator */
			 long *ordering_ptr,      /* 0-based ordering of times */
			 long *rankmin_ptr,       /* 0-based ranking with min tie breaking */
			 long ncase               /* how many subjects / times */
			 )       
{
  long idx;
  long order_idx, rankmin_idx;
  double linear_pred;
  long censoring;
  double *weight, *weight_denom;
  double cur_val = 0;

  // reversed reverse cumsum of exp(eta)

  for (idx=0; idx<ncase; idx++) {
    order_idx = *((long *) ordering_ptr + (ncase - 1 - idx));
    linear_pred = *((double *) linear_pred_ptr + order_idx);
    cur_val = cur_val + exp(linear_pred);
    weight = ((double *) weight_ptr + (ncase - 1 - idx));
    *weight = cur_val;
    fprintf(stderr, "%d %d, %f\n", order_idx, ncase - 1 - idx, exp(linear_pred));
  }

  // accumulate inverse cumsums at rankmin
  // i-th value is \sum_{j=1}^i 1 / W(r[o[j]])  where r is rankmin of times so r[o]
  // is rankmin of ordered times

  cur_val = 0;
  for (idx=0; idx<ncase; idx++) {
    order_idx = *((long *) ordering_ptr + idx);
    rankmin_idx = *((long *) rankmin_ptr + order_idx);
    weight = ((double *) weight_ptr + rankmin_idx);
    censoring = *((long *) censoring_ptr + order_idx);
    cur_val = cur_val + censoring / *weight;
    weight_denom = ((double *) weight_denom_ptr + idx);
    *weight_denom = cur_val;
  }

}

// Objective value

double _cox_objective(double *linear_pred_ptr, /* Linear term in objective */
		      double *weight_ptr,      /* accumulated weights */
		      double *weight_denom_ptr,/* accumulated inverse weights */
		      long *censoring_ptr,     /* censoring indicator */
		      long *ordering_ptr,      /* 0-based ordering of times */
		      long *rankmin_ptr,       /* 0-based ranking with min tie breaking */
		      long *rankmax_ptr,       /* 0-based ranking with max tie breaking */
		      long ncase               /* how many subjects / times */
		      )       
{
  long idx, rankmin_idx;
  double linear_pred, weight;
  long censoring;
  double cur_val = 0;

  _update_cox_weights(linear_pred_ptr,
		      weight_ptr,
		      weight_denom_ptr,
		      censoring_ptr,
		      ordering_ptr,
		      rankmin_ptr,
		      ncase);

  for (idx=0; idx<ncase; idx++) {
    rankmin_idx = *((long *) rankmin_ptr + idx);
    weight = *((double *) weight_ptr + rankmin_idx);
    censoring = *((long *) censoring_ptr + idx);
    linear_pred = *((double *) linear_pred_ptr + idx);
    cur_val += (censoring) * (linear_pred - log(weight));
  }

  return(cur_val);

}

void _cox_gradient(double *gradient_ptr,    /* Where gradient is stored */
		   double *linear_pred_ptr, /* Linear term in objective */
		   double *weight_ptr,      /* accumulated weights */
		   double *weight_denom_ptr,/* accumulated inverse weights */
		   long *censoring_ptr,     /* censoring indicator */
		   long *ordering_ptr,      /* 0-based ordering of times */
		   long *rankmin_ptr,       /* 0-based ranking with min tie breaking */
		   long *rankmax_ptr,       /* 0-based ranking with max tie breaking */
		   long ncase               /* how many subjects / times */
		   )       
{
  long idx, rankmax_idx;
  double linear_pred, weight_denom;
  double *gradient;
  long censoring;

  _update_cox_weights(linear_pred_ptr,
		      weight_ptr,
		      weight_denom_ptr,
		      censoring_ptr,
		      ordering_ptr,
		      rankmin_ptr,
		      ncase);

  // fill in entries of gradient

  for (idx=0; idx<ncase; idx++) {
    censoring = *((long *) censoring_ptr + idx);
    rankmax_idx = *((long *) rankmax_ptr + idx);
    weight_denom = *((double *) weight_denom_ptr + rankmax_idx);
    linear_pred = *((double *) linear_pred_ptr + idx);
    gradient = ((double *) gradient_ptr + idx);
    *gradient = censoring - weight_denom * exp(linear_pred);
  }

}

