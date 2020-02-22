#include <math.h> // for exp, log

// This provides (as a function of linear predictor), the Cox partial likelihood, 
// its gradient and its Hessian acting on an n \times k matrix on the right.

// This forms the vector of weights
// W_i = \sum_{j: T_j \geq T_{(i)}} \exp(\eta_j)
// so it is ordered according to argsort(T)

void _update_cox_exp(double *linear_pred_ptr, /* Linear term in objective */
		     double *exp_accum_ptr,   /* inner accumulation vector */
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
  double *exp_accum, *outer_accum;
  double cur_val = 0;

  // reversed reverse cumsum of exp(eta)

  for (idx=0; idx<ncase; idx++) {
    order_idx = *((long *) ordering_ptr + (ncase - 1 - idx));
    linear_pred = *((double *) linear_pred_ptr + order_idx);
    cur_val = cur_val + exp(linear_pred);
    exp_accum = ((double *) exp_accum_ptr + (ncase - 1 - idx));
    *exp_accum = cur_val;
  }

}

void _update_cox_expZ(double *linear_pred_ptr,  /* Linear term in objective */
		      double *right_vector_ptr, /* Linear term in objective */
		      double *expZ_accum_ptr,   /* inner accumulation vector */
		      long *censoring_ptr,      /* censoring indicator */
		      long *ordering_ptr,       /* 0-based ordering of times */
		      long *rankmin_ptr,        /* 0-based ranking with min tie breaking */
		      long ncase                /* how many subjects / times */
		      )       
{
  long idx;
  long order_idx, rankmin_idx;
  double linear_pred, right_vector;
  long censoring;
  double *expZ_accum, *outer_accum;
  double cur_val = 0;

  // reversed reverse cumsum of exp(eta)

  for (idx=0; idx<ncase; idx++) {
    order_idx = *((long *) ordering_ptr + (ncase - 1 - idx));
    linear_pred = *((double *) linear_pred_ptr + order_idx);
    right_vector = *((double *) right_vector_ptr + order_idx);
    cur_val = cur_val + right_vector * exp(linear_pred);
    expZ_accum = ((double *) expZ_accum_ptr + (ncase - 1 - idx));
    *expZ_accum = cur_val;
  }

}

void _update_outer_1st(double *linear_pred_ptr,     /* Linear term in objective */
		       double *exp_accum_ptr,       /* inner accumulation vector */
		       double *outer_1st_accum_ptr, /* outer accumulation vector */
		       long *censoring_ptr,         /* censoring indicator */
		       long *ordering_ptr,          /* 0-based ordering of times */
		       long *rankmin_ptr,           /* 0-based ranking with min tie breaking */
		       long ncase                   /* how many subjects / times */
		       )       
{
  long idx;
  long order_idx, rankmin_idx;
  double linear_pred;
  long censoring;
  double *exp_accum, *outer_1st_accum;
  double cur_val = 0;

  // accumulate inverse cumsums at rankmin
  // i-th value is \sum_{j=1}^i 1 / W(r[o[j]])  where r is rankmin of times so r[o]
  // is rankmin of ordered times

  cur_val = 0;
  for (idx=0; idx<ncase; idx++) {
    order_idx = *((long *) ordering_ptr + idx);
    rankmin_idx = *((long *) rankmin_ptr + order_idx);
    exp_accum = ((double *) exp_accum_ptr + rankmin_idx);
    censoring = *((long *) censoring_ptr + order_idx);
    cur_val = cur_val + censoring / *exp_accum;
    outer_1st_accum = ((double *) outer_1st_accum_ptr + idx);
    *outer_1st_accum = cur_val;
  }

}

void _update_outer_2nd(double *linear_pred_ptr,     /* Linear term in objective */
		       double *exp_accum_ptr,       /* inner accumulation vector e^{\eta} */
		       double *expZ_accum_ptr,      /* inner accumulation vector  Ze^{\eta} */
		       double *outer_2nd_accum_ptr, /* outer accumulation vector */
		       long *censoring_ptr,         /* censoring indicator */
		       long *ordering_ptr,          /* 0-based ordering of times */
		       long *rankmin_ptr,           /* 0-based ranking with min tie breaking */
		       long ncase                   /* how many subjects / times */
		       )       
{
  long idx;
  long order_idx, rankmin_idx;
  double linear_pred, right_vector;
  long censoring;
  double *expZ_accum, *exp_accum, *outer_2nd_accum;
  double cur_val_num = 0;
  double cur_val_den = 0;
  double cur_val = 0;

  // accumulate inverse cumsums at rankmin
  // i-th value is \sum_{j=1}^i 1 / W(r[o[j]])  where r is rankmin of times so r[o]
  // is rankmin of ordered times

  for (idx=0; idx<ncase; idx++) {
    order_idx = *((long *) ordering_ptr + idx);
    rankmin_idx = *((long *) rankmin_ptr + order_idx);
    expZ_accum = ((double *) expZ_accum_ptr + rankmin_idx);
    exp_accum = ((double *) exp_accum_ptr + rankmin_idx);
    censoring = *((long *) censoring_ptr + order_idx);
    cur_val = cur_val + censoring * (*expZ_accum) / ((*exp_accum) * (*exp_accum));
    outer_2nd_accum = ((double *) outer_2nd_accum_ptr + idx);
    *outer_2nd_accum = cur_val;
  }

}

// Objective value

double _cox_objective(double *linear_pred_ptr,     /* Linear term in objective */
		      double *inner_accum_ptr,     /* inner accumulation vector */
		      double *outer_1st_accum_ptr, /* outer accumulation vector */
		      long *censoring_ptr,         /* censoring indicator */
		      long *ordering_ptr,          /* 0-based ordering of times */
		      long *rankmin_ptr,           /* 0-based ranking with min tie breaking */
		      long *rankmax_ptr,           /* 0-based ranking with max tie breaking */
		      long ncase                   /* how many subjects / times */
		      )       
{
  long idx, rankmin_idx;
  double linear_pred, inner_accum;
  long censoring;
  double cur_val = 0;

  // ensure you have updated the inner / outer accumulation
  // vectors with current linear predictors
  // this can be done in the wrapper _update_cox_weights

  for (idx=0; idx<ncase; idx++) {
    rankmin_idx = *((long *) rankmin_ptr + idx);
    inner_accum = *((double *) inner_accum_ptr + rankmin_idx);
    censoring = *((long *) censoring_ptr + idx);
    linear_pred = *((double *) linear_pred_ptr + idx);
    cur_val += (censoring) * (log(inner_accum) - linear_pred);
  }

  return(cur_val);

}

void _cox_gradient(double *gradient_ptr,        /* Where gradient is stored */
		   double *linear_pred_ptr,     /* Linear term in objective */
		   double *outer_1st_accum_ptr, /* outer accumulation vector */
		   long *censoring_ptr,         /* censoring indicator */
		   long *ordering_ptr,          /* 0-based ordering of times */
		   long *rankmin_ptr,           /* 0-based ranking with min tie breaking */
		   long *rankmax_ptr,           /* 0-based ranking with max tie breaking */
		   long ncase                   /* how many subjects / times */
		   )
{
  long idx, rankmax_idx;
  double linear_pred, outer_1st_accum;
  double *gradient;
  long censoring;

  // ensure you have updated the inner / outer accumulation
  // vectors with current linear predictors
  // this can be done in the wrapper _update_cox_weights

  // fill in entries of gradient

  for (idx=0; idx<ncase; idx++) {
    censoring = *((long *) censoring_ptr + idx);
    rankmax_idx = *((long *) rankmax_ptr + idx);
    outer_1st_accum = *((double *) outer_1st_accum_ptr + rankmax_idx);
    linear_pred = *((double *) linear_pred_ptr + idx);
    gradient = ((double *) gradient_ptr + idx);
    *gradient = outer_1st_accum * exp(linear_pred) - censoring;
  }

}

void _cox_hessian(double *hessian_ptr,          /* Where hessian is stored */
		  double *linear_pred_ptr,      /* Linear term in objective */
		  double *outer_1st_accum_ptr,  /* outer accumulation vector used in outer prod "mean"*/
		  double *outer_2nd_accum_ptr,  /* outer accumulation vector used in "2nd" moment*/
		  long *censoring_ptr,          /* censoring indicator */
		  long *ordering_ptr,           /* 0-based ordering of times */
		  long *rankmax_ptr,            /* 0-based ranking with max tie breaking */
		  long ncase                    /* how many subjects / times */
		  )
{
  long idx, rankmax_idx;
  double linear_pred, outer_1st_accum, outer_2nd_accum;
  double *hessian;
  long censoring;

  // ensure you have updated the inner / outer accumulation
  // vectors with current linear predictors
  // this can be done in the wrapper with _update_cox_hessian

  // fill in entries of hessian

  for (idx=0; idx<ncase; idx++) {
    censoring = *((long *) censoring_ptr + idx);
    rankmax_idx = *((long *) rankmax_ptr + idx);
    outer_1st_accum = *((double *) outer_1st_accum_ptr + rankmax_idx);
    outer_2nd_accum = *((double *) outer_2nd_accum_ptr + rankmax_idx);
    linear_pred = *((double *) linear_pred_ptr + idx);
    hessian = ((double *) hessian_ptr + idx);
    *hessian =  exp(linear_pred) * (outer_1st_accum - outer_2nd_accum);
  }

}

