#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

double _cox_objective(double *linear_pred_ptr, /* Linear term in objective */
		      double *weight_ptr,      /* accumulated weights */
		      double *weight_denom_ptr,/* accumulated inverse weights */
		      long *censoring_ptr,     /* censoring indicator */
		      long *ordering_ptr,      /* 0-based ordering of times */
		      long *rankmin_ptr,       /* 0-based ranking with min tie breaking */
		      long *rankmax_ptr,       /* 0-based ranking with max tie breaking */
		      long ncase               /* how many subjects / times */
		      );       

void _cox_gradient(double *gradient_ptr,    /* Where gradient is stored */
		   double *linear_pred_ptr, /* Linear term in objective */
		   double *weight_ptr,      /* accumulated weights */
		   double *weight_denom_ptr,/* accumulated inverse weights */
		   long *censoring_ptr,     /* censoring indicator */
		   long *ordering_ptr,      /* 0-based ordering of times */
		   long *rankmin_ptr,       /* 0-based ranking with min tie breaking */
		   long *rankmax_ptr,       /* 0-based ranking with max tie breaking */
		   long ncase               /* how many subjects / times */
		   );

void _update_cox_weights(double *linear_pred_ptr, /* Linear term in objective */
			 double *weight_ptr,      /* accumulated weights */
			 double *weight_denom_ptr,/* accumulated denom terms */
			 long *censoring_ptr,     /* censoring indicator */
			 long *ordering_ptr,      /* 0-based ordering of times */
			 long *rankmin_ptr,       /* 0-based ranking with min tie breaking */
			 long ncase               /* how many subjects / times */
			 );       

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */
