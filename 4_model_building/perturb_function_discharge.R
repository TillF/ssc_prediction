# perturbation function for discharge
# perturbs the given observations with assumed errors from A. Zimmermann's observations of measured and recorded water stage
# perturbed discharge is a linear function of measured with coefficients derived from 5 different instruments (see error.R)

perturb_observations = function(mydata_training, replicate_no)
{
  #read specified functions representing discharge error
  perturbation_coefficients=read.table(header=TRUE, file="error_coefs.txt", sep="\t"	)                                      
  
  discharge_perturbed = mydata_training$discharge * perturbation_coefficients$slope[mod_no] + perturbation_coefficients$slope[mod_no] 
  
  if (any(discharge_perturbed < 0)) warning(paste0(sum(discharge_perturbed<0), " negative discharge produced due to perturbations. Values truncated to 0."))
  discharge_perturbed = pmax(discharge_perturbed, 0)  #discard negative values
  
  return(discharge_perturbed)
}