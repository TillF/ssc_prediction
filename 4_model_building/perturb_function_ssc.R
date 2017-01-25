# perturbation function for SSC
# perturbs the given observations with assumed errors from USDA study
# errors are normal, with their standard deviation being a linear function of the measured value

perturb_observations = function(mydata_training, replicate_no)
{
  set.seed(replicate_no)
  sd_modelled = 1/1000 * (0.0175* mydata_training$ssc*1000 +2.9611) #apply regression derived by A. Zimmermann from USDA study, convert from mg to g
  ssc_perturbed = rnorm(n = nrow(mydata_training), mean = mydata_training$ssc, sd = sd_modelled)
  #plot(mydata_training$ssc, ssc_perturbed)
  
  if (any(ssc_perturbed < 0)) warning(paste0(sum(ssc_perturbed<0), " negative SSC produced due to perturbations. Values truncated to 0."))
  ssc_perturbed = pmax(ssc_perturbed, 0)  #discard negative values
  
  return(ssc_perturbed)
}