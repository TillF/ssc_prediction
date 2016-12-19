# build RF/QRF-model for SSC / discharge (based on training data sets generated in 2_predictor_generation)
# perform n-fold cross-validation

# clear workspace
	rm(list = ls())
	source("../settings.R")

# load necessary packages
	library(randomForest)
	library(quantregForest)
  source("build_validate_model.R")
  

#SSC
  graphics.off()

	#quick version without cross validition
	build_validate_model(response_variable="ssc", tres = tres, model_type="RF", importance=importance, dontuse_cols=dontuse_cols, n_valid_periods = 0, doplot=doplot, do_identify=TRUE)
	
  #replace NA values with finite number to check if model benefits from (potentially) additional training data
	build_validate_model(response_variable="ssc", tres = tres, model_type="RF", importance=importance, dontuse_cols=dontuse_cols, n_valid_periods = 0, doplot=doplot, replaceNAs=TRUE)
	
	#quick version without cross validition - for plot in paper
	build_validate_model(response_variable="ssc", tres = tres, model_type="QRF", importance=importance, dontuse_cols=dontuse_cols, n_valid_periods = 0, doplot=doplot, plotWhiskers=TRUE, OOB=FALSE)
	
  
  
  #build SSC model (RF)
  build_validate_model(response_variable="ssc", tres = tres, model_type="RF", importance=importance, dontuse_cols=dontuse_cols, n_valid_periods = n_valid_periods, doplot=doplot, replaceNAs=replaceNAs)

  #build SSC model (QRF, prediction from mean instead of median)
  build_validate_model(response_variable="ssc", tres = tres, model_type="QRFm", dontuse_cols=dontuse_cols, n_valid_periods = n_valid_periods, doplot=doplot, replaceNAs=replaceNAs)

  if (do_continous_validation)
  {
    source("build_validate_model_continuous.R")               
    build_validate_model_continuous(percentages4test = 10, model_type="RF", blocks2test = 10)
  }  


  cat("finished.\7\nPlease proceed with ../5_model_application (see readme.txt)")

