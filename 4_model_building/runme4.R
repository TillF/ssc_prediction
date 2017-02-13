# build RF/QRF-model for SSC / discharge (based on training data sets generated in 2_predictor_generation)
# perform n-fold cross-validation

# for discharge and SSC data, disturb observations based on provided uncertainty in discharge and set up replicate models

#application_time_scale = "daily"

# clear workspace
	rm(list = ls())
	source("../settings.R")

# load necessary packages
	library(randomForest)
	library(quantregForest)
  source("build_validate_model2.R")
  

#SSC
  graphics.off()

	#quick version without cross validition
	#build_validate_model(response_variable="ssc", tres = tres, model_type="QRF", importance=importance, dontuse_cols=dontuse_cols, n_valid_periods = 10, doplot=TRUE, do_identify=FALSE, mtrybest=mtrybest)
	#build_validate_model(response_variable="ssc", tres = tres, model_type="RF", importance=importance, dontuse_cols=dontuse_cols, n_valid_periods = 10, doplot=TRUE, do_identify=FALSE, mtrybest=mtrybest)
	#build_validate_model(response_variable="ssc", tres = tres, model_type="QRFm", importance=importance, dontuse_cols=dontuse_cols, n_valid_periods = 10, doplot=TRUE, do_identify=FALSE, mtrybest=mtrybest)
	
  #permute training data
  source("build_validate_model.R")  #version with yearly cross validation
  
  for (target_var in c("ssc", "discharge"))
  {
    # load training data
    org_data_file = paste0("../3_predictor_generation/ancillary_data_",tres,"_train_", target_var,".RData") 
    load(org_data_file)     #load training data
    print("predictor data loaded"); flush.console()
    mydata_training_org = mydata_training
    file.rename(from = org_data_file, paste0(org_data_file,"_org"))  #keep the original data
    #read specified functions representing error
    source(paste0("perturb_function_", target_var,".R")) #include function for generating the error realisations
  
    model_variants=list() #collect one model per error realisation
    
    #loop thru all specified realisations of uncertainty
    for (mod_no in 1: n_perturbations[target_var]) 
    {
      mydata_training = mydata_training_org  #start from original data
      
      #perturb target variable in training data
      mydata_training[, target_var] = perturb_observations(mydata_training, replicate_no = mod_no)
      
      #store modified training data
      save(file=org_data_file, list = "mydata_training")
      
      #call actual model builing
      build_validate_model(response_variable=target_var, tres = tres, model_type="QRFm")
      
      load(paste0("model_", target_var, ".RData")) #load generated model ("best_model")
      model_variants[[mod_no]] = best_model #collect generated model
    }
    
    save(file=paste0("model_", target_var, "_multi.RData"), list="model_variants") #save all model variants in single file
  
    file.rename(to = org_data_file, from=paste0(org_data_file,"_org"))  #restore the original data
  }
  

  cat("finished.\7\nPlease proceed with ../5_model_application (see readme.txt)")
  