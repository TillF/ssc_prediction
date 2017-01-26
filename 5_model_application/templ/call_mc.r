# load models, ancillary data and config-file and call MC-routine (apply_model_mc)
  
  source("../../settings.R")
  library(quantregForest)
  if (installed.packages()["quantregForest","Version"]  < "0.3-5") stop("quantregForest Version >= 0.3-5 required")

  source("../apply_model_mc_ssc_q.R")

  for (target_var in c("ssc", "discharge")) #load models for all target variables
  {
    model_file = paste0("../../5_model_building/model_", target_var,".RData")
    multi_model_file = sub(model_file, pattern = "\\.RData", repl="_multi.RData")
    model_name = paste0("qrf_model_", target_var)  #name for model of target variable
    
    if (file.exists(multi_model_file))   #try to load multi-model file first
    {
      load(file =multi_model_file)       #load multi discharge-model from file (must be a list of models named "model_variants" or adjust name of model in next line)
      assign(x = model_name, value = model_variants)
      rm(model_variants)
	  print(paste0(multi_model_file," loaded."))
    } else
    if (file.exists(model_file))
    {
      load(file = model_file)       #load discharge-model from file (must be named "best_model_discharge" or adjust name of model in next line)
      assign(x = model_name, value = list(best_model))
      rm(best_model)
	  print(paste0(multi_model_file," loaded."))
    } else
	{
		assign(x = model_name, value = NULL) #no model for discharge available
		print(paste0(target_var,": no model loaded."))
	}	
  }
    
    
  if (file.exists("conf.txt"))
  { #read config file
    conf=read.table("conf.txt", header=TRUE)
    seed = conf$seed
    flood_period=eval(parse(text=as.character(conf$flood_period)))  #allow for specification of multiple period in the style of "3:5"
  } else
  {
    seed = 1
    flood_period=1
  }


  set.seed(seed) #for making results reproducable independent of MC-randomness
  ##call MC-model
  xx = apply_model_mc(gauge_name=gauge_name, nrealisations=nrealisations, q_conf=q_conf, subset_periods=c(flood_period),individual=FALSE, verbose=verbose,
      predictordata=paste("../../3_predictor_generation/ancillary_data_",tres,".RData",sep=""), iqr=iqr, write_dist=write_dist,overwrite=overwrite, do_interfloods=do_interfloods)


