# load models, ancillary data and config-file and call MC-routine (apply_model_mc)

  source("../../settings.R")
  library(quantregForest)
  if (installed.packages()["quantregForest","Version"]  < "0.3-5") stop("quantregForest Version >= 0.3-5 required")

  source("../apply_model_mc_ssc_q.R")

  load(file  ="../../5_model_building/model_ssc.RData")       #load ssc-model from file (must be named "best_model_ssc" or adjust name of model in next line)
  qrf_model_ssc  =best_model
  rm(best_model)

  discharge_model_file = "../../5_model_building/model_discharge.RData"
  discharge_multi_model_file = sub(discharge_model_file, pattern = "\\.RData", repl="_multi.Rdata")
  
  if (file.exists(discharge_multi_model_file))   #try to load multi-model file first
  {
    load(file =discharge_multi_model_file)       #load mulit discharge-model from file (must be a list of models named "model_variants" or adjust name of model in next line)
    qrf_model_discharge  =model_variants
    rm(model_variants)
  } else
  if (file.exists(discharge_model_file))
  {
    load(file =discharge_model_file)       #load discharge-model from file (must be named "best_model_discharge" or adjust name of model in next line)
    qrf_model_discharge  =list(best_model)
    rm(best_model)
  } else
  qrf_model_discharge = NULL               #no model for discharge available

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


