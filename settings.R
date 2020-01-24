##user settings

##input
	tres= as.difftime(15,units="mins") ##temporal resolution of input files  (Q, P) and desired output [min]


	precip_data_file = "0_data\\TESP_15minres_2010_2011.txt"    #input file containing (mutiple) rainfall series
	fmt_precip = "d.%m.%Y %H:%M" #input format of date column, eg "%Y-%m-%d %H:%M" or "%d.%m.%Y %H:%M"

	discharge_data_file = "0_data\\H1_discharge_puebla.txt" #input file containing discahrge time series (single gauge)
	fmt_discharge = "%d.%m.%Y %H:%M" #input format of date column

	ssc_data_file    = "0_data\\H1_ssc_puebla.txt" #input file containing intermittent SSC samples (single gauge)
	fmt_ssc = "%d.%m.%Y %H:%M"  #input format of date column

  gauge_name="h1" #name of discharge gauge
  
##4_model_building
  mtrybest = "tune" #RF-parameter, set to "tune" to auto-adjust or a fixed value (1 to number of predictors)
  doplot = FALSE    #enable/disable validation plots
  importance = TRUE #estimate (and plot) variable importance (only for RF)
  dontuse_cols = c("datenum", "anothercolumn")  #specify column names that should not be used as predictors
  do_continous_validation = FALSE #do "continuous" validation with moving blocks of test fractions (more suitable for discharge prediction)

  n_valid_periods = 5  ##number of validation periods  ( n-fold cross-validation by excluding every n-th part from the training data)
  n_perturbations = NULL #number of realisations of perturbed training data to create for target variables ()
  n_perturbations["ssc"]       = 1 
  n_perturbations["discharge"] = 1 

##5_model_application
  n_slaves = 1 #number of slaves available (that many replicate directories will be generated)
  flood_periods = 1:16 ##selection of flood periods that will be treated (according to flood_numbering.txt)
  do_interfloods = FALSE #should the periods between the floods be treated, too?
  nrealisations=250 #number of MC-realisations to make for each period
  q_conf=1.00 #quantile range used for drawing the random quantile
  verbose=TRUE #enable progress messages
  iqr=TRUE # optional output of IQR-values of all realisations of each timestep
  write_dist=FALSE # optional output of entire sampled distribution of each timestep (call with write_dist=TRUE)
  overwrite=TRUE #overwrite (TRUE) or resume (FALSE) previous runs