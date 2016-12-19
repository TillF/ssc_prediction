build_validate_model=function(response_variable, tres, model_type, importance = FALSE, ntree = 1000, mtrybest="tune", doplot=FALSE,
            dontuse_cols="datenum", n_valid_periods = 5, replaceNAs=FALSE, plotWhiskers=FALSE, OOB=TRUE, do_identify=FALSE)
{
#  response_variable="ssc"         #name of response variable
#  tres = 5                        #temporal resolution of time series [min]
#  model_type="RF"                 # model type: RF or QRF

set.seed(1)                        #enable reproducability


`nashS` <- function (modelled, measured, weight=NA) { #nash sutcliffe index
    t.na <- is.na(measured) | is.na(modelled)
    t.meas <- measured[!t.na]
    t.model <- modelled[!t.na]
    t.mean <- mean(t.meas)
    if(is.na(weight[1])){
        weight <- rep(1,NROW(t.meas))
    } else {
        weight <- weight[!t.na]
    }
    t.a <- (t.meas - t.model)*weight
    t.b <- (t.meas - t.mean)*weight
    ns <- 1- ( t.a %*% t.a / t.b %*% t.b)
    return(ns)
}



# load training data
  load(file=paste("../2_predictor_generation/ancillary_data_",tres,"_train_",response_variable,".RData",sep=""))     #load training data
  print("predictor data loaded"); flush.console()
  if (!identical(dontuse_cols,"datenum"))
    mydata_training = mydata_training[,!grepl(paste(setdiff(dontuse_cols,"datenum"),collapse="|"), names(mydata_training))]  #remove disabled columns
  
#   if (file.exists("outliers.RData"))
#   {
#     load("outliers.RData");
#     mydata_training=mydata_training[!(mydata_training$datenum %in% outliers),] #remove critical points
#   }  

  print(paste(nrow(mydata_training), "rows contained in raw training data"))
  print("nodata cases per column:")
  print(sort((apply(is.na(mydata_training), MARGIN=2, sum))))

  if (replaceNAs)
  {
  	print("replaceNAs enabled, replacing NAs by -999")
  	for (i in 1:ncol(mydata_training))
  	{
  		if (names(mydata_training)[i]=="datenum") next
  		mydata_training[is.na(mydata_training[,i]),i]=-999
	  }	
  }
  
  mydata_training = na.omit(mydata_training) #remove nodata cases
  print(paste(nrow(mydata_training), "rows remaining in training data after NA-removal"))
  response_var_index=which(toupper(names(mydata_training))==toupper(response_variable))  #get index to column of predictor variable
  names(mydata_training)[response_var_index]=response_variable #ensure correct spelling
  dontuse=c(response_var_index, which(names(mydata_training) %in% dontuse_cols))  #get indices to columns not to use in regression (i.e. usually date and response variable)
  if (mtrybest=="tune")   #tune mtry-parameter for randomForest
  {
    windows()
  	mtry_default=if (!is.null(mydata_training[,response_var_index]) && !is.factor(mydata_training[,response_var_index]))
             max(floor(ncol(mydata_training[,-dontuse])/3), 1) else floor(sqrt(ncol(mydata_training[,-dontuse]))) #default mtry of RF
    tune_r=tuneRF(x=mydata_training[,-dontuse],y=mydata_training[,response_var_index], ntreeTry=1000, stepFactor=mtry_default/(mtry_default-max(1,mtry_default/4)), improve = 0, trace=TRUE, plot=doplot)
    mtrybest=tune_r[which.min(tune_r[,2]),1]
  # Hastie et al., 2009, p.592: "...In practice the best value for these parameters will depend on the problem,
  # and they should be treated as tuning parameters.
  }


  print("building model (full training period)"); flush.console()
  
  if (model_type=="RF")
  {
    rf_model = randomForest(mydata_training[,response_var_index] ~ ., data = mydata_training[,-dontuse], importance = importance, mtry=mtrybest, ntree = ntree) 
    print(round(importance(rf_model), 2)[ sort.int(importance(rf_model)[,1], index.return=TRUE)$ix,])

    if (doplot)
    {
      windows()
      varImpPlot(rf_model, sort=TRUE, type=1)
      savePlot(filename=paste("../saved_plots/",gauge_name,"_importance.wmf",sep=""),type="wmf")
    }  
  }   else #build actual model
	rf_model=  quantregForest(mydata_training[,-dontuse], mydata_training[,response_var_index], mtry = mtrybest, ntree = ntree, nodesize=5)

  print("saving model"); flush.console()
  best_model=rf_model # rename before saving
  if (model_type!="RF") #only QRF model can be used in MC-routine
    save(file=paste("model_",response_variable,".RData",sep=""),list="best_model") # save model for later use


# plot model and its performance
  print("apply model (full training period)"); flush.console()
  if (OOB)
  {
    newdata=NULL 
    print("OOB-data")
  }  else
  {  newdata=mydata_training[,-dontuse]
     print("all data (no OOB)")
  }

  if (model_type=="RF")
  {  
    if (is.null(newdata))
      
      prediction_OOB=predict(rf_model) else
      prediction_OOB=predict(rf_model, newdata=newdata)
  } else
  {
    if (model_type=="QRFm")		#mean of QRF quantile values (instead of median)
    {                            
  	 tt=predict(rf_model,quantiles=seq(from=0,to=1,by=1/250) )  #predict 250 quantile values
  	 prediction_OOB= apply(tt,1,mean)
    }   else
    {  
      #prediction_OOB=predict(rf_model,quantiles=0.5) 
      tt=predict(rf_model,quantiles=c(0.05,0.5, 0.95), newdata=newdata) #also predict quantiles
      prediction_OOB=tt[,2]
      tt=tt[,-2]
    }  
  }  
  nash_OOB_full = nashS(prediction_OOB, mydata_training[,response_var_index])        #Nash-Sutcliffe-index of OOB-predictions of full dataset
  rmse_OOB_full = sqrt(mean((prediction_OOB - mydata_training[,response_var_index])^2)) #store rmse
  print(paste("NS_oob:",nash_OOB_full,"; RMSE_OOB:", rmse_OOB_full))
	if (doplot)
	{
    windows(height=5,width=5)
    maxssc= max(max(mydata_training[,response_var_index]), max(prediction_OOB))
    plot(mydata_training[,response_var_index],prediction_OOB,xlab="measured [g/l]",ylab="predicted [g/l]",main=paste("OOB-predictions for full training dataset",response_variable),
         xlim=c(0, maxssc),     ylim=c(0, maxssc))
  	if (plotWhiskers & exists("tt"))
  	{
      for (i in 1: nrow(tt))
        lines(rep(mydata_training[,response_var_index][i],2),tt[i,], col="grey")
      points(mydata_training[,response_var_index],prediction_OOB)
  	  legend("bottomright", legend=c("90 % CI", "1:1"), lty=c("solid", "dotted"), col=c("grey", "black"), bg="white")
  	}  
    abline(0,1, lty="dotted")
    
    if (do_identify)
    {  
      print("Please click points in plot to identify. Rightclick to finish.")
      id=identify(mydata_training[,response_var_index],prediction_OOB, 1:nrow(mydata_training))     #identify records by mouse click
      print(paste("identified point numbers:",paste(id, collapse=", ")))
      print(mydata_training$datenum[id])
      outliers=mydata_training$datenum[id]
      save(file="outliers.RData", list="outliers")
    }  
  }
#  stop()

	if (n_valid_periods==0) return()
# crosswise validation of model
  print("loading ancillary data"); flush.console()
  load(file=paste("../2_predictor_generation/ancillary_data_",tres,".RData",sep=""))     #load entire dataset

  obsolete_columns  = which(!(names(mydata) %in% c(dimnames(best_model$importance)[[1]],response_variable,"datenum")))
  for (i in rev(obsolete_columns))
    mydata[, i] = NULL                  #discard unnecessary columns

  if (replaceNAs)
  {
    print("replaceNAs enabled, replacing NAs by -999")
    for (i in 1:ncol(mydata))
    {
      if (names(mydata)[i]=="datenum") next
      mydata[is.na(mydata[,i]),i]=-999
    }	
  } else

  mydata = na.omit(mydata)                        #remove records that contain no-data values

  if (!(response_variable %in% names(mydata))) #add column of response var, if not already there
   mydata = merge(mydata, mydata_training[,c("datenum",response_variable )], by="datenum")         #


  
  response_var_index_validation =which(names(mydata)==response_variable)  #get index to column containing response variable
  dontuse_validation=c(response_var_index_validation, which(names(mydata)=="datenum"))  #get indeces to columns not to use in regression (i.e. date and response variable)
  nash_period=array(NA,n_valid_periods)     #for storing the Nash-Sutcliffe indices of n-fold cross-validation 
  rmse_period=array(NA,n_valid_periods)     #for storing the RMSE values of n-fold cross-validation


  test_period_dates=data.frame()

  maxssc= max(mydata[,response_var_index_validation], na.rm=TRUE) #for scaling plots

  for (test_period in 1:n_valid_periods)
  {
    print(paste("test period:",test_period));    flush.console()
  
    #test data
    test_index_from =  1 + round((test_period - 1)/n_valid_periods * nrow(mydata)) #index of start of test period
    test_index_to   =  round((test_period    )/n_valid_periods * nrow(mydata)) #index of end of test period
    test_index      =  array(FALSE, nrow(mydata)) #logical index vector
    test_index[test_index_from:test_index_to] = TRUE   #define test period
    mydata_test_x  = mydata[ test_index ,]            
    test_period_dates = rbind(test_period_dates, data.frame(from=mydata$datenum[test_index_from], to=mydata$datenum[test_index_to]))
    
    
    #training data
    test_index_from =  round((test_period - 1)/n_valid_periods * nrow(mydata_training)) #index of start of test period
    test_index_to   =  round((test_period    )/n_valid_periods * nrow(mydata_training)) #index of end of test period
    test_index      =  array(FALSE, nrow(mydata_training)) #logical index vector
    test_index[test_index_from:test_index_to] = TRUE   #define test period
    mydata_train_x = mydata_training[!test_index ,]            #use non-test data to train model


    if (nrow(mydata_train_x) == 0 | nrow(mydata_test_x) == 0)
      next                                                    #skip this year if traing or test dataset is empty

    max_records=Inf #maximum number of records to use in validation (may take too long otherwise) #rr
#    print("remove me!");flush.console()
    mydata_test_x = mydata_test_x[1:(min(max_records,nrow(mydata_test_x))),]

    #build model
    print("build model...");    flush.console()

    if (model_type=="RF")
    rf_model = randomForest(mydata_train_x[,response_var_index] ~ ., data = mydata_train_x[,-dontuse], mtry=mtrybest, ntree = ntree) else #build cross validation model
  	rf_model=  quantregForest(mydata_train_x[,-dontuse], mydata_train_x[,response_var_index], mtry = mtrybest, ntree = ntree, nodesize=5)

    print("apply model to test period...");    flush.console()
    stepsize=40000       #number of records that are predicted in one loop (reduce this number if you get memory allocation errors)
    predictions=array(NA,nrow(mydata_test_x))                #pre-allocate array
    for (i in seq(from=1, to=nrow(mydata_test_x),by=stepsize))      #do loops because otherwise memory may become short
    {
      print(paste(i,"/",nrow(mydata_test_x)))
      flush.console()
      if (model_type=="QRFm")
      {                            #mean of QRF quantile values (instead of median)
         tt=predict(rf_model,newdata=mydata_test_x[i:(min(i+stepsize,nrow(mydata_test_x))),-dontuse_validation],quantiles=seq(from=0,to=1,by=0.01) )  #predict 100 quantile values
         predictions[i:(min(i+stepsize,nrow(mydata_test_x)))]= apply(tt,1,mean)
      }     else
      predictions[i:(min(i+stepsize,nrow(mydata_test_x)))]=predict(rf_model,newdata=mydata_test_x[i:(min(i+stepsize,nrow(mydata_test_x))),-dontuse_validation],quantiles=0.5)
    
    }

    nash_period[test_period] = nashS(predictions, mydata_test_x[,response_var_index_validation]) #store Nash-Sutcliffe index
    rmse_period[test_period] = sqrt(mean((predictions - mydata_test_x[,response_var_index_validation])^2)) #store rmse

  	if (doplot)
  	{
      windows()
      #maxssc= max(max(mydata_test_x[,response_var_index_validation]), max(predictions))
      plot(mydata_test_x[,response_var_index_validation],predictions,col="red",pch=20,xlab="observed",ylab="predicted",main=paste(response_variable,"/",model_type,"( test period",test_period,";", 
              paste(format(test_period_dates[test_period,],"%d.%m.%Y"),collapse="-"),")\nNSE =",round(nash_period[test_period],3),
              "\nn_train/n_test = ",nrow(mydata_train_x),"/",nrow(mydata_test_x)),     xlim=c(0, maxssc),     ylim=c(0, maxssc))
      
      
      abline(0,1)
    }
  }

  performance_vals=data.frame(nash_OOB_full,t(nash_period), mean(nash_period),rmse_OOB_full,t(rmse_period), mean(rmse_period))
  names(performance_vals)=c("nash_OOB_full",paste("nash",1:n_valid_periods,sep=""), "meanNS","rmse_OOB_full",paste("rmse",1:n_valid_periods,sep=""), "meanRMSE")
  print(performance_vals)
  performance_vals$gauge_name=gauge_name
  performance_vals$model_type=model_type
  write.table(file="performance.txt",x=performance_vals,col.names=!file.exists("performance.txt"), append=TRUE, quote=FALSE, sep="\t", row.names=FALSE)
}

