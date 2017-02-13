build_validate_model=function(response_variable, tres, model_type)
{
#  response_variable="ssc"         #name of response variable
#  tres = 5                        #temporal resolution of time series [min]
#  model_type="RF"                 # model type: RF or QRF

set.seed(42)                        #enable reproducability
ntree = 500
doplot=FALSE

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
  load(file=paste("../3_predictor_generation/ancillary_data_",tres,"_train_",response_variable,".RData",sep=""))     #load training data
  print("predictor data loaded"); flush.console()
  response_var_index=which(names(mydata_training)==response_variable)  #get index to ssc or discharge-column
  dontuse=c(response_var_index, which(names(mydata_training)=="datenum"))  #get indeces to columns not to use in regression (i.e. date and response variable)

  #tune mtry-parameter for randomForest
#  windows()
#	tune_r=tuneRF(x=mydata_training[,-dontuse],y=mydata_training[,response_var_index], ntreeTry=1000, stepFactor=1.5, improve = 0, trace=TRUE, plot=T)
#  mtrybest=tune_r[which.min(tune_r[,2]),1]
   mtrybest=8
  # Hastie et al., 2009, p.592: "...In practice the best value for these parameters will depend on the problem,
  # and they should be treated as tuning parameters.

  print("building model (full training period)"); flush.console()
  if (model_type=="RF")
  rf_model = randomForest(mydata_training[,response_var_index] ~ ., data = mydata_training[,-dontuse], importance = TRUE, mtry=mtrybest, ntree = ntree) else #build actual model
	rf_model=  quantregForest(mydata_training[,-dontuse], mydata_training[,response_var_index], mtry = mtrybest, ntree = ntree, nodesize=5, keep.inbag=TRUE)

  print("saving model"); flush.console()
  best_model=rf_model # rename before saving
  if (model_type!="RF") #only QRF model can be used in MC-routine
    save(file=paste("model_",response_variable,".RData",sep=""),list="best_model") # save model for later use


# plot model and its performance
  print("apply model (full training period)"); flush.console()
  if (model_type=="RF")
  prediction_OOB=predict(rf_model) else 
  if (model_type=="QRFm")		#mean of QRF quantile values (instead of median)
  {                            
	 tt=predict(rf_model,quantiles=seq(from=0,to=1,by=0.01) )  #predict 100 quantile values
	 prediction_OOB= apply(tt,1,mean)
  } 
  else
  prediction_OOB=predict(rf_model,quantiles=0.5) 
  nash_OOB_full = nashS(prediction_OOB, mydata_training[,response_var_index])        #Nash-Sutcliffe-index of OOB-predictions of full dataset
	if (doplot)
	{
    windows(height=5,width=5)
    plot(mydata_training[,response_var_index],prediction_OOB,xlab="measured",ylab="predicted",main=paste("OOB-predictions for full training dataset",response_variable))
  	abline(0,1)
  # identify(mydata_training[,response_var_index],prediction_OOB, mydata_training$datenum)     #identify records by mouse click
  }

# 3-year crosswise validation of model
  print("loading ancillary data"); flush.console()
  load(file=paste("../3_predictor_generation/ancillary_data_",tres,".RData",sep=""))     #load entire dataset

  obsolete_columns  = which(!(names(mydata) %in% c(dimnames(best_model$importance)[[1]],response_variable,"datenum")))
  mydata[,obsolete_columns]=NULL                  #discard unnecessary columns
  mydata = na.omit(mydata)                        #remove records that contain no-data value

  response_var_index_validation =which(names(mydata)==response_variable)  #get index to column containing response variable
  dontuse_validation=c(response_var_index_validation, which(names(mydata)=="datenum"))  #get indeces to columns not to use in regression (i.e. date and response variable)
  nash_year=array(0,3)     #for storing the Nash-Sutcliffe indices of threefold cross-validation  (2008-2010)

  for (test_year in 2008:2010)
  {
    print(paste("test year:",test_year));    flush.console()
    to_year   = test_year +1         #end year
    from_year = test_year            #start year of test period
    if (test_year == 2008)
     from_year = 2007 else  #for 2008, add rest of 2007
     from_year = test_year

    test_index =  (mydata_training$datenum >= ISOdatetime(from_year,1,1,0,0,0,tz="GMT")) & (mydata_training$datenum<=ISOdatetime(to_year,1,1,0,0,0,tz="GMT"))  #define test period
    mydata_train_x= mydata_training[!test_index ,]            #use 2 years from training data to train

    test_index =  (mydata$datenum >= ISOdatetime(from_year,1,1,0,0,0,tz="GMT")) & (mydata$datenum<=ISOdatetime(to_year,1,1,0,0,0,tz="GMT"))  #define test period
    mydata_test_x = mydata        [ test_index ,]            #use complement 1 year from entire dataset to validate

    max_records=Inf #maximum number of records to use in validation (may take too long otherwise) #rr
#    print("remove me!");flush.console()
    mydata_test_x=mydata_test_x[1:(min(max_records,nrow(mydata_test_x))),]

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

    nash_year[test_year-2007] = nashS(predictions, mydata_test_x[,response_var_index_validation]) #store Nash-Sutcliffe index

  	if (doplot)
  	{
      windows()
      plot(mydata_test_x[,response_var_index_validation],predictions,col="red",pch=".",xlab="observed",ylab="predicted",main=paste(response_variable,"/",model_type,"( test period",test_year,")\nNSE =",nash_year[test_year-2007]))
      abline(0,1)
    }
  }

  nash_vals=data.frame(nash_OOB_full,t(nash_year))
  names(nash_vals)=c("nash_OOB_full",paste("nash",2008:2010,sep=""))
  print(nash_vals)
}