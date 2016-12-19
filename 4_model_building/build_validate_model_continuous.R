build_validate_model_continuous=function(response_variable="ssc", tres=5, ntree = 1000, gauge_name="H1", percentages4test = c(0,5,10,15,20,25,30,40,50), #percentage of test data)
            dontuse_cols="datenum", model_type="RF", blocks2test=NULL)
{
#  response_variable="ssc"         #name of response variable
#  tres = 5                        #temporal resolution of time series [min]
#  model_type="RF"                 # model type: RF or QRF

set.seed(1)                        #enable reproducability

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
  load(file=paste("../2_predictor_generation/ancillary_data_",tres,"_train_",response_variable,".RData",sep=""))     #load training data
  print("predictor data loaded"); flush.console()
  mydata_training = mydata_training[,!(names(mydata_training) %in% dontuse_cols)]  #remove disabled columns
  print(paste(nrow(mydata_training),"rows contained in raw training data"))
  mydata_training = na.omit(mydata_training) #remove nodata cases
  print(paste(nrow(mydata_training),"rows remaining in training data after NA-removal"))  
  response_var_index=which(names(mydata_training)==response_variable)  #get index to column of predictor variable
  dontuse=c(response_var_index, which(names(mydata_training)=="datenum"))  #get indeces to columns not to use in regression (i.e. date and response variable)
  

####

  source("../5_model_application/apply_model_mc_ssc_q.r")

  predict_mc = function(imydata)      ##apply mc-routine on test-data set
  {
    flood_scheme=data.frame(no=1:nrow(imydata),begin=imydata$datenum,end=imydata$datenum) 
    predictions=apply_model_mc(gauge_name=gauge_name,nrealisations=100,q_conf=1.,flood_scheme=flood_scheme,individual=F,
        predictordata=paste("../2_predictor_generation/ancillary_data_",tres,".RData",sep=""),iqr=FALSE, 
        write_dist=FALSE, write_files=FALSE, qrf_model_discharge_arg=qrf_model_discharge, qrf_model_ssc_arg = qrf_model_ssc)
    return(predictions)
  }


  mydata_training=na.omit(mydata_training)
 # mydata_training=mydata_training[1:20,]; cat("testoption ausschalten")  
 print(paste(nrow(mydata_training),"total rows"))

 all_goodness=data.frame()       #create empty dataframe for collecting all validation results across multiple percentages
 
 if (is.null(blocks2test))
   block_start=1:nrow(mydata_training) else #test all possible blockstart
     block_start=round(seq(from=1, to=nrow(mydata_training), length.out=min(nrow(mydata_training), blocks2test))) #only test specified number of blockstarts
 
 
for (j in 1:length(percentages4test))
{
 # sink("log.txt")
  ntest= round(nrow(mydata_training) * percentages4test[j]/100)       #number of records to be included in test data
  print(paste(ntest,"rows used for testing"))
  goodness=data.frame(trainingNS=array(NA,nrow(mydata_training)),testNS=array(0,nrow(mydata_training)), trainingRMSE=array(0,nrow(mydata_training)),testRMSE=array(0,nrow(mydata_training)))
  
  for (i in block_start)
  {   
     print(paste("test set", i,"of",length(block_start), " (at ",j," of ",length(percentages4test), " percentages)"))
     if (ntest > 0)
     {
      testrecords=((i-1+(1:ntest))-1) %% nrow(mydata_training) + 1    #define test records
      mydata_training_selected=mydata_training[-testrecords,]         #assemble dataset for training
      mydata_test             =mydata_training[ testrecords,]         #assemble dataset for testing
     } else      {
      mydata_training_selected=mydata_training
      mydata_test             =NULL
     }
     
    
      #search for optimum mtry
    	a=tuneRF(x=mydata_training_selected[,-response_var_index],y=mydata_training_selected[,response_var_index], ntreeTry=1000, stepFactor=1.15, improve=0.002,trace=FALSE, plot=FALSE)
    	mtrybest=a[which.min(a[,2]),1]
      # Hastie et al., 2009, p.592: "...In practice the best value for these parameters will depend on the problem,
      # and they should be treated as tuning parameters.
    
    	if (model_type=="RF")
    	 { 
    	  rf_model <- randomForest(mydata_training_selected[,response_var_index] ~ ., data = mydata_training_selected[,-response_var_index], importance = FALSE, mtry=mtrybest, ntree=ntree)
    	  ssc_pred=predict(rf_model)      #do prediction on training data
  	  }
    	else
    	{  
    	  qrf_model_ssc=quantregForest(mydata_training_selected[,-response_var_index], mydata_training_selected[,response_var_index], mtry = mtrybest, ntree = ntree, nodesize=5)
        qrf_model_discharge = NULL
        ssc_pred=predict_mc(mydata_training_selected[,-response_var_index])      #do prediction on training data
      }
    	
      goodness$trainingNS[i]=nashS(ssc_pred,mydata_training_selected$ssc)
    	goodness$trainingRMSE[i]=sqrt(mean((ssc_pred-mydata_training_selected$ssc)^2,na.rm=TRUE))

#    	plot(ssc_pred,mydata_training_selected$ssc)
#      windows()
#      plot(predict(rf_model),mydata_training_selected$ssc)
     	if (is.null(mydata_test))
     	{
    	 goodness$trainingNS     =goodness$trainingNS[1]
     	 goodness$trainingRMSE[i]=goodness$trainingRMSE[1]
    	 break
      }else
      {
        #do prediction on test data
        if (model_type=="RF")
          ssc_pred=predict(rf_model, newdata = mydata_test[,-response_var_index])      
        else
          ssc_pred=predict_mc(mydata_test[,-response_var_index])               
        goodness$testNS[i]    =nashS(ssc_pred,mydata_test$ssc)     	
       	goodness$testRMSE[i]  =sqrt(mean((ssc_pred-mydata_test$ssc)^2,na.rm=TRUE))

        print(paste("ns_ntraining =",goodness$trainingNS[i],"; ns_test =",goodness$testNS[i])); flush.console()
      }

  }
  
  all_goodness = na.omit(all_goodness)
  all_goodness=rbind(all_goodness, cbind(percentage4test=percentages4test[j], test_block_start=1:nrow(mydata_training), goodness))
  
  write.table(goodness, paste(percentages4test[j],"_validation_NS.txt",sep=""), sep="\t", row.names = F,quote=FALSE,col.names=TRUE) #write goodness measures for current percentage
  print(paste("mean(ns_ntraining) =",mean(goodness$trainingNS),"; mean(ns_test) =",mean(goodness$testNS)))
}

write.table(all_goodness, paste("validation_NS.txt",sep=""), sep="\t", row.names = F,quote=FALSE,col.names=TRUE) #write all goodness measures 

windows()
plot(all_goodness$test_block_start,all_goodness$trainingNS,xlab="index of start of test block [-]",ylab="NS", col=all_goodness$percentage4test, pch=20, ylim=c(quantile(all_goodness$testNS,0.1), 1))
points(all_goodness$test_block_start,all_goodness$testNS, col=all_goodness$percentage4test, pch=21)
legend("bottomright",legend=c("training","test"),pch=20:21)
savePlot(filename = "moving_Nash.emf", type = "emf") 
#sink()

}

