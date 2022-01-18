#compute ancillary predictors from primary predictors stored in global variable mydata according to settings in anc_predictors_settings
#compute_ancillary_predictors=function(anc_predictors_settings, ancillary_data_filename)
{
    for (i in 1:nrow(anc_predictors_settings))              #create backup-columns of primary predictors
    {
      col_name=anc_predictors_settings$name[i]
      mydata[,paste(col_name,"_org",sep="")] =  mydata[,col_name]
    }

    datenum=mydata$datenum  #save date column and convert data frame to matrix (speed up computations)
    mydata$datenum=NULL
    mydata=as.matrix(mydata)

    nlag=1                #shift by single timestep  (to match Alex' scheme: ancillary predictor of level 0 is just shifted)
    anc_predictors_settings$n_aggregation_levels = anc_predictors_settings$n_aggregation_levels-1 #level 0 aggregation is already done with the next step
    col_names = anc_predictors_settings$name
    mydata[(nlag+1):nrow(mydata), col_names]          = mydata[1:(nrow(mydata)-nlag),col_names]

    n_datacolumns= ncol(mydata)          #number of data columns in dataset
    datacolumns4derivatives = which(dimnames(mydata)[[2]] %in% anc_predictors_settings$name)    ##which datacolumns will be used to construct derived predictors (use those mentioned in anc_predictors_settings)
    n_datacolumns4derivatives = length(datacolumns4derivatives)               #number of datacolumns used for building ancillary predictors

    mydata=cbind(mydata, array(NA,c(nrow(mydata), sum(anc_predictors_settings$n_aggregation_levels) ))) #enlarge matrix

    for (i in 1:nrow(anc_predictors_settings))         #loop over all primary predictors
    {
      source_column = which(dimnames(mydata)[[2]] == anc_predictors_settings$name[i])       #column index of current primary predictor
      if (length(source_column)==0) stop(paste("column",anc_predictors_settings$name[i],"not found"))
      p_prev  =-1
      cum_period_real=1              #number of timesteps to be included in cumulated sum
      cum_start = 0                  #start of aggregation window with respect to current timestep
      print(paste("primary predictor",i,"of",nrow(anc_predictors_settings),":",anc_predictors_settings$name[i]))

      for (aggregation_step in 1:anc_predictors_settings$n_aggregation_levels[i])          #loop through all data aggregation levels to construct ancillary predictors
      {
        p=aggregation_step/max(anc_predictors_settings$n_aggregation_levels[i])*100
        if (p>p_prev+ 0.5)   #do printout of percentage in 0.5 % steps (progress output)
        {
          cat(paste(p,"%\n"))
          flush.console()
          p_prev=p
        }
  
        cum_period_real=cum_period_real*anc_predictors_settings$aggregation_steps_factor[i]              #number of timesteps to be included in cumulated sum (real number)
        cum_period = round(cum_period_real)       #round real number to define aggregation window
        cum_start = cum_start - cum_period
  
        target_columns = n_datacolumns + sum(anc_predictors_settings$n_aggregation_levels[1:i-1]) +
                             +aggregation_step     #destination column where to put the aggregated data
        dimnames(mydata)[[2]][target_columns] = paste(anc_predictors_settings$name[i],cum_period*tres,sep="_")
  
        if (-cum_start+1 > nrow(mydata)) next
        if(cum_period==1) #no actual cumulation, just offset by cumperiod
          mydata[(cum_period+1):nrow(mydata),target_columns]          = mydata   [1:(nrow(mydata)-cum_period), source_column] else
        for (j in ((-cum_start+1):nrow(mydata)))      #actual aggregation
          mydata[j,target_columns] = apply(as.matrix(mydata[j + (cum_start : (cum_start+cum_period-1)), source_column],ncol=1), 2 ,sum,na.rm=FALSE)
  
      }

      dimnames(mydata)[[2]][dimnames(mydata)[[2]] == anc_predictors_settings$name[i]                     ] = 
            paste(anc_predictors_settings$name[i],tres,sep="_")              #denote primary predictor as first shifted time series
      dimnames(mydata)[[2]][dimnames(mydata)[[2]] == paste(anc_predictors_settings$name[i],"_org",sep="")] = 
            anc_predictors_settings$name[i]                                #restore column name of primary predictor

    }

    #compute rate-of-change in discharge (if discharge is among primary predictors)
    if ("discharge" %in% anc_predictors_settings$name)
    {
      print("computing rate of change in discharge...")
      time_lag=3;               #time to use for observing discharge before and after (in timesteps)
      mean_q_before  = filter(mydata[,"discharge"], filter=c(rep(0, time_lag), rep(1/(time_lag+1), time_lag+1)), sides = 2)
      mean_q_after   = c(mean_q_before[-(1:3)], rep(NA,3))
      limb_dec       = (mean_q_after-mean_q_before)/time_lag;  #numerical expression for rate of change in discharge (>0, rising, <0 falling) change per time step (5-min)
#      mydata = cbind(mydata, as.data.frame(limb_dec))
      mydata = cbind(mydata, data.frame(limb_dec=limb_dec))
    }

    mydata = data.frame(datenum=datenum, mydata) #convert back to dataframe
    mydata$julian_day=as.POSIXlt(mydata$datenum)$yday                                    #compute julian day
    trig_seq=2*pi*mydata$julian_day
    trig_seq=trig_seq/365
    mydata$sinDOY=sin(trig_seq)
    mydata$cosDOY=cos(trig_seq)
    mydata$julian_day=NULL  

    if (ancillary_data_filename!="")                                              #write file containing ancillary data
    {
      save(file=paste(ancillary_data_filename,".RData",sep=""),list="mydata")     #as binary
      write.table(mydata, paste(ancillary_data_filename,".txt",sep=""), sep="\t", row.names = F, quote=FALSE)   #as ASCII (for checking only)
    }

    #extract training data for ssc
    mydata_training=merge(ssc_data, mydata, all.x=T)
    #mydata_training$discharge = NULL
    #mydata_training = na.omit(mydata_training)
    
    save(file=paste(ancillary_data_filename,"_train_ssc.RData",sep=""),list="mydata_training")  #save as binary
    write.table(mydata_training, paste(ancillary_data_filename,"_train_ssc.txt",sep=""), sep="\t", row.names = FALSE, quote=FALSE) #save as ASCII (for checking only)

    #extract training data for discharge
    #see 2_predictor_generation/5min/record_selection_discharge_model/readme.txt
  }
  