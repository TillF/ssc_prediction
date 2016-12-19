# apply selected model in Monte-Carlo-mode: perform [nrealisations] by randomly choosing the prediction interval quantiles from pred_quantile

# 17.9.2013
# proper date formatting in output files (possibly DST problem before)

# 30.8.2013
#suppress irritating warnings of appending header to file

# 28.5.2013
# multiple cosmetics, output of NA-fraction

# 23.2.2012
# enable passing qrf_model_discharge and qrf_model_ssc as arguments via *_arg to resolve difficult scoping situations

# 22.2.2012
# adaptions to allow returning values as a function
# modified checking of flood_period coverage, disabling of output files

# 25.5.2011
# minor simplifications, comments

# 24.11.2010
# minor corrections for plotting, comments

# 18.11.2010
# include model for discharge in simulation
# improved retrieval of quantiles (speedup)
# improved file handling (overwrite option)

# 7.6.2010
# allow direct specification of flood-scheme (skips reading from file) via parameter flood_scheme
# write_files=FALSE disables all file output

# 31.5.2010
# discard columns not needed by the model


# 20.4.2010
# optional output of IQR-values of all realisations of each timestep (call with iqr=TRUE)
# optional output of entire sampled distribution of each timestep (call with write_dist=TRUE)

# 29.3.2010
# additional output of mean and median value of all realisations of each timestep

# 16.3.2007
# corrected the estimation of CI (former: l´CI of mean freight, now 95 % CI for all values)

# Till Francke, 29.3.2010


# call like:
#  apply_model_mc(gauge_name="h1",nrealisations=100,q_conf=1.00,subset_periods=c(2:2),individual=F,predictordata="ancillary_data2.RData",iqr=TRUE, write_dist=FALSE,overwrite=TRUE)


# gauge_name: "B1", "Villacarli" or "Cabecera"
# subset_preds: index to subset of input data
# nrealisations: number of realisations for estimating confidence interval
# q_conf:	confidence level
# subset_periods: index (i.e. lines), which periods read from the flood_numbering shall be treated
# individual: T:	use individual, F: use general flood numbers from flood_numbering.txt
# predictordata: Rworkspace file containing the predictordata
# iqr: T: include interquartile range columns in reconstructed hydro-/sedigraphs   
# write_dist: T: write entire sampled distribution of each timestep to file (slow, large file, don't use for computing budgets (correlated quantiles)!)
# overwrite: T: overwrite any existing files; F: append to existing time series files

apply_model_mc=function(gauge_name="B1",subset_preds=NULL,nrealisations=100,q_conf=0.95, subset_periods=NULL, individual=T, flood_scheme=NULL, predictordata=NULL, iqr=FALSE, 
  write_dist=FALSE,overwrite=FALSE, write_files=TRUE, qrf_model_discharge_arg=NA, qrf_model_ssc_arg=NA, verbose=TRUE, base_dir="../",
  do_interfloods=TRUE, replaceNA=FALSE)
{
	do_plot=F
  complete_view=NULL   # ALEX: set back to NULL from TRUE
  saves_per_period = 2  # how often memory snapshots (for resuming broken runs) are saved to disk
	
	
  if (is.null(predictordata))
  { # load predictordata from ASCII file   #load predictordata of specified gauge_name prepared with gen_data_4_apply.m
	datain = read.table(paste("../predictors_",gauge_name,".txt",sep = ""), header = TRUE, sep = "\t", quote = "\"'", dec = ".", na.strings = "NaN",strip.white=T,stringsAsFactors=F)
  } else
  {
    load(predictordata) #load from workspace
    datain=mydata
#    datain$ssc      = NULL
   	rm(mydata)
  }
  
	if (replaceNAs)
	{
	  print("replaceNAs enabled, replacing NAs by -999")
	  for (i in 1:ncol(datain))
	  {
	    if (names(datain)[i]=="datenum") next
	    datain[is.na(datain[,i]),i]=-999
	  }  
	}
  
  
  if (is.null(qrf_model_discharge_arg) || 
      !is.na(qrf_model_discharge_arg[1])) qrf_model_discharge = qrf_model_discharge_arg     #if function arguments are given, use them instead of global variables
  if (!is.na(qrf_model_ssc_arg[1]))       qrf_model_ssc       = qrf_model_ssc_arg

  #discard unneeded columns 
  unneeded_columns=which(!(names(datain) %in% c("datenum", dimnames(qrf_model_ssc$importance)[[1]], dimnames(qrf_model_discharge$importance)[[1]])))
  if (is.null(qrf_model_discharge)) unneeded_columns=setdiff(unneeded_columns, which(names(datain) == "discharge")) #preserve discharge, if it is not to be modelled
  for (i in sort(unneeded_columns, decreasing=TRUE))
    datain[,i]=NULL 

	#remove nodata cases
	mydata_raw=na.omit(datain)		
	rm(datain)

	if (is.null(subset_preds))
	{										#use the settings below if not specified from outside
		subset_preds=1:length(mydata_raw[,1])		#use all
	}
		
	subset_preds=subset_preds[which(subset_preds<=nrow(mydata_raw))]		#trim subset_preds to maximum length of time series
	
  mydata_predict=mydata_raw[subset_preds,]	#use subset_preds
	rm(mydata_raw)
	rm(subset_preds)
	
  date_vec=mydata_predict$datenum		


  if (verbose) cat("predictor data loaded\n")
	flush.console()

	if (is.data.frame(flood_scheme))
	{
	 flood_numbering=flood_scheme      #flood-numbering is given directly
	 flood_scheme = TRUE
	} else
	{
    source(paste(base_dir,"fload_flood_numbering.R",sep=""))		#load flood numbering scheme
  	flood_numbering=fload_flood_numbering(gauge_name=gauge_name,individual=individual, do_interfloods=do_interfloods)
  	if (do_interfloods)
  	{
      flood_numbering$begin[1]=min(date_vec)	#set begin of first pre-flood period to begin of time-series 
      if (!is.null(subset_periods))
        subset_periods=sort(c(2*subset_periods,2*subset_periods-1))     #add interfloods to selected periods
    }
    if (verbose) 
    {
      cat("flood numbering scheme loaded")
      if (individual)
    		cat(" (individual flood numbers)\n")  else
    		cat(" (general flood numbers)\n") 
     	flush.console()
  	}
	}	
  return_val=array(0,nrow(flood_numbering))
	
	if (is.null(subset_periods))
		subset_periods=1:nrow(flood_numbering) #use all periods, if not specified from otherwise
	
#	subset_periods=subset_periods[which((subset_periods>=1) & (subset_periods<=nrow(flood_numbering)))]	#limit period selection to available periods
#  subset_periods=which(flood_numbering$no %in% subset_periods)       #convert flood-ids to index of flood_numbering


	mod_name="quantForest"

  #reference vectors to columns to match the same order as in the training data for the forest-models (needed by *Forest :-[)
	ssc_model_column_indices=which(names(mydata_predict) %in% dimnames(qrf_model_ssc$importance)      [[1]]) 
	q_model_column_indices  =which(names(mydata_predict) %in% dimnames(qrf_model_discharge$importance)[[1]]) 


  discharge=mydata_predict$discharge    #save discharge in any case

  # *************************************************************************************************************************************
	# important note****** A: q_conf determines the random quantile!	  *******************************************************************
	# *************************************************************************************************************************************
	quantile_range=c(0.5-q_conf/2,0.5+q_conf/2)	#compute quantile range that is used for drawing the random quantile
	
	file_prefix=paste("flood_periods_",mod_name,"_",gauge_name,sep="")
#	if (individual==F) file_prefix=paste(file_prefix,"_general",sep="")  #mark filename according to flood numbering used

  #determine existence of old files/overwrite behaviour
  if (overwrite)
    append_periodfile=append_timeseriesfile  =append_timestepdistsscfile=append_timestepdistqfile  =TRUE else
  {  
    append_periodfile      =file.exists(file = paste(file_prefix,".txt",sep=""))
    append_timeseriesfile  =file.exists(file = paste(file_prefix,"_ssc_q",".txt",sep=""))
    append_timestepdistsscfile=file.exists(file = paste(file_prefix,"_hist_timestep_ssc",".txt",sep=""))
    append_timestepdistqfile  =file.exists(file = paste(file_prefix,"_hist_timestep_q",  ".txt",sep=""))
  }

	dt=as.numeric(mean(difftime(date_vec[2:100],date_vec[1:99],units="secs")))
	if (verbose) print(paste("timestep resolution:",dt)); flush.console()
	for (j in 1:length(subset_periods))	#treat all desired periods
	{
    if (flood_numbering$end    [subset_periods[j]] < min(date_vec) |
        flood_numbering$begin  [subset_periods[j]] > max(date_vec)   )
		{
      cat(paste("Period",flood_numbering$no[subset_periods[j]],"not covered by timeseries, skipped.\n"))
			next									#no valid data for current period
		}
		
    start_record_ix=which.min(abs(as.numeric(date_vec-flood_numbering$begin[subset_periods[j]])))	#find record in timeset corresponding to start of period
		stop_record_ix =which.min(abs(as.numeric(date_vec-flood_numbering$end  [subset_periods[j]])))	#find record in timeset corresponding to end of period

    #compute fraction of missing data
    expected_timesteps = as.numeric(difftime(flood_numbering$end[subset_periods[j]], flood_numbering$begin[subset_periods[j]], units="sec")) / as.numeric(dt)
    nodata_frac =  (expected_timesteps - (stop_record_ix - start_record_ix)) / expected_timesteps
    if (nodata_frac > 0)
      warning(round(nodata_frac*100), "% datagaps encontered in predictor data for period ",flood_numbering$no[subset_periods[j]],".")

		ssc       =array(0,c(stop_record_ix-start_record_ix+1,nrealisations))	#store simulated values for each timestep and each realisation
 		mean_ssc  =rep  (0, (stop_record_ix-start_record_ix+1))	#for storing mean of simulated values for each timestep and all realisations
 		median_ssc=rep  (0, (stop_record_ix-start_record_ix+1))	#for storing median of simulated values for each timestep and all realisations
    if (iqr) 
      iqr_ssc   =cbind(median_ssc,median_ssc)        #for storing interquartile-range of simulated values for each timestep and all realisations
  	dist_ssc_timestep=rep(0,nrealisations)					#for collecting distribution of ssc of timestep 

		q_sim       =array(0,c(stop_record_ix-start_record_ix+1,nrealisations))	#store simulated values for each timestep and each realisation
 		mean_q  =rep  (0, (stop_record_ix-start_record_ix+1))	#for storing mean of simulated values for each timestep and all realisations
 		median_q=rep  (0, (stop_record_ix-start_record_ix+1))	#for storing median of simulated values for each timestep and all realisations
    if (iqr) 
      iqr_q   =cbind(median_q,median_q)        #for storing interquartile-range of simulated values for each timestep and all realisations
  	dist_q_timestep=rep(0,nrealisations)					#for collecting distribution of ssc of timestep 



		#apply model
		for (i in start_record_ix:stop_record_ix)	#do for all records of current period
		{
			if (verbose) cat(paste("timestep",i-start_record_ix+1,"/",stop_record_ix-start_record_ix+1,"in period",j,"/",length(subset_periods),"\n"))		#display loop counter
			flush.console()

		  if (write_dist)       
      rand_quant=rep(seq(from=quantile_range[1], to=quantile_range[2],length.out=nrealisations)   ,2)	else	#step through quantile range systematically
      rand_quant=runif(2*nrealisations,min=quantile_range[1],max=quantile_range[2])		#draw 2*nrealisations random number from quantile range

			dist_ssc_timestep=as.vector(predict(qrf_model_ssc,      newdata=mydata_predict[i,ssc_model_column_indices],quantiles =rand_quant[ 1:nrealisations]))	#do prediction for this record
      if (!is.null(qrf_model_discharge))
        dist_q_timestep  =as.vector(predict(qrf_model_discharge,newdata=mydata_predict[i,  q_model_column_indices],quantiles =rand_quant[(1:nrealisations)+nrealisations]))	#do prediction for this record
      else
       dist_q_timestep  =rep (discharge[i], nrealisations)	#fake prediction, just use recorded discharge for this record
    
      ssc       [i-start_record_ix+1,]= dist_ssc_timestep #store timestep distribution of SSC
      mean_ssc  [i-start_record_ix+1] = mean  (dist_ssc_timestep)
      median_ssc[i-start_record_ix+1] = median(dist_ssc_timestep)
      
     	q_sim     [i-start_record_ix+1,]= dist_q_timestep		#store timestep distribution of discharge
      mean_q    [i-start_record_ix+1] = mean  (dist_q_timestep)
      median_q  [i-start_record_ix+1] = median(dist_q_timestep)

      if (iqr) 
      {
        iqr_ssc[i-start_record_ix+1,1:2] = quantile(dist_ssc_timestep,c(0.25,0.75))  #interquartile-range 
        iqr_q  [i-start_record_ix+1,1:2] = quantile(dist_q_timestep,  c(0.25,0.75))  #interquartile-range 
      }  

          
      if (write_dist)       #write entire sampled distribution to file
      {
        ww=getOption("warn")
        options(warn=-1) #disable warnings that will arise when trying to convert strings to float
        write.table(file = paste(file_prefix,"_hist_timestep_ssc",".txt",sep=""), data.frame(format(date_vec[i], "%Y-%m-%d %H:%M", tz="GMT"),t(dist_ssc_timestep)), append = append_timestepdistsscfile, quote = F,row.names=F,sep="\t",col.names=!append_timestepdistsscfile)
        write.table(file = paste(file_prefix,"_hist_timestep_q",".txt",sep=""),   data.frame(format(date_vec[i], "%Y-%m-%d %H:%M", tz="GMT"),t(dist_q_timestep)),   append = append_timestepdistqfile,   quote = F,row.names=F,sep="\t",col.names=!append_timestepdistqfile)
#        write.table(file = paste(file_prefix,"_hist_timestep_ssc",".txt",sep=""), data.frame(date_vec[i],t(dist_ssc_timestep)), append = append_timestepdistsscfile, quote = F,row.names=F,sep="\t",col.names=!append_timestepdistsscfile)
#        write.table(file = paste(file_prefix,"_hist_timestep_q",".txt",sep=""),   data.frame(date_vec[i],t(dist_q_timestep)),   append = append_timestepdistqfile,   quote = F,row.names=F,sep="\t",col.names=!append_timestepdistqfile)
        options(warn=ww) #reactivate warnings
        append_timestepdistsscfile=append_timestepdistqfile=TRUE
      }
      ##visualisation of single day
#  			windows()
#  			hist(dist_ssc_timestep,main=paste("SSC at ",date_vec[i]))
#  			abline(v=mean_ssc  [i-start_record_ix+1], col="orange")
#    	  abline(v=median_ssc[i-start_record_ix+1], col="blue")
#  		  legend(x="topright",c("mean","median"),col=c("orange","blue"),lty=1)

      #points(mydata_predict$julian_day[i],mean_ssc_timestep,pch=".",col="red")	#plot mean (for display only)
		
      if ((saves_per_period > 0) && !(i %in% c(stop_record_ix,start_record_ix)))    #dump memory image during period to allow for resuming run
         if ( (i-start_record_ix) %% ((stop_record_ix-start_record_ix)%/% min(stop_record_ix-start_record_ix ,saves_per_period+1)) == 0)
           save.image(file = "resume.RData")
    }	
		sed_flux=ssc * q_sim/1000		#compute sediment flux [kg/s]
		
		sed_sums=apply(sed_flux,2,sum)*dt		#compute mass exported  [kg]
		
		if (is.null(complete_view))
		{
		headerline=c("period_no","mean_yield[kg]", "median_yield[kg]", "CI_lo[kg]","CI_hi[kg]","CI_lo_ok[kg]","CI_hi_ok[kg]","CI_losd[kg]","CI_hisd[kg]","nodata_frac[-]")		#first loop: write headerline
		} else {
		headerline=F
		}
# write.table(cbind(flood_numbering$no[subset_periods[j]],mean(sed_sums),t(t.test(sed_sums, conf.level = q_conf)$conf.int[1:2])), file = paste(file_prefix,".txt",sep=""), append = T, quote = F,row.names=F,col.names=headerline,sep="\t")
# Alex		a=1-q_conf  # former version, now:
    a=0.05            # use always 0.95 conf. limits of the mean! (this step has nothing to do with the step that determines the random quantile!	
		ci_lo=mean(sed_sums)-qt(1-a/2,length(sed_sums)-1)*sd(sed_sums)		#lower confidence limit
		ci_hi=mean(sed_sums)+qt(1-a/2,length(sed_sums)-1)*sd(sed_sums)		#upper confidence limit
		tt=quantile(sed_sums, probs=c(a/2,0.5,1-a/2))		#quantile values and median
	
    ci_lo2=tt[1]
		ci_hi2=tt[3]
		ci_0_5=tt[2] # median

		ci_losd=mean(sed_sums)-sd(sed_sums)		# -1sd
		ci_hisd=mean(sed_sums)+sd(sed_sums)		# +1sd

# ALEX: ci_lo=mean(sed_sums)-qt(1-a/2,length(sed_sums)-1)*(sd(sed_sums)/length(sed_sums)^0.5)		
# ALEX: ci_hi=mean(sed_sums)+qt(1-a/2,length(sed_sums)-1)*(sd(sed_sums)/length(sed_sums)^0.5)	
# ALEX: confidence limits and quantiles are two different concepts, the former describes
#       the uncertainty of an estimator, whereas the latter tells something about the
#       distribution of the data.      
    

##file output
        
    #write flood-based sediment-yield budgets for all periods
   if(!is.null(return_val))
      return_val[j]=mean_ssc[1]   #intended for validation use, so only a single timestep is contained in mean_ssc
    if (write_files)
    {
      ww=getOption("warn")
      options(warn=-1) #disable warnings that will arise when trying to convert strings to float
      write.table(cbind(flood_numbering$no[subset_periods[j]], mean(sed_sums), ci_0_5, ci_lo,ci_hi,ci_lo2,ci_hi2,ci_losd,ci_hisd, nodata_frac), file = paste(file_prefix,".txt",sep=""), append = append_periodfile, quote = F,row.names=F,col.names=headerline,sep="\t")
      options(warn=ww) #reactivate warnings
    }


    #write reconstructed hydro/sedigraph data
    ww=getOption("warn")
    options(warn=-1) #disable warnings that will arise when trying to convert strings to float
    if (write_files)
    if (iqr)
      write.table(data.frame(format(date_vec[start_record_ix:stop_record_ix], "%Y-%m-%d %H:%M", tz="GMT"), mean_q,median_q,iqr_q, mean_ssc,median_ssc,iqr_ssc), file = paste(file_prefix,"_ssc_q",".txt",sep=""), append = append_timeseriesfile, quote = F,row.names=F,col.names=c("date","mean_q","median_q","iqr_lo_q","iqr_hi_q","mean_ssc","median_ssc","iqr_lo_ssc","iqr_hi_ssc"),sep="\t")
    else
      write.table(data.frame(format(date_vec[start_record_ix:stop_record_ix], "%Y-%m-%d %H:%M", tz="GMT"), mean_q,median_q,       mean_ssc,median_ssc        ), file = paste(file_prefix,"_ssc_q",".txt",sep=""), append = append_timeseriesfile, quote = F,row.names=F,col.names=c("date","mean_q","median_q",                      "mean_ssc","median_ssc"                          ),sep="\t")
    options(warn=ww) #reactivate warnings
 
    #write SY-histogram data for current period
    if (write_files)
    write.table(sed_sums, file = paste(file_prefix,"_hist",flood_numbering$no[subset_periods[j]],".txt",sep=""), append = F, quote = F,row.names=F,col.names=c("yield[MC simulations]"),sep="\t")
    
#A: suppress plots
		
		#plot results in complete view
#A		if (is.null(complete_view))
#A		{	
#A			windows()
#A 			plot(0,0,type="l",ylim=c(0,300),xlim=c(min(mydata_predict$julian_day), 	max(mydata_predict$julian_day)),xlab="julian day",ylab="SSC [g/l]")	#prepare plot window
#A			plot(range(date_vec),c(0,0),type="n",ylim=c(0,3),xlab="date",ylab="SSC [g/l]")	#prepare plot window
#A
#A			complete_view=dev.cur()
#A		} else
#A		{
#A			bringToTop(complete_view)	#put focus on overall view
#A			dev.set(complete_view)
#A		}
#		xplot=mydata_predict$julian_day[start_record_ix:stop_record_ix]
#A		xplot=date_vec[start_record_ix:stop_record_ix]

#A		points(rep(xplot,ncol(ssc)),ssc,pch=".")
    #points(rep(mydata_predict$julian_day[start_record_ix:stop_record_ix],ncol(ssc)),ssc,pch=".")
    #		lw_plot=lowess(xplot,apply(ssc,1,mean),f = 0.05)
    #		lines(lw_plot, col="red")
#A		lines(xplot,apply(ssc,1,median), col="red")
		#points(xplot,apply(ssc,1,mean),pch=".",col="red")
		
		if (do_plot)
		{
#plot results of period in separate view
    #  		windows()
    #  		plot(rep(xplot,ncol(ssc)),ssc,pch=".",xlim=mydata_predict$julian_day[c(start_record_ix,stop_record_ix)],xlab="julian day",ylab="SSC[g/l]",main=paste("Period",flood_numbering$no[subset_periods[j]]))
#A  		plot(rep(xplot,ncol(ssc)),ssc,pch=".",xlim=date_vec[c(start_record_ix,stop_record_ix)],xlab="date",ylab="SSC[g/l]",main=paste("Period",flood_numbering$no[subset_periods[j]]))

#	lines(lw_plot, col="red")
#A      lines(xplot,apply(ssc,1,median), col="red")  	
  		#if ((stop_record_ix-start_record_ix)<300)
  		#{
  	#		lines(mydata_predict$julian_day[start_record_ix:stop_record_ix],apply(ssc,1,mean),pch=".",col="red")
  	#		
  	#	}else  {	#draw lines if only few points available
  	#		points(mydata_predict$julian_day[start_record_ix:stop_record_ix],apply(ssc,1,mean),pch=".",col="red")
  	#	}
#A  		savePlot(filename = paste(file_prefix,"_p",flood_numbering$no[subset_periods[j]],sep=""),type = "wmf",device = dev.cur(),restoreConsole = TRUE)	#save plot
  	
      #plot histogram of period in separate view
  		mean_yield  =mean(sed_sums)
  		median_yield=median(sed_sums)
      
      windows()
  		hist(sed_sums,main=paste("Period",flood_numbering$no[subset_periods[j]]),xlab="cum. sediment freight [kg]",xlim=range(c(sed_sums,mean_yield,mean_yield,ci_lo,ci_hi))) #plot histogram of cumulated sediment export
  		lines(c(ci_lo,ci_lo),c(0,1000),col="blue",lty=2)		#add lines of CI limits
  		lines(c(ci_hi,ci_hi),c(0,1000),col="blue",lty=2)		#add lines of CI limits

      	abline(v=mean_yield  , col="orange")
    	  abline(v=median_yield, col="blue")
  		  legend(x="topright",c("mean","median"),col=c("orange","blue"),lty=1)
    
    	par(new=TRUE)
  		xplot=seq(min(sed_sums),max(sed_sums),length.out=100)
  		if (mean(sed_sums)!=0)	plot(xplot,dnorm(xplot,mean=mean(sed_sums),sd=sd(sed_sums)),col = "red",xlab="",ylab="",axes = FALSE,type="l")
  		savePlot(filename = paste(file_prefix,"_p",flood_numbering$no[subset_periods[j]],"hist",sep=""),type = "wmf",device = dev.cur(),restoreConsole = TRUE)	#save plot
		}
	}

#A savePlot(filename = paste(file_prefix,"_complete",sep=""),type = "wmf",device = complete_view,restoreConsole = TRUE)	#save complete view
	
	if(!is.null(return_val))        #return value if run in function mode
    return(return_val)
}
