# generate ancillary predictors from rainfall data
 
source("../settings.R")

precip_data_file = sub(precip_data_file, pattern = "\\..*$", replacement = "")
precip_data_file = sub(precip_data_file, pattern = ".*/([^/]+)$", replacement = "\\1", perl = TRUE)
precip_data_file = paste0("1_rainfall_correction/",precip_data_file,"_corrected.txt")

precip_data = read.table(paste("../",precip_data_file, sep=""), header = TRUE, sep = "\t", dec = ".", na.strings = c("NaN","NA"), stringsAsFactors=F)
precip_data$datenum    = as.POSIXct(strptime(precip_data$date, format=fmt_precip, tz="GMT"))  #convert string to date
precip_data$date       = NULL                                    #discard string colummn

ssc_data = read.table(paste("../",ssc_data_file, sep=""), header = TRUE, sep = "\t", dec = ".", na.strings = c("NaN","NA"), stringsAsFactors=F)
ssc_data$datenum    = as.POSIXct(strptime(ssc_data$date, format=fmt_ssc, tz="GMT"))  #convert string to date
ssc_data$date       = NULL                                      #discard string colummn
#move SSC-samples to nearest sampling point in time corresponding to time series
ssc_data$datenum = 
  min(precip_data$datenum) + round(as.numeric(ssc_data$datenum - as.numeric(min(precip_data$datenum)))/(as.numeric(tres*60)))*(as.numeric(tres*60))


discharge_data = read.table(paste("../",discharge_data_file, sep=""), header = TRUE, sep = "\t", dec = ".", na.strings = c("NaN","NA"), stringsAsFactors=F)
discharge_data$datenum    = as.POSIXct(strptime(discharge_data$date, format=fmt_discharge, tz="GMT"))  #convert string to date
discharge_data$date       = NULL                                      #discard string colummn
#discharge_data=NULL #don't use any discharge data

rain_gauges = head(names(precip_data),-1)

if (!is.null(discharge_data))
  mydata=merge(precip_data,discharge_data,all = TRUE) else#merge precip-data with discharge-data  
  mydata=precip_data
names(mydata)=c("datenum",rain_gauges,"discharge")[1:ncol(mydata)]                                               #rename columns of data frame
rm(discharge_data)
rm(precip_data)
mydata$julian_day=as.POSIXlt(mydata$datenum)$yday+1


contiguous_time = seq(from=min(mydata$datenum), to=max(mydata$datenum), by=tres)  #create contiguous time line (opposed to intermittent bucket data)
mydata = merge(data.frame(datenum=contiguous_time), mydata, all = TRUE)                        
names(mydata) = c("datenum",rain_gauges,"discharge", "julian_day")                                               #rename columns of data frame
rm(contiguous_time)

#mydata=mydata[min(which(!is.na(mydata[,"discharge"]))) +(1:1000),]     ##rr test option


##settings for ancillary predictors to be derived
#name:                      name of primary predictor
#n_aggregation_levels:       number of progressive accumulation levels  
#aggregation_steps_factor:  progression factor of aggregation time between levels

anc_predictors_settings= data.frame(
  name=c(rain_gauges, NULL),
  n_aggregation_levels = c(rep(9, length(rain_gauges)),  NULL),
  aggregation_steps_factor = c(rep(3, length(rain_gauges)), NULL),
  stringsAsFactors = FALSE     
)
comp_time_range=function(aggregation_steps_factor, n_aggregation_levels) #compute number of timesteps the ancillary predictors look into the past
   return(sum(aggregation_steps_factor^(0:n_aggregation_levels)))   

  paste("length of timespan looked into the past by ancillary predictors:") 
  (cbind(pred=anc_predictors_settings$name, range_minutes=tres*mapply(comp_time_range,anc_predictors_settings$aggregation_steps_factor, anc_predictors_settings$n_aggregation_levels)))

  ancillary_data_filename = paste("ancillary_data_",tres,sep="") #name of file to store the generated data to

  source("compute_ancillary_predictors.R") #don't call as a function to avoid memory replication

  
   cat("finished.\7\nPlease proceed with ../3_record_selection_discharge_model (see readme.txt)\n")

   
