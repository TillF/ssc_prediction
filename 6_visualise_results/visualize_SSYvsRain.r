# visualize SSY (results of MC-runs) vs. annual rainfall

#load rainfall data
source("../settings.R")
load(file=paste("../2_predictor_generation/ancillary_data_",tres,sep=""))           
  mydata=mydata[,c("datenum","rain")]  #use only rainfall records

mc_dir = "../5_model_application/"    #directories containing the results of MC-sumulation
#load flood number scheme
source(paste(mc_dir,"fload_flood_numbering.r",sep=""))	#load flood numbering scheme
flood_numbering=fload_flood_numbering(gauge_name,individual=T,base_dir=mc_dir)
flood_numbering=flood_numbering[flood_numbering$no>0,]   #discard interflood periods


ssyields=data.frame()         #prepare dataframe
for (i in 6:nrow(flood_numbering))
{
  source_dir=paste(mc_dir,"MC_",formatC(flood_numbering$no[i],flag="0",digits=1),sep="")  #set source directory
  if (!file.exists(source_dir))
  {
    warning(paste(source_dir,"not found, skipped."))
    next
  }
  tt=read.table(file=paste(source_dir,"/flood_periods_quantForest_",gauge_name,"_hist",flood_numbering$no[i],".txt",sep=""),header=T,sep="\t")  #read realisations of MC-simulation
  annual_rain=sum(mydata$rain[ mydata$datenum>=flood_numbering$begin[i] &  mydata$datenum < flood_numbering$end[i]])     #compute annual rainfall

  ssyields=rbind(ssyields,data.frame(year=as.POSIXlt(flood_numbering$begin[i])$year+1900,annual_rain= annual_rain, yield=tt[,1]))

}

ssyields$yield=log(ssyields$yield)

lin_regression <- lm(yield ~ annual_rain, data=ssyields)        # do linear regression

#plot(lin_regression)                 # sequence of diagnostic plots. See ?plot.lm

# Create prediction values and confidence limits
# using a new dataframe of annual_rain values
newData = data.frame(annual_rain = seq(min(ssyields$annual_rain), max(ssyields$annual_rain), by = (max(ssyields$annual_rain) - min(ssyields$annual_rain)) / 49))
pred_lim = predict(lin_regression, newdata = newData, interval = "prediction", level=0.95)
conf_lim = predict(lin_regression, newdata = newData, interval = "confidence", level=0.95)

conf_lim[conf_lim < 0] = NA
pred_lim[pred_lim < 0] = NA

plot(ssyields$annual_rain, exp(ssyields$yield), xlab = "annual_rain (mm)", ylab = "SSY (kg)", ylim = c(0,max(exp(ssyields$yield), pred_lim, na.rm = TRUE)),pch=".") # Create the plot

matlines(newData$annual_rain, exp(pred_lim), lty = c(1, 4, 4), lwd = 1) # Draw the fitted regression line and the prediction and confidence intervals
matlines(newData$annual_rain, exp(conf_lim), lty = c(1, 3, 3), lwd = 2)

legend("topleft", max(pred_lim, na.rm = TRUE), legend = c("Fitted Line", "Confidence Bands", "Prediction Bands"),
lty = c(1, 3, 4), lwd = c(2,1,2) , horiz = FALSE, cex = 0.9) # Draw the legend

# use mean annual SSY instead of MC-distribution
# ssyieldsMC = ssyields
# ssyields = ssyieldsMC
# ssyields = aggregate(ssyieldsMC,by=list(ssyields$year),FUN=mean)
#windows()
#hist(log(ssyieldsMC$yield))