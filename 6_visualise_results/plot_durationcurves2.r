# merge MC-results and plot duration curves of water and sediment transport

source("../settings.R")
mc_dir = "../5_model_application/"    #directories containing the results of MC-simulation
#load flood number scheme
source(paste(mc_dir,"fload_flood_numbering.R",sep=""))	#load flood numbering scheme
flood_numbering=fload_flood_numbering(gauge_name,individual=T,base_dir=mc_dir)
flood_numbering=flood_numbering[flood_numbering$no>0,]   #discard interflood periods

if(file.exists(file=paste(mc_dir,"merged_series.RData",sep="")))  #load aggregated results, if already present
{
  load(file=paste(mc_dir,"merged_series.RData",sep="")) 
} else
{
  timeseries=data.frame()         #prepare dataframe
  for (i in 6:nrow(flood_numbering))
  {
    source_dir=paste(mc_dir,"MC_",formatC(flood_numbering$no[i],flag="0",digits=1),sep="")  #set source directory
    if (!file.exists(source_dir))
    {
      warning(paste(source_dir,"not found, skipped."))
      next
    }
    print(paste("reading",source_dir)); flush.console()
    tt=read.table(file=paste(source_dir,"/flood_periods_quantForest_",gauge_name,"_ssc_q.txt",sep=""),header=T,sep="\t")  #read reconstructed hydro- and sedigraphs
    tt$datenum=as.POSIXct(strptime(tt$date,"%Y-%m-%d %H:%M:%S"),tz="GMT")
    tt$date =NULL
    timeseries=rbind(timeseries,tt)
  }
  # write.table(file=paste(mc_dir,"merged_series.txt",sep=""),timeseries, sep="\t",row.names=FALSE, quote=FALSE) #save merged data as txt
  
  save(file=paste(mc_dir,"merged_series.RData",sep=""),list="timeseries") #save merged data in binary
}

plot(1,1,xlim=c(0,1),ylim=c(0,1),type="n", axes=F, xlab="", ylab="")
axis(side=2)
axis(side=3)
mtext(side=2,"fraction sediment flux [-]",padj = -3)
mtext(side=3,"fraction time [-]",padj = -3)
abline(0,1)
par(new=TRUE)
plot(1,1,xlim=c(1,0),ylim=c(1,0),type="n", axes=F, xlab="", ylab="")
axis(side=4)
axis(side=1)
mtext(side=4,"fraction water flux [-]",padj = 3)
mtext(side=1,"fraction time [-]",padj = 3)

q90min = ssc90min =  Inf
q90max = ssc90max = -Inf

aucvar=NULL
timespan=as.POSIXlt(range(timeseries$datenum))  #get start and end of time series
for (year in 1900+(timespan[1]$year:timespan[2]$year))
{
  yearstart=ISOdatetime(year,  1,1,0,0,0,tz="GMT")
  yearend  =ISOdatetime(year+1,1,1,0,0,0,tz="GMT")
  valid_records= (timeseries$datenum >= yearstart) & (timeseries$datenum < yearend)
  ssy = sum(timeseries$mean_ssc[valid_records])    #compute sediment yield  [relative units]
  xplot =  (1:sum(valid_records))/sum(valid_records)
  yplot =  cumsum(sort(timeseries$mean_ssc[valid_records])) / ssy
  points(1-yplot, 1-xplot, col="red", pch=".") 
#  if (ssc90min > xplot[min(which(yplot >= 0.9))])
#  {
#    yearssc90min = year
#    ssc90min = xplot[min(which(yplot >= 0.9))]
#    points(1-xplot[min(which(yplot >= 0.9))],1-yplot[min(which(yplot >= 0.9))])
#  }
#
  wy = sum(timeseries$mean_q[valid_records])     #compute water yield     [relative units]
  yplot =  cumsum(sort(timeseries$mean_q[valid_records]  )) / wy
  points( yplot,    xplot, col="blue", pch=".") 
  
  ssy = sum(timeseries$mean_q[valid_records] * timeseries$mean_ssc[valid_records]) * 5*60/1000   #SSY in kg
  qy  = sum(timeseries$mean_q[valid_records]                                     ) * 5*60/1000   #total discharge in m³
  auc=sum(yplot)/(length(yplot)-1)       #compute area under curve
  vari=var(timeseries$mean_ssc[valid_records])                  #compute variance
  cv  =vari/mean(timeseries$mean_ssc[valid_records])           #compute coefficient of variation
  
  aucvar=rbind(aucvar,c(year=year, auc=auc,var=vari,cv=cv))
}
ssy = sum(timeseries$mean_ssc)    #compute sediment yield  [relative units]
wy = sum(timeseries$mean_q)     #compute water yield     [relative units]
points(1-cumsum(sort(timeseries$mean_ssc)) / ssy, 1-(1:nrow(timeseries))/nrow(timeseries), col="black", pch=".", lwd=2) 
points(  cumsum(sort(timeseries$mean_q)) /  wy,     (1:nrow(timeseries))/nrow(timeseries), col="black", pch=".", lwd=2) 

#identify()
#legend("center",col=c("blue","red","black"),lty=1, lwd=c(1,1,2), legend=c("water", "sediment", "all years"), bg="white")

windows()
#plot(rank(aucvar[,1]),rank(aucvar[,2]),xlab="rank auc",ylab="rank var")
plot(aucvar[,"auc"],rank(aucvar[,"var"]),xlab="rank auc",ylab="rank var")

windows()
plot(aucvar[,"auc"],sqrt(aucvar[,"var"]),xlab="auc",ylab="var")

windows()
plot(aucvar[,"auc"],sqrt(aucvar[,"cv"]),xlab="auc",ylab="cv")

