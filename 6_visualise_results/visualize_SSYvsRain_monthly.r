# visualize SSY (results of MC-runs) vs. monthly rainfall

#load rainfall data
source("../settings.R")
paste("../2_predictor_generation/ancillary_data_",tres,sep="")
  mydata=mydata[,c("datenum","rain")]  #use only rainfall records

mc_dir = "../5_model_application_monthly/"    #directories containing the results of MC-sumulation
#load flood number scheme
source(paste(mc_dir,"fload_flood_numbering.r",sep=""))	#load flood numbering scheme
flood_numbering=fload_flood_numbering(gauge_name,individual=T,base_dir=mc_dir)
flood_numbering=flood_numbering[flood_numbering$no>0,]   #discard interflood periods


ssyields=data.frame()         #prepare dataframe
mc_threads=dir(path=mc_dir,pattern="MC_")
for (mc_thread in mc_threads)
{
  source_dir=paste(mc_dir,mc_thread,sep="")  #set source directory
  if (!file.exists(source_dir))
  {
    warning(paste(source_dir,"not found, skipped."))
    next
  }
  period_files=dir(path=source_dir,pattern="flood_periods_quantForest_",gauge_name,"_hist")
  
  for (period_file in period_files)
  {
    tt=read.table(file=paste(source_dir,"/",period_file,sep=""),header=T,sep="\t")  #read realisations of MC-simulation
    rr=regexpr(pattern="[0-9]+\\.",text=period_file)
    flood_period=as.numeric(substr(period_file,start=rr,rr+attr(rr,"match.length")-2)) #extract flood number form file name
    i = flood_numbering$no[flood_numbering$no==flood_period]                          #get row in flood_numbering
    period_rain=sum(mydata$rain[ mydata$datenum>=flood_numbering$begin[i] &  mydata$datenum < flood_numbering$end[i]])     #compute rainfall of respective period
  
    ssyields=rbind(ssyields,data.frame(flood_period=flood_period,rain= period_rain, yield=tt[,1]))  
  }
}

save(list="ssyields", file="SSYmonthly.RData") #save for later use
#load(file="SSYmonthly.RData")                 #restore from file

ssyields$yield=log(ssyields$yield)

lin_regression <- lm(yield ~ rain, data=ssyields)        # do linear regression

#plot(lin_regression)                 # sequence of diagnostic plots. See ?plot.lm

# Create prediction values and confidence limits
# using a new dataframe of period_rain values
newData = data.frame(rain = seq(min(ssyields$rain), max(ssyields$rain), by = (max(ssyields$rain) - min(ssyields$rain)) / 49))
pred_lim = predict(lin_regression, newdata = newData, interval = "prediction", level=0.95)
conf_lim = predict(lin_regression, newdata = newData, interval = "confidence", level=0.95)

conf_lim[conf_lim < 0] = NA
pred_lim[pred_lim < 0] = NA

plot(ssyields$rain, exp(ssyields$yield), xlab = "rain (mm)", ylab = "SSY (kg)", ylim = c(0,max(exp(ssyields$yield), pred_lim, na.rm = TRUE)),pch=".") # Create the plot

matlines(newData$rain, exp(pred_lim), lty = c(1, 4, 4), lwd = 1) # Draw the fitted regression line and the prediction and confidence intervals
matlines(newData$rain, exp(conf_lim), lty = c(1, 3, 3), lwd = 2)

legend("topleft", max(pred_lim, na.rm = TRUE), legend = c("Fitted Line", "Confidence Bands", "Prediction Bands"),
lty = c(1, 3, 4), lwd = c(2,1,2) , horiz = FALSE, cex = 0.9) # Draw the legend



# use mean annual SSY instead of MC-distribution
ssyields_aggr = aggregate(ssyields,by=list(ssyields$flood_period),FUN=mean)
windows()
hist(ssyields_aggr$yield, main="log yield")

load(file="SSYmonthly.RData")                 #restore from file
#2D kernel-density plots
rx=range(ssyields$rain)
ry=range(ssyields$yield)
poly=cbind(rx[c(1,1,2,2)] ,ry[c(1,2,2,1)])#needed by kernel2d to limit the area of the density function

#using splancs
library(splancs)
windows()
k2d=kernel2d(as.points(x=ssyields$rain,y=ssyields$yield), poly, h0=20, nx=50, ny=50)  #adjust h0, nx, ny
plot(poly, type="n")
image(k2d, add=TRUE, col=terrain.colors(20))

#using MASS
library(MASS)
windows()
k2d=kde2d(x=ssyields$rain,y=ssyields$yield, n=c(25,50),h=1000) #adjust n, h
plot(poly, type="n")
image(k2d, add=TRUE, col=terrain.colors(20))

#using KernSmooth
library(KernSmooth)
est <- bkde2D( ssyields[,c("rain","yield")], gridsize=c(25,50), bandwidth=c(100,100)) #adjust gridsize and bandwith
contour(est$x1, est$x2, est$fhat)
persp(est$fhat)

#selfmade - estimate 2D-Kernel-density for each x-bin and transfer to regular grid (nxbins x nybins)
x=ssyields$rain
y=ssyields$yield
nxbins=20 #number of bins in x-direction
nybins=50 #number of bins in x-direction
dens_matrix=array(0,c(nxbins,nybins))
yvec=seq(from=min(y), to=max(y), length.out=nybins)
for (i in 1:nxbins)
{
  xmin=min(x) + (i-1) * diff(range(x))/nxbins
  xmax=min(x) + (i  ) * diff(range(x))/nxbins
  valid_records= (x >= xmin) & (x <= xmax)
  if (sum(valid_records)==0) next
  dens=density(x=y[valid_records])
  dens_matrix[i,]=approx(x=dens$x, y=dens$y, xout=yvec,yleft=0,yright=0)$y
}


mat=dens_matrix#no transformation

#mat=mat / apply (mat,1,sum) #normalize columns to sum up to 1
mat=log(mat)#transform to improve visibility
mat[mat==-Inf]=NA #mask out Inf values

contour(mat)
persp(mat)
image(mat)

