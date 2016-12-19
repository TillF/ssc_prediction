#loads flood numbering scheme from flood_numbering.txt
#returns individual/general flood number, begin, end
#call like fload_flood_numbering("B1",F)

fload_flood_numbering=function(gauge_name,individual=T,base_dir="../", do_interfloods=TRUE)
{
	floodfile=paste(base_dir,"flood_numbering.txt",sep="")
	flood_data = read.table(floodfile, header = TRUE, sep = "\t", quote = "\"'", dec = ".", na.strings = "NaN",strip.white=T,stringsAsFactors=F,skip=1,fill=T)
	
	valid_data=flood_data[which(flood_data$Location==gauge_name | flood_data$Location=="all"),]	#select desired location
	rm(flood_data)
	if (individual){
		numbering_col=which(names(valid_data)=="individual_flood_number")	#use column with individual flood numbering scheme
	} else {	
		numbering_col=which(names(valid_data)=="general_flood_number")		#use column with general flood numbering scheme
	}
	
	nrows=nrow(valid_data)
	
	i=2
	while (i < nrows)			 #union entries having the same flood number (only relevant for merging periods with the same ID)
	{
		if (valid_data[i,numbering_col]==valid_data[i-1,numbering_col])
		{
			valid_data$end[i-1]=valid_data$end[i]		#previous entry has the same flood number - merge into single entry
			valid_data=valid_data[-i,]				#removecurrent entry
			nrows=nrows-1
		} else
		i=i+1
	}
	
	nrows=nrow(valid_data)

	tt=data.frame(no=as.integer(valid_data[,numbering_col]),begin=strptime(valid_data$begin,"%d.%m.%Y %H:%M",tz="GMT"),end=strptime(valid_data$end,"%d.%m.%Y %H:%M",tz="GMT"))		#convert date-strings to date-class
  
  if (do_interfloods)
	{
    tt2 = rbind(tt, data.frame(no=-tt$no, begin=c(tt$end[1]+NA, tt$end[-nrows]), end= tt$begin)) #append pre-flood periods (interfloods)
#	  order_ix=rbind((1:nrows)+nrows,1:nrows)[1:(2*nrows)]		#order chronologically
	} else
  tt2 = tt  
	
  order_ix=sort.int(tt2$end, index.return=TRUE)$ix
  
	flood_info=data.frame(no=tt2$no[order_ix],begin=tt2$begin[order_ix],end=tt2$end[order_ix])	
	return(flood_info)
}
		