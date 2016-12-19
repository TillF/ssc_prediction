#create replicate directories of model runs for call_mc.R

rm(list = ls())

seeds=1

dest_dirs=NULL

source("../settings_daily.R")

for (dir_no in 1:n_slaves)
  for (seed in seeds)
  {
    from = ceiling (length(flood_periods) * (dir_no-1)/n_slaves  ) +1
    to=    ceiling (length(flood_periods) * dir_no    /n_slaves  )   
    flood_period=paste(flood_periods[from],":", flood_periods[to], sep="")
    dest_dir=paste("MC_",formatC(dir_no,flag="0",digits=1),sep="")
    dest_dirs = c(dest_dirs,dest_dir)
    dir.create(dest_dir, showWarnings = TRUE)
    file.copy(from=paste("templ",dir("templ"),sep="/"), to=dest_dir)
    write.table(file=paste(dest_dir,"conf.txt",sep="/"),data.frame(seed=seed, flood_period=flood_period),sep="\t",quote=F,row.names=F)
    unknown_keywords=NULL
    for (tfile in dir("templ/dynamic"))
    {
      content = scan(paste("templ/dynamic",tfile,sep="/"),what="character",sep="\n",strip.white=T,quiet=T)

      for (i in 1:length(content))
      {
         rr=gregexpr("\\$\\w*",content[i])   #look for strings starting with $
         if (rr[[1]]!= -1)
         for (j in 1:length(rr[[1]]))
         {
           sstring=substr(content[i],rr[[1]][j]+1,rr[[1]][j]+attr(rr[[1]],"match.length")[j]-1)
           if (exists(sstring))
             content[i] = gsub(patt=paste("\\$",sstring,sep=""),repl=eval(parse(text=sstring)),x=content[i]) else
             unknown_keywords = c(unknown_keywords,sstring)
         }
      }
      write(content, file = paste(dest_dir,"/",tfile,sep=""),append=F)
    }
    
  }
  
  if (length(unknown_keywords)>0)
  warning(paste("Keywords not replaced in template files:",paste(unknown_keywords,collapse=", ")))

  write("#!/bin/bash", file = "start_jobs.sh",append=F)
  write(paste("cd",paste(dest_dirs,"\nqsub job.pbs","\ncd ..\nsleep 5",sep="")), file = "start_jobs.sh",append=T)
  if(.Platform$OS.type != "windows")  #call chmod to make start-script executable
  system("chmod u+x start_jobs.sh")
