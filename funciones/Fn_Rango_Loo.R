Rango_Loo<-function(admb,dir.0,dir.1,Carpeta,rango_Loo){
  
  dir<-paste(dir.0,Carpeta,sep="")
  
  dat_admb<-paste(admb,".dat",sep="")
  tpl_admb<-paste(admb,".tpl",sep="")
  rep_admb<-paste(admb,".rep",sep="")
  std_admb<-paste(admb,".std",sep="")
  
  
  unlink(dir,recursive=T) #borra la carpeta "CBA2016"
  dir.create(file.path(dir.0,Carpeta))#crea la carpeta "CBA2016"" nuevamente
  setwd(dir.1);file.copy(c(dat_admb,tpl_admb),dir) #copia los archivos de la carpeta MAE0316
  setwd(dir)
  
 # system(paste("~/admb-12.2/admb",admb,sep=" "))
#  system(paste("./",admb,sep=""))
  
 # rep.0       <- reptoRlist(rep_admb)
  data        <- lisread(paste(dir,dat_admb, sep='/'))
  names(data) <- str_trim(names(data), side="right")
  data.1      <- data
  
  
  rangoLoo <- rango_Loo
  casos    <- length(rangoLoo)
  #==================================================
  #  CREA LOS ARCHIVOS .DAT DE CADA CASO
  #======================================================
  for(i in 1:casos){
    data.1$Loo_k_Lo_cv_M[1]     <- rangoLoo[i]
    
    writeData(paste(admb,"s",i,".dat",sep=""), data.1, append=F)
    
    setwd(dir.1)
    file.copy(c(paste(admb,".tpl",sep="")),dir)
    setwd(dir)
    file.rename(paste(admb,".tpl",sep=""),paste(admb,"s",i,".tpl",sep="")) 
    
    system(paste("~/admb-12.2/admb ",admb,"s",i,sep=""))
    system(paste("./",admb,"s",i,sep="")) 
    
    
    file.remove(paste(admb,"s",i,".htp", sep=""),
                paste(admb,"s",i,".cpp", sep=""),
                paste(admb,"s",i,".obj", sep=""),
                paste(admb,"s",i,".p01", sep=""),
                paste(admb,"s",i,".b01", sep=""),
                paste(admb,"s",i,".r01", sep=""),
                paste(admb,"s",i,".p02", sep=""),
                paste(admb,"s",i,".b02", sep=""),
                paste(admb,"s",i,".r02", sep=""),
                paste(admb,"s",i,".p03", sep=""),
                paste(admb,"s",i,".b03", sep=""),
                paste(admb,"s",i,".r03", sep=""),
                paste(admb,"s",i,".p04", sep=""),
                paste(admb,"s",i,".b04", sep=""),
                paste(admb,"s",i,".r04", sep=""),
                paste(admb,"s",i,".p05", sep=""),
                paste(admb,"s",i,".b05", sep=""),
                paste(admb,"s",i,".r05", sep=""),
                paste(admb,"s",i,".p06", sep=""),
                paste(admb,"s",i,".b06", sep=""),
                paste(admb,"s",i,".r06", sep=""),
                paste(admb,"s",i,".p07", sep=""),
                paste(admb,"s",i,".b07", sep=""),
                paste(admb,"s",i,".r07", sep=""),
                paste(admb,"s",i,".p08", sep=""),
                paste(admb,"s",i,".b08", sep=""),
                paste(admb,"s",i,".r08", sep=""),
                paste(admb,"s",i,".p09", sep=""),
                paste(admb,"s",i,".b09", sep=""),
                paste(admb,"s",i,".r09", sep=""),
                paste(admb,"s",i,".p10", sep=""),
                paste(admb,"s",i,".b10", sep=""),
                paste(admb,"s",i,".r10", sep=""),
                paste(admb,"s",i,".par", sep=""),
                paste(admb,"s",i,".bar", sep=""),
                paste(admb,"s",i,".eva", sep=""),
                paste(admb,"s",i,".cor", sep=""),
                paste(admb,"s",i,".log", sep=""),
                paste(admb,"s",i,".tpl", sep=""),
                paste(admb,"s",i,".exe", sep=""))
    
    
  }
  
}

