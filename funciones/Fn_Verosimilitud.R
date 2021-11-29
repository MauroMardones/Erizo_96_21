
Verosimilitud<-function(admb,dir.0,dir.1,Carpeta,logRo,priorRo){
  
dir<-paste(dir.0,Carpeta,sep="")

dat_admb<-paste(admb,".dat",sep="")
tpl_admb<-paste(admb,".tpl",sep="")
rep_admb<-paste(admb,".rep",sep="")
std_admb<-paste(admb,".std",sep="")


unlink(dir,recursive=T) #borra la carpeta "CBA2016"
dir.create(file.path(dir.0,Carpeta))#crea la carpeta "CBA2016"" nuevamente
setwd(dir.1);file.copy(c(dat_admb,tpl_admb),dir) #copia los archivos de la carpeta MAE0316
setwd(dir)

system(paste("~/admb-12.2/admb",admb,sep=" "))
system(paste("./",admb,sep=""))

rep.0       <- reptoRlist(rep_admb)
data        <- lisread(paste(dir,dat_admb, sep='/'))
names(data) <- str_trim(names(data), side="right")
data.1      <- data

log_Ro    <- logRo
prior_Ro <- priorRo
casos    <- length(prior_Ro)
#==================================================
#  CREA LOS ARCHIVOS .DAT DE CADA CASO
#======================================================
for(i in 1:casos){
  data.1$log_priorRo     <- prior_Ro[i]
  data.1$opt_Ro          <- -1
  data.1$Fase_desvRt     <- 1
  
  writeData(paste(admb,"s",i,".dat",sep=""), data.1, append=F)
  
  setwd(dir.1)
  file.copy(c(paste(admb,".tpl",sep="")),dir)
  setwd(dir)
  file.rename(paste(admb,".tpl",sep=""),paste(admb,"s",i,".tpl",sep="")) 
  
  system(paste("~/admb-12.2/admb ",admb,"s",i,sep=""))
  system(paste("./",admb,"s",i,sep="")) 
  
}

}

