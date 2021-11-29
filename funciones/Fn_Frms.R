
fn_Frms<-function(admb,dir.0,dir.1,Carpeta,Dat){

  dir<-paste(dir.0,Carpeta,sep="")
  
  tpl_admb<-paste(admb,".tpl",sep="")
  
  
  unlink(dir,recursive=T) #borra la carpeta "CBA2016"
  dir.create(file.path(dir.0,Carpeta))#crea la carpeta "CBA2016"" nuevamente
  setwd(dir.1);file.copy(c(tpl_admb),dir) #copia los archivos de la carpeta MAE0316
  setwd(dir)
  
  writeData(paste(admb,".dat",sep=""), Dat, append=F)
  system(paste("~/admb-12.2/admb ",admb,sep=""))
  system(paste("./",admb,sep="")) 
  
}