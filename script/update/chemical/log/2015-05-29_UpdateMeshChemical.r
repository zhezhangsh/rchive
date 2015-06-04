library(devtools);
install_github("zhezhangsh/rchive");
library(rchive);
options(stringsAsFactors=FALSE);

##########################
mesh<-ParseMeshChemical();
##########################

##############################################################################################################
UpdateLog(mesh, paste(Sys.getenv("RCHIVE_HOME"), 'data/chemical/public/mesh', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/chemical/UpdateMeshChemical.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/chemical/log/', tm, '_UpdateMeshChemical.r' , sep='');
file.copy(fn0, fn1)
