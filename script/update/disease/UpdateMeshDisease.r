library(devtools);
install_github("zhezhangsh/rchive");
library(rchive);
options(stringsAsFactors=FALSE);

mesh<-ParseMeshDisease();

##############################################################################################################
UpdateLog(mesh, paste(Sys.getenv("RCHIVE_HOME"), 'data/disease/public/mesh', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/disease/UpdateMeshDisease.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/disease/log/', tm, '_UpdateMeshDisease.r' , sep='');
file.copy(fn0, fn1)
