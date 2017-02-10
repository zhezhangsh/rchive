library(devtools);
install_github("zhezhangsh/rchive");
library(rchive);

omim<-ParseOmim();

##############################################################################################################
UpdateLog(omim, paste(Sys.getenv("RCHIVE_HOME"), 'data/disease/public/omim', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/disease/UpdateOmim.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/disease/log/', tm, '_UpdateOmim.r' , sep='');
file.copy(fn0, fn1, overwrite = TRUE);
