library(devtools);
install_github("zhezhangsh/rchive");
library(rchive);



##############################################################################################################
# UpdateLog(ids, paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/unigene', sep='/'), just.new=FALSE);
# 
# tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
# fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene/UpdateUnigene.r', sep='');
# fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene/log/', tm, '_UpdateUnigene.r' , sep='');
# file.copy(fn0, fn1); 
