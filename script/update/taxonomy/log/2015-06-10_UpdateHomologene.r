devtools::install_github("zhezhangsh/rchive");
library(rchive);
options(stringsAsFactors=FALSE);

ids<-ParseTaxonomy();

##############################################################################################################
UpdateLog(ids, paste(Sys.getenv('RCHIVE_HOME'), 'data/taxonomy/public/ncbi', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/taxonomy/UpdateHomologene.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/taxonomy/log/', tm, '_UpdateHomologene.r' , sep='');
file.copy(fn0, fn1)
