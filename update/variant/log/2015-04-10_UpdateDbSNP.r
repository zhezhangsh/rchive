library(devtools);
source_url("https://raw.githubusercontent.com/zhezhangsh/rchive/master/load.r");

# Save run
tm<-as.character(Sys.Date());
fn0<-paste(RCHIVE_HOME, 'update/variant/UpdateDbSNP.r', sep='/');
fn1<-paste(RCHIVE_HOME, '/update/variant/log/', tm, '_UpdateDbSNP.r' , sep='');
file.copy(fn0, fn1)
