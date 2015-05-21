args<-commandArgs(TRUE); # YAML file name

source(args[1]);

dt<-Sys.Date();
ln<-paste(dt, args[1], sep='\t');
writeLines(dt, paste(Sys.getenv('RCHIVE_HOME', 'source/script/update/update.log')), append=TRUE);

