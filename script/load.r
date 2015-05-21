# load all functions in the RCHIVE_HOME/source/dev directory
dev.path<-paste(Sys.getenv("RCHIVE_HOME"), 'source/dev', sep='/');
dev.fn<-dir(dev.path, recursive=TRUE);
dev.fn<-paste(dev.path, dev.fn, sep='/');
dev.fn<-dev.fn[grep('.r$', dev.fn, ignore.case=TRUE)];

for (i in 1:length(dev.fn)) source(dev.fn[i]);