dev.path<-paste(Sys.getenv('RCHIVE_HOME'), 'source/dev', sep='/');

dev.fn<-dir(dev.path, recursive=TRUE);
dev.fn<-paste(dev.path, dev.fn, sep='/');
dev.fn<-dev.fn[grep('.r', dev.fn, ignore.case=TRUE)];

loaded<-sapply(dev.fn, function(f) {
  cat('Loading', f, '\n');
  source(f);
})