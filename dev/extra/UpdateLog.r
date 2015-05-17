# Update log
UpdateLog<-function(updates, path, just.new=FALSE) {
  # updates   A list of updates
  # path      Path to the log file
  # just.new  Only save the new entries
  
  fn<-paste(path, 'log.rds', sep='/');
  
  if (!file.exists(fn)) log<-updates else {
    log<-readRDS(fn);
    if (!just.new) log[[as.character(Sys.Date())]]<-updates else {
      lst<-lapply(names(updates), function(nm) unlist(lapply(log, function(log) log$nm), use.names=FALSE));
      names(lst)<-names(updates);
      up<-lapply(names(updates), function(nm) setdiff(updates[[nm]], lst[[nm]]));
      up$removed<-lapply(names(updates), function(nm) setdiff(lst[[nm]], updates[[nm]]));
      log[[as.character(Sys.Date())]]<-up;
    }
  }
  
  saveRDS(log, fn);
}