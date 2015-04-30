# Parse source files from PubTator
ParsePubtator<-function(ftp.files, path=paste(RCHIVE_HOME, 'data/literature/public/pubtator', sep='/'), download.new=TRUE) {
  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  if(!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'));
  if(!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'));
  
  fn<-paste(path, '/src/', names(ftp.files), '.gz', sep='');
  names(fn)<-names(ftp.files);
  fn.out<-paste(path, '/r/original_pubmed2', names(fn), '.rds', sep='')
  names(fn.out)<-names(fn);
  
  sapply(names(ftp.files), function(nm) if (download.new | !file.exists(fn[nm])) download.file(ftp.files[nm], fn[nm]) else NA)->x;
  
  n<-lapply(names(fn), function(nm) {
    cat(nm, '\n');
    d<-scan(fn[nm], sep='\n', what='', flush=TRUE);
    hd<-strsplit(d[[1]], '\t')[[1]];
    d<-strsplit(d[-1], '\t');
    df<-lapply(d, function(d) d[1:length(hd)]);
    df<-do.call('rbind', df);
    df<-data.frame(df, stringsAsFactors=FALSE);
    df<-df[df[,1]!='',];
    colnames(df)<-hd;
    saveRDS(df, fn.out[nm]);
    id1<-strsplit(df[, 2], ',');
    id2<-rep(df[[1]], sapply(id1, length));
    id1<-unlist(id1, use.names=FALSE);
    mp<-sapply(split(id1, id2), unique);
    saveRDS(mp, file=paste(path, '/r/', nm, '2pubmed.rds', sep=''));
    mp<-sapply(split(id2, id1), unique);
    saveRDS(split(df[, 2], df[, 1]), file=paste(path, '/r/pubmed2', nm, '.rds', sep=''));
    #save(df, file=paste(path, '/src/', nm, '.rdata', sep=''))
    nrow(df);
  });
  names(n)<-names(fn);
  
  fn.in<-paste(path, '/src/', names(fn), '.rdata', sep='');
  names(fn.in)<-names(fn);
   
  # update
  if (download.new) {
    log<-readRDS(paste(path, 'log.rds', sep='/'));
    #log<-c(log, list(n));
    #names(log)[length(log)]<-as.character(Sys.Date());
    log[[as.character(Sys.Date())]]<-n;
    saveRDS(log, file=paste(path, 'log.rds', sep='/'));
  }
};

