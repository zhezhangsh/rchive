# Download and parse the MSigDb database
path<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/msigdb', sep='/');
  
options(stringsAsFactors=FALSE);

if (!file.exists(path)) dir.create(path, recursive=TRUE);
if (!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'));
if (!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'));

library(rchive);

fn<-dir(paste(path, 'src', sep='/'));
fn<-fn[grep('.gmt$', fn, ignore.case=TRUE)];
fn<-paste(path, 'src', fn, sep='/');

nm<-sapply(strsplit(fn, '/'), function(x) x[length(x)]);
nm<-sub('.gmt$', '', nm);
names(fn)<-nm;

meta<-lapply(names(fn), function(nm) {
  gs<-scan(fn[nm], what='', flush=TRUE, sep='\n');
  gs<-strsplit(gs, '\t');
  gs.nm<-sapply(gs, function(gs) gs[1]);
  gs.url<-sapply(gs, function(gs) gs[2]);
  gs<-lapply(gs, function(gs) gs[3:length(gs)]);
  meta<-data.frame(N=sapply(gs, length), URL=gs.url, row.names=gs.nm, stringsAsFactors=FALSE);
  names(gs)<-gs.nm;
  saveRDS(list(meta=meta, gene.sets=gs), paste(path, '/r/', nm, '.rds', sep=''));
  meta;
});

UpdateLog(meta, path);
tm<-as.character(Sys.Date());
fn0<-paste(Sys.getenv('RCHIVE_HOME'), 'source/script/update/gene.set/UpdateMSigDb.r', sep='/');
fn1<-paste(Sys.getenv('RCHIVE_HOME'), '/source/script/update/gene.set/log/', tm, '_UpdateMSigDb.r' , sep='');
file.copy(fn0, fn1)


