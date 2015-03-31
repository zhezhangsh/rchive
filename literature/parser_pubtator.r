# Parse source files from PubTator
path<-paste(RCHIVE_HOME, 'data/literature/public/pubtator', sep='/'); # home directory of data

if (!file.exists(path)) dir.create(path, recursive=TRUE);
if(!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'));
if(!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'));

# download source files from PubTator FTP server
ftp.files<-c(
    disease="ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator//disease2pubtator.gz",
    chemical="ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator//chemical2pubtator.gz",
    gene="ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator//gene2pubtator.gz",
    species="ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator//species2pubtator.gz",
    mutation="ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator//mutation2pubtator.gz"
)
fn<-paste(path, '/src/', names(ftp.files), '.gz', sep='');
names(fn)<-names(ftp.files);

sapply(names(ftp.files), function(nm) if (!file.exists(fn[nm])) download.file(ftp.files[nm], fn[nm]) else NA)->x;

df<-lapply(names(fn), function(nm) {
  print(nm);
  d<-scan(fn[nm], sep='\n', what='', flush=TRUE);
  hd<-strsplit(d[[1]], '\t')[[1]];
  d<-strsplit(d[-1], '\t');
  df<-lapply(d, function(d) d[1:length(hd)]);
  df<-do.call('rbind', df);
  df<-data.frame(df, stringsAsFactors=FALSE);
  df<-df[df[,1]!='',];
  colnames(df)<-hd;
  save(df, file=paste(path, '/src/', nm, '.rdata', sep=''))
});

    
fn.in<-paste(path, '/src/', names(fn), '.rdata', sep='');
names(fn.in)<-names(fn);

fn.out<-paste(path, '/r/original_pubmed2', names(fn), '.rds', sep='')

cnm<-list();

for (i in 1:length(fn.in)) {
    print(names(fn.in)[i]);
    load(fn.in[i]);
    saveRDS(df, file=fn.out[i]);
    cnm[[i]]<-colnames(df);
}
n<-sapply(fn.out, function(fn) nrow(readRDS(fn))); # number of entries
names(n)<-names(fn.in); 

# update
log<-readRDS(paste(path, 'log.rds', sep='/'));
log<-c(log, list(n));
names(log)[length(log)]<-as.character(Sys.Date());
save(log, file=paste(path, 'log.rds', sep='/'));

