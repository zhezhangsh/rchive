#############################################################
download.new<-FALSE; # whether to re-download source files###
#############################################################

# Parse source files from dbSNP
path<-paste(RCHIVE_HOME, 'data/variant/public/dbsnp', sep='/'); # home directory of data
if (!file.exists(path)) dir.create(path, recursive=TRUE);
if(!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'));
if(!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'));

# source file path
vcfs<-c(
  'GRCh37' = "ftp://ftp.ncbi.nlm.nih.gov/snp//organisms/human_9606_b142_GRCh37p13/VCF/All_20150217.vcf.gz",
  'GRCh38' = "ftp://ftp.ncbi.nlm.nih.gov/snp//organisms/human_9606_b142_GRCh38/VCF/All_20150218.vcf.gz"
);


library(VariantAnnotation);
library(dplyr);

# create SQLite database or connect to existing one
db.fn<-paste(RCHIVE_HOME, '/data/variant/db/variant_position.sqlite', sep='');
if (!file.exists(db.fn)) db<-src_sqlite(db.fn, create=TRUE);

log0<-sapply(names(vcfs), function(nm) {
  # Download from source
  fn<-paste(path, '/src/', nm, '_All.vcf.gz', sep='');
  if (download.new | !file.exists(fn)) {
    cat('Downloading ', nm, '...\n');
    download.file(vcfs[nm], fn);
    download.file(paste(vcfs[nm], '.tbi', sep=''), paste(fn, '.tbi', sep=''));
  }
  
  # load VCF
  cat('Loading ', nm, '...\n'); 
  dbsnp<-scanVcf(fn, param=ScanVcfParam(fixed='ALT', info=NA, geno=NA));
  saveRDS(dbsnp[[1]], file=paste(path, '/r/full_vcf_', nm, '.rds', sep=''));
  
  # SNP positions
  gr<-dbsnp[[1]]$rowData;
  fn.pos<-paste(path, '/r/position_', nm, '.rds', sep='');
  if (file.exists(fn.pos)) id.old<-names(readRDS(fn.pos)) else id.old<-c();
  log<-list(N=length(gr), Source=vcfs[nm], Added=setdiff(names(gr), id.old), Removed=setdiff(id.old, names(gr))); # Get update info.
  saveRDS(gr, file=fn.pos);
  
  # create database table
  db<-src_sqlite(db.fn);  
  cat('Creating table', nm, '...\n');
  tbl<-data.frame(id=names(gr), chr=as.vector(seqnames(gr)), pos=start(gr), width=width(gr), stringsAsFactors=FALSE);
  t<-copy_to(db, tbl, temporary=FALSE, names=paste('dbsnp', nm, sep='_'), indexes=list('id', 'chr', 'pos'));
  
  log;
})
