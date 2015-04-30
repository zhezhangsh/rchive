#############################################################
download.new<-TRUE; # whether to re-download source files###
#############################################################

# Parse source files from NCBI Taxonomy
path<-paste(RCHIVE_HOME, 'data/taxonomy/public/ncbi', sep='/'); # home directory of data

if (!file.exists(path)) dir.create(path, recursive=TRUE);
if(!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'));
if(!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'));

# Download source file
ftp.file<-"ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip";
fn<-paste(path, 'src', 'ncbi_taxdmp.zip', sep='/');
if (download.new) {
  download.file(ftp.file, fn);
  wd<-getwd();
  setwd(paste(path, 'src', sep='/'));
  try(unzip(fn));
  setwd(wd);
}

# read in species names
nms.fn<-paste(path, 'src', 'names.dmp', sep='/');
fll<-scan(nms.fn, sep='\n', what='', flush=TRUE);
fld<-strsplit(fll, '\t');
tbl<-do.call(rbind, fld)[, seq(1, 7, 2)];
tbl[is.na(tbl)]<-'';
tbl<-data.frame(tbl, stringsAsFactors=FALSE);
names(tbl)<-c('ID', 'Name', 'Unique_Name', 'Class');

# split
id2name<-sapply(split(tbl[,2], tbl[,1]), unique);
name2id<-sapply(split(tbl[,1], tbl[,2]), unique);
cls<-split(tbl[, 1:3], tbl[[4]]);
cmm<-cls[['common name']];
id2name.cmm<-sapply(split(cmm[,2], cmm[, 1]), unique);
name2id.cmm<-sapply(split(cmm[,1], cmm[, 2]), unique);

# Save data
saveRDS(tbl, file=paste(path, 'r', 'full_table.rds', sep='/'));
saveRDS(id2name, file=paste(path, 'r', 'id2name.rds', sep='/'));
saveRDS(name2id, file=paste(path, 'r', 'name2id.rds', sep='/'));
saveRDS(id2name.cmm, file=paste(path, 'r', 'id2name_common.rds', sep='/'));
saveRDS(name2id.cmm, file=paste(path, 'r', 'name2id_common.rds', sep='/'));
saveRDS(cls, file=paste(path, 'r', 'table_by_class.rds', sep='/'));

##################################################################################
# existing full taxonomy table 
if (download.new) {
  if (file.exists(paste(path, 'r', 'full_table.rds'))) {
    tbl0<-readRDS(paste(path, 'r', 'full_table.rds'));
    id.old<-as.vector(tbl0[,1]);
    nm.old<-as.vector(tbl0[,2]);
  } else {
    id.old<-c();
    nm.old<-c();
  }
  
  # updates
  up<-list(
    N = list(ID = length(unique(tbl[[1]])), Name=length(unique(tbl[[2]]))),
    Added = list(ID = setdiff(tbl[[1]], id.old), Name = setdiff(tbl[[2]], nm.old)),
    Removed = list(ID = setdiff(id.old, tbl[[1]]), Name = setdiff(nm.old, tbl[[2]]))
  )
  
  # update logs
  log<-readRDS(paste(path, 'log.rds', sep='/'));
  log<-c(log, list(up));
  names(log)[length(log)]<-as.character(Sys.Date());
  saveRDS(log, file=paste(path, 'log.rds', sep='/'));
}
