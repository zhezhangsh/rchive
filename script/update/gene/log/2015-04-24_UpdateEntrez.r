library(devtools);
source_url("https://raw.githubusercontent.com/zhezhangsh/rchive/master/load.r");

ftp.file<-"ftp://ftp.ncbi.nih.gov/gene//DATA/GENE_INFO/All_Data.gene_info.gz";

# Selected species
species<-c(
  'human' = '9606',
  'mouse' = '10090',
  'rat' = '10116',
  'chimp' = '9598',
  'pig' = '9823',
  'chicken' = '9031',
  'dog' = '9615',
  'cow' = '9913',
  'worm' = '6239',
  'fly' = '7227',
  'zebrafish' = '7955',
  'ecoli' = '511145'
);

path=paste(RCHIVE_HOME, 'data/gene/public/entrez', sep='/');
ParseEntrez(ftp.file, species, TRUE, path);

##############################################################################################################
tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(RCHIVE_HOME, 'source/update/gene/UpdateEntrez.r', sep='/');
fn1<-paste(RCHIVE_HOME, '/source/update/gene/log/', tm, '_UpdateEntrez.r' , sep='');
file.copy(fn0, fn1)
