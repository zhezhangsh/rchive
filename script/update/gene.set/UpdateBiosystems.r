library(devtools);
source_url("https://raw.githubusercontent.com/zhezhangsh/rchive/master/load.r");

options(stringsAsFactors=FALSE);

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
  'yeast' = '4932',
  'ecoli' = '511145'
);

# Path to output files
out<-paste(RCHIVE_HOME, 'data/gene.set/public/biosystems', sep='/'); 

# Download data and save results
ids<-ParseBiosystemsGeneral(species, download.new=TRUE); 

# Mapping BSID to gene, protein, etc.
cat("Mapping BSID to gene, etc.\n");
bsid<-readRDS(paste(RCHIVE_HOME, 'data/gene.set/public/biosystems/r/biosystem_all.rds', sep='/'))
Map2Biosystems(bsid, species, path=paste(RCHIVE_HOME, 'data/gene.set/public/biosystems', sep='/'));

# Save run
tm<-as.character(Sys.Date());
fn0<-paste(RCHIVE_HOME, 'source/update/gene.set/UpdateBiosystems.r', sep='/');
fn1<-paste(RCHIVE_HOME, '/source/update/gene.set/log/', tm, '_UpdateBiosystems.r' , sep='');
file.copy(fn0, fn1)

