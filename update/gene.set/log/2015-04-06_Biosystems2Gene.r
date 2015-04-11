library(devtools);
source_url("https://raw.githubusercontent.com/zhezhangsh/rchive/master/load.r");

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
ids<-ParseBiosystems(species); 

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(RCHIVE_HOME, 'update/gene.set/Biosystems2Gene.r', sep='/');
fn1<-paste(RCHIVE_HOME, '/update/gene.set/log/', tm, '_Biosystems2Gene.r' , sep='');
file.copy(fn0, fn1)

