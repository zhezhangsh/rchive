library(devtools);
install_github("zhezhangsh/rchive");
library(rchive);

options(stringsAsFactors=FALSE);

RCHIVE_HOME<-Sys.getenv("RCHIVE_HOME");

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
ids<-ParseBiosystemsGeneral(species, download.new=FALSE); 

# Mapping BSID to gene, protein, etc.
cat("Mapping BSID to gene, etc.\n");
bsid<-readRDS(paste(RCHIVE_HOME, 'data/gene.set/public/biosystems/r/biosystem_all.rds', sep='/'))
Map2Biosystems(bsid, species, path=paste(RCHIVE_HOME, 'data/gene.set/public/biosystems', sep='/'));

# Mapping GO to species specific genes
go.lst<-readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/biosystems/r/biosystem2gene_conserved_GO.rds', sep='/'));
go2gn.by.species<-MapSet2Species(go.lst, names(species));
fn.prefix<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/biosystems/r/', sep='/');
go.anno<-readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/biosystems/r/biosystem_GO.rds', sep='/'))[, -5];
go2gn.by.species<-lapply(go2gn.by.species, function(lst) lst[names(lst) %in% rownames(go.anno)]);
f<-sapply(names(go2gn.by.species), function(nm) 
  saveRDS(go2gn.by.species[[nm]], file=paste(fn.prefix, 'biosystem2gene_', nm, '_GO.rds', sep='')));
f<-sapply(names(go2gn.by.species), function(nm) 
  saveRDS(go.anno[names(go2gn.by.species[[nm]]), , drop=FALSE], file=paste(fn.prefix, 'biosystem_', nm, '_GO.rds', sep='')));

# Save run
tm<-as.character(Sys.Date());
fn0<-paste(RCHIVE_HOME, 'source/update/gene.set/UpdateBiosystems.r', sep='/');
fn1<-paste(RCHIVE_HOME, '/source/update/gene.set/log/', tm, '_UpdateBiosystems.r' , sep='');
file.copy(fn0, fn1)

