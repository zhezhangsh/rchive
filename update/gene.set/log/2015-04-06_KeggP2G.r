library(devtools);
source_url("https://raw.githubusercontent.com/zhezhangsh/rchive/master/load.r");

# Species to map
species<-c('human'='hsa', 'mouse'='mmu', 'rat'='rno', 'fly'='dme', 'worm'='cel', 'zebrafish'='dre', 'yeast'='sce', 'ecoli'='eco');

# Path to output files
out<-paste(RCHIVE_HOME, 'data/gene.set/public/kegg/r', sep='/'); 

# Download data and save results
ids<-MapKeggPath2Gene(species, out); 

# save log, all pathway IDs
log<-readRDS(paste(out, 'log.rds', sep='/'));
tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
log[[tm]]<-ids;
saveRDS(log, file=paste(out, 'log.rds', sep='/'));

fn0<-paste(RCHIVE_HOME, 'update/gene.set/KeggP2G.r', sep='/');
fn1<-paste(RCHIVE_HOME, '/update/gene.set/log/', tm, '_KeggP2G.r' , sep='');