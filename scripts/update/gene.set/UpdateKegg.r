library(devtools);
source_url("https://raw.githubusercontent.com/zhezhangsh/rchive/master/load.r");

# Species to map
species<-c('human'='hsa', 'mouse'='mmu', 'rat'='rno', 'fly'='dme', 'worm'='cel', 'zebrafish'='dre', 'yeast'='sce', 'ecoli'='eco');

# Path to output files
out<-paste(RCHIVE_HOME, 'data/gene.set/public/kegg', sep='/'); 

# Download data and save results
ids<-MapKeggPath2Gene(species, out); 

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(RCHIVE_HOME, 'source/update/gene.set/KeggP2G.r', sep='/');
fn1<-paste(RCHIVE_HOME, '/source/update/gene.set/log/', tm, '_KeggP2G.r' , sep='');
file.copy(fn0, fn1)

