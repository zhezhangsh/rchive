library(devtools);
install_github("zhezhangsh/rchive");
library(rchive);

options(stringsAsFactors=FALSE);

# Types of KEGG databases
types<-c('pathway', 'module', 'compound', 'glycan', 'reaction', 'enzyme', 'disease', 'drug', 'dgroup', 'environ');

# Species to map
species<-c('human'='hsa', 'mouse'='mmu', 'rat'='rno', 'fly'='dme', 'worm'='cel', 'zebrafish'='dre', 'yeast'='sce', 'ecoli'='eco');

# Path to output files
out<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/kegg', sep='/'); 

# Download data and save results
ids<-MapKegg2Gene(types, species, out); 


##############################################################################################################
UpdateLog(ids, paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/kegg', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/UpdateKegg.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/log/', tm, '_UpdateKegg.r' , sep='');
file.copy(fn0, fn1)
