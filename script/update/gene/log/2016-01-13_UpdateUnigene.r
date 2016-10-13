library(devtools);
install_github("zhezhangsh/rchive");
library(rchive);
options(stringsAsFactors=FALSE);

# Selected species
ftp.files<-c(
  'human' = 'ftp://ftp.ncbi.nih.gov/repository/UniGene//Homo_sapiens/Hs.data.gz',
  'worm' = 'ftp://ftp.ncbi.nih.gov/repository/UniGene//Caenorhabditis_elegans/Cel.data.gz',
  'mouse' = 'ftp://ftp.ncbi.nih.gov/repository/UniGene//Mus_musculus/Mm.data.gz',
  'rat' = 'ftp://ftp.ncbi.nih.gov/repository/UniGene//Rattus_norvegicus/Rn.data.gz',
  'chimp' = 'ftp://ftp.ncbi.nih.gov/repository/UniGene//Pan_troglodytes/Ptr.data.gz',
  'pig' = 'ftp://ftp.ncbi.nih.gov/repository/UniGene//Sus_scrofa/Ssc.data.gz',
  'chicken' = 'ftp://ftp.ncbi.nih.gov/repository/UniGene//Gallus_gallus/Gga.data.gz',
  'dog' = 'ftp://ftp.ncbi.nih.gov/repository/UniGene//Canis_lupus_familiaris/Cfa.data.gz',
  'cow' = 'ftp://ftp.ncbi.nih.gov/repository/UniGene//Bos_taurus/Bt.data.gz',
  'fly' = 'ftp://ftp.ncbi.nih.gov/repository/UniGene//Drosophila_melanogaster/Dm.data.gz',
  'zebrafish' = 'ftp://ftp.ncbi.nih.gov/repository/UniGene//Danio_rerio/Dr.data.gz'
);

path=paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/unigene', sep='/');
parsed<-lapply(names(ftp.files)[1], function(nm) ParseUnigene(ftp.files[nm], nm, TRUE, path=path)); 
ids<-lapply(parsed, function(x) list(id=rownames(x[[1]]), tissue=colnames(x[[2]]), homolog=nrow(x[[3]])));

##############################################################################################################
UpdateLog(ids, paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/unigene', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene/UpdateUnigene.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene/log/', tm, '_UpdateUnigene.r' , sep='');
file.copy(fn0, fn1); 
