require(rchive);
require(RCurl); 

# Download and parse the MSigDb database
path<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/trrust', sep='/');
dir.create(path, showWarnings=FALSE);
dir.create(paste(path, 'r', sep='/'), showWarnings=FALSE);
dir.create(paste(path, 'src', sep='/'), showWarnings=FALSE);

sym <- readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/entrez/r/human_genes_symbol2id.rds', sep='/'));

tbl <- read.csv('http://www.grnpedia.org/trrust/trrust_rawdata.txt', sep='\t', header = FALSE, stringsAsFactors = FALSE);
names(tbl) <- c('TF', 'Gene', 'Regulation', "PubMed"); 

mp <- sym[tbl[, 2]];
mp <- split(as.vector(mp[!is.na(mp)]), tbl[!is.na(mp), 1]); 
nm <- paste('Genes regulated by', names(mp));
id <- as.vector(sym[names(mp)]); 
id[is.na(id)] <- names(mp)[is.na(id)]; 
an <- data.frame(Name=nm, URL=awsomics::UrlEntrezGene(id), stringsAsFactors = FALSE); 
rownames(an) <- names(mp); 

saveRDS(mp, paste(path, 'r', 'tf2gene.rds', sep='/'));
saveRDS(an, paste(path, 'r', 'tf_anno.rds', sep='/'));

##############################################################################################################
UpdateLog(out, paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/trrust', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/UpdateTrrust.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/log/', tm, '_UpdateTrrust.r' , sep='');
file.copy(fn0, fn1, overwrite = TRUE)

