# Download gene-level read counts of RNA-seq data
q <- c("Assay"='RNA-seq', "File.format"='tsv', "Assembly"='GRCh38', "Output.type"='gene quantifications', "Project"='ENCODE', 
       "Library.made.from"='RNA', "Library.size.range"='>200', "Library.depleted.in"='rRNA')
p <- paste(Sys.getenv("RCHIVE_HOME"), 'data/encode/public/rnaseq/src/grch38/gene', sep='/'); 
o <- paste(Sys.getenv("RCHIVE_HOME"), 'data/encode/public/rnaseq/r/grch38/gene', sep='/'); 

if (!file.exists(p)) dir.create(p, recursive = TRUE);
if (!file.exists(o)) dir.create(o, recursive = TRUE);

devtools::install_github("zhezhangsh/rchive");

require(rchive);
require(awsomics);
require(RoCA);
require(DEGandMore);

download.all <- FALSE;

# Download metadata table
url.meta <- "https://www.encodeproject.org/metadata/type=Experiment/metadata.tsv";
fnt.meta  <- paste(Sys.getenv("RCHIVE_HOME"), 'data/encode/src/metadata.tsv', sep='/'); 
fnr.meta  <- paste(Sys.getenv("RCHIVE_HOME"), 'data/encode/r/metadata.tsv', sep='/'); 

if (download.all | !file.exists(fnr.meta)) {
  f <- TrimPath(url.meta);
  if (file.exists(f)) file.remove(f); 
  system(paste('wget', url.meta)); 
  file.rename(f, fnt.meta); 
  tbl <- read.csv(fnt.meta, sep='\t', stringsAsFactors = FALSE); 
  rownames(tbl) <- tbl[[1]]; 
  tbl <- tbl[, -1]; 
  saveRDS(tbl, fnr.meta); 
} else tbl <-readRDS(fnr.meta); 

sub <- tbl; 
for (i in 1:length(q)) sub <- sub[sub[, names(q)[i]]==q[i], , drop=FALSE]; 

url <- sub$File.download.URL;
fns <- paste(f1, '/', rownames(sub), '.tsv', sep=''); 
for (i in 1:length(url)) {
  if (!file.exists(fns[i])) {
    cat(i, '\n'); 
    system(paste('wget', url[i]));
    file.rename(TrimPath(fns[i]), fns[i]); 
  }
};

##################################
cnt <- MergeEncodeReadCount(fns); 
##################################

##############################################################################################################
UpdateLog(sub, paste(Sys.getenv("RCHIVE_HOME"), 'data/encode/public/', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/encode/DownloadEncode.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/encode/log/', tm, '_DownloadEncode.r' , sep='');
file.copy(fn0, fn1)