library(devtools);
install_github("zhezhangsh/rchive");
library(rchive);

path <- paste(Sys.getenv("RCHIVE_HOME"), "data/mirna/public/mirbase", sep = "/");

download.file('ftp://mirbase.org/pub/mirbase/CURRENT/database_files/mirna.txt.gz',
              paste(path, 'src', 'mirna.txt.gz', sep='/'));
download.file('ftp://mirbase.org/pub/mirbase/CURRENT/database_files/mirna_mature.txt.gz',
              paste(path, 'src', 'mirna_mature.txt.gz', sep='/'));

tbl1 <- read.csv(paste(path, 'src', 'mirna.txt.gz', sep='/'), sep='\t', header = FALSE, 
                row.names = 3, stringsAsFactors = FALSE);
tbl1 <- tbl1[, 2:6];
colnames(tbl1) <- c('Accession', 'Previous_ID', 'Description', 'Sequence', 'Reference');
tbl1$URL <- paste('http://mirbase.org/cgi-bin/mirna_entry.pl?acc=', tbl$Accession, sep='');
saveRDS(tbl1, file=paste(path, 'r', 'anno_mirna_pre.rds', sep='/'));

tbl2 <- read.csv(paste(path, 'src', 'mirna_mature.txt.gz', sep='/'), sep='\t', header = FALSE, 
                 stringsAsFactors = FALSE);
tbl2 <- tbl2[!duplicated(tbl2[, 2]), ];
rownames(tbl2) <- tbl2[, 2];
tbl2 <- tbl2[, c(4, 3, 5, 7, 6)];
colnames(tbl2) <- c('Accession', 'Previous_ID', 'Evidence', 'Similarity', 'Technology');
tbl2$URL <- paste('http://mirbase.org/cgi-bin/mature.pl?mature_acc=', tbl2$Accession, sep='');
tbl2 <- cbind(Name=paste('Mature sequence', rownames(tbl2)), tbl2);
saveRDS(tbl1, file=paste(path, 'r', 'anno_mirna_mature.rds', sep='/'));

##############################################################################################################
UpdateLog(ids, paste(Sys.getenv("RCHIVE_HOME"), 'data/mirna/public/mirbase', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/mirna/UpdateMirbase.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/mirna/log/', tm, '_UpdateMirbase.r' , sep='');
file.copy(fn0, fn1, overwrite = TRUE)

