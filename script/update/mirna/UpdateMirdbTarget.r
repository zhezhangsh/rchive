library(devtools);
install_github("zhezhangsh/rchive");
library(rchive);

path <- paste(Sys.getenv("RCHIVE_HOME"), "data/mirna/public/mirdb", sep = "/");

download.file('http://www.mirdb.org/miRDB/download/miRDB_v5.0_prediction_result.txt.gz',
              paste(path, 'src', 'miRDB_v5.0_prediction_result.txt.gz', sep='/'));

spe <- rbind(
  c('human', 'hsa', 'Hs'),
  c('mouse', 'mmu', 'Mm'), 
  c('rat', 'rno', 'Rn'),
  c('chicken', 'gga', 'Gg'),
  c('dog', 'cfa', 'Cf')); 

tbl <- read.csv(paste(path, 'src', 'miRDB_v5.0_prediction_result.txt.gz', sep='/'), sep='\t', header=FALSE, 
                stringsAsFactors = FALSE);

ids <- apply(spe, 1, function(s) {
  cat(s, '\n'); 
  x <- tbl[grep(s[2], tbl[, 1]), , drop=FALSE];
  require(paste('org.', s[3], '.eg.db', sep=''), character.only = TRUE);
  y <- get(paste('org.', s[3], '.egACCNUM', sep=''));
  y <- as.list(y[mappedkeys(y)]);
  z <- unlist(y, use.names=FALSE);
  y <- rep(names(y), sapply(y, length));
  names(y) <- z;
  y <- y[x[, 2]];
  m <- split(as.vector(y), x[, 1]);
  m <- lapply(m, function(m) sort(unique(m[m!=''&!is.na(m)])));
  m <- m[sapply(m, length)>0]; 
  f <- paste(path, 'r', paste(s[1], 'mir2target.rds', sep='_'), sep='/'); 
  saveRDS(m, f); 
  names(m); 
}); 

##############################################################################################################
UpdateLog(ids, paste(Sys.getenv("RCHIVE_HOME"), 'data/mirna/public/mirdb', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/mirna/UpdateMirdbTarget.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/mirna/log/', tm, '_UpdateMirdbTarget.r' , sep='');
file.copy(fn0, fn1, overwrite = TRUE)

