library(devtools);
install_github("zhezhangsh/rchive");
library(rchive);
library(GenomicRanges);

options(stringsAsFactors=FALSE);

##################################################################################################################################
path <- paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/dfam', sep='/');
path.r <- paste(path, 'r', sep='/');
path.s <- paste(path, 'src', sep='/');
if (!dir.exists(path.r)) dir.create(path.r);
if (!dir.exists(path.s)) dir.create(path.s);

species <- c('human'='hg38', 'mouse'='mm10', 'zebrafish'='danRer10', 'fly'='dm6', 'worm'='ce10');

# Non-redundant hits
fn <- sapply(names(species), function(nm) {
  fn <- paste(path.s, paste(species[nm], '_dfam.nrph.hits.gz', sep=''), sep='/');
  cat(fn, '\n');
  tb <- read.csv(fn, sep='\t', comment.char = '#', header=FALSE);
  tb <- tb[tb[,2]!='', ];
  p1 <- pmin(tb[, 11], tb[, 12]);
  p2 <- pmax(tb[, 11], tb[, 12]);
  gr <- GRanges(tb[[1]], IRanges(p1, p2), strand=tb[[10]], Bit_Score=tb[, 4], E_Value=tb[, 5],
                Model_Start=tb[[8]], Model_End=tb[[9]], Envelope_Start=tb[[13]], Envelope_End=tb[[14]]);
  names(gr) <- 1:length(gr);
  gr <- split(gr, tb[[2]]);
  
  nms <- sapply(split(tb[[3]], tb[[2]]), unique);
  dsc <- sapply(split(tb[[17]], tb[[2]]), unique);
  len <- sapply(gr, length); 
  
  an <- data.frame(Name=nms, Num_Nonredundant=len, Description=dsc);
  f1 <- paste(path.r, paste(nm, '_anno.rds', sep=''), sep='/');
  f2 <- paste(path.r, paste(nm, '_nonredundant.rds', sep=''), sep='/');
  rownames(an) <- names(gr); 
  saveRDS(an, f1);
  saveRDS(gr, f2);
  f1; 
}); 

# ln <- paste('gunzip', paste(path.s, paste(species, '_dfam.hits.gz', sep=''), sep='/'));
# writeLines(ln, paste(path.s, 'unzip.sh', sep='/'));
# ##############
# 
# # Split all hits
# fn <- sapply(names(fn), function(nm) {
#   an <- readRDS(fn[nm]);
#   id <- rownames(an);
#   f1 <- paste(path.s, paste(species[nm], '_dfam.hits', sep=''), sep='/');
#   f2 <- lapply(id, function(i) {
#     cat(i); 
#     cmmd <- paste("grep '", i, "' ", f1, " > temp.txt", sep=''); 
#     system(cmmd); 
#     tb <- read.csv('temp.txt', sep='\t', comment.char = '#', header=FALSE);
#     tb <- tb[tb[,2]!='', ];
#     p1 <- pmin(tb[, 11], tb[, 12]);
#     p2 <- pmax(tb[, 11], tb[, 12]);
#     gr <- GRanges(tb[[1]], IRanges(p1, p2), strand=tb[[10]], Bit_Score=tb[, 4], E_Value=tb[, 5],
#                   Model_Start=tb[[8]], Model_End=tb[[9]], Envelope_Start=tb[[13]], Envelope_End=tb[[14]]);
#     names(gr) <- 1:length(gr); 
#     fn <- paste(path.r, '/', nm, '_', i, '.rds', sep='');
#     saveRDS(gr, fn); 
#     fn;
#   });
#   f2;
# });
# if (file.exists('temp.txt')) file.remove('temp.txt'); 

##############################################################################################################
tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), 'source/update/gene/UpdateDfam.r', sep='/');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/update/gene/log/', tm, '_', '_UpdateDfam.r' , sep='');
file.copy(fn0, fn1, overwrite = TRUE)

