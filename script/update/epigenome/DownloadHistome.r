devtools::install_github("zhezhangsh/rchive");

require(rchive);
require(RoCA);

options(stringsAsFactors=FALSE);

gn1 <- readRDS('/zhangz/rchive/data/gene/public/entrez/r/human_genes_synonyms2id.rds');
gn2 <- readRDS('/zhangz/rchive/data/gene/public/entrez/r/human_genes_id2symbol.rds');

path <- paste(Sys.getenv("RCHIVE_HOME"), 'data/epigenome/public/histome', sep='/');
path.r <- paste(path, 'r', sep='/');
path.s <- paste(path, 'src', sep='/');
if (!dir.exists(path)) dir.create(path, recursive = TRUE);
if (!dir.exists(path.r)) dir.create(path.r, recursive = TRUE);
if (!dir.exists(path.s)) dir.create(path.s, recursive = TRUE);

writer <- "Arginine deiminases, Arginine methyltransferases, Lysine acetyltransferases, Lysine biotinases, Lysine methyltransferases, Lysine ribosylases, Lysine ubiquitinases, Serine/threonine/tyrosine kinases";
eraser <- "Arginine demethylases, Lysine deacetylases, Lysine demethylases, Lysine deribosylase, Lysine deubiquitinases, Serine/threonine/tyrosine phosphatases";

writer <- strsplit(writer, ', ')[[1]];
eraser <- strsplit(eraser, ', ')[[1]];

writer <- gsub(' ', '_', writer);
eraser <- gsub(' ', '_', eraser); 

fn.w <- sapply(writer, function(x) {
  f1 <- paste('http://www.actrec.gov.in/histome/enzymes.php?enzyme=', x, sep='');
  f2 <- paste(path.s, paste('Writer_', gsub('/', '_', x), '.html', sep=''), sep='/');
  download.file(f1, f2); 
  f2; 
});

fn.r <- sapply(eraser, function(x) {
  f1 <- paste('http://www.actrec.gov.in/histome/enzymes.php?enzyme=', x, sep='');
  f2 <- paste(path.s, paste('Eraser_', gsub('/', '_', x), '.html', sep=''), sep='/');
  download.file(f1, f2); 
  f2; 
});

tbl1 <- lapply(fn.w, function(f) ImportTable(f, rownames = FALSE));
tbl2 <- lapply(fn.r, function(f) ImportTable(f, rownames = FALSE));

tbl1 <- data.frame(Type=rep(names(tbl1), sapply(tbl1, nrow)), do.call('rbind', tbl1), stringsAsFactors = FALSE);
tbl2 <- data.frame(Type=rep(names(tbl2), sapply(tbl2, nrow)), do.call('rbind', tbl2), stringsAsFactors = FALSE);

rownames(tbl1) <- rownames(tbl2) <- NULL;
colnames(tbl1) <- colnames(tbl2) <- c('Type', 'Enzyme', 'Gene', 'Site');

saveRDS(tbl1, paste(path.r, 'writer.rds', sep='/')); 
saveRDS(tbl2, paste(path.r, 'eraser.rds', sep='/')); 

tbl <- data.frame(Function=rep(c('Writer', 'Eraser'), c(nrow(tbl1), nrow(tbl2))), rbind(tbl1, tbl2), stringsAsFactors = FALSE);
sit <- sort(unique(unlist(strsplit(tbl$Site, ', '))));
gns <- lapply(sit, function(s) {
  x <- tbl[grep(s, tbl$Site), , drop=FALSE];
  x <- list(Writer=x[x[,1]=='Writer', 'Gene'], Eraser=x[x[,1]=='Eraser', 'Gene']); 
  x <- lapply(x, function(x) unique(unlist(strsplit(x, ', '))));
  x <- lapply(x, function(x) unique(unlist(gn1[x])));
  x <- lapply(x, function(x) x[!is.na(x) & x!='']); 
  x <- lapply(x, function(x) gn2[x]); 
  x;
});
names(gns) <- sit; 
saveRDS(gns, paste(path.r, 'site.rds', sep='/')); 

##############################################################################################################
# UpdateLog(, paste(Sys.getenv("RCHIVE_HOME"), 'data/epigenome/public/', sep='/'), just.new=FALSE);

tm  <-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0 <-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/epigenome/DownloadHistome.r', sep='');
fn1 <-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/epigenome/log/', tm, '_DownloadHistome.r' , sep='');
file.copy(fn0, fn1, overwrite = TRUE)
