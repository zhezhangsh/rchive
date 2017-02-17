require(rchive);
require(RCurl); 

# Download and parse the MSigDb database
path<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/interactome', sep='/');
dir.create(path, showWarnings=FALSE);
dir.create(paste(path, 'r', sep='/'), showWarnings=FALSE);
dir.create(paste(path, 'src', sep='/'), showWarnings=FALSE);

gns <- readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/entrez/r/human_genes_id2symbol.rds', sep='/'));
sym <- names(gns);
names(sym) <- gns;

url <- c(
  'HI-I-05' = 'http://interactome.dfci.harvard.edu/H_sapiens/download/HI-I-05.tsv', 
  'HI-II-14' = 'http://interactome.dfci.harvard.edu/H_sapiens/download/HI-II-14.tsv', 
  'Venkatesan-09' = 'http://interactome.dfci.harvard.edu/H_sapiens/download/Venkatesan-09.tsv', 
  'Yu-11' = 'http://interactome.dfci.harvard.edu/H_sapiens/download/Yu-11.tsv', 
  'Lit-BM-13' = 'http://interactome.dfci.harvard.edu/H_sapiens/download/Lit-BM-13.tsv'
);
fns <- paste(path, 'src', sapply(strsplit(url, '/'), function(x) rev(x)[1]), sep='/');
for (i in 1:length(fns)) download.file(url[i], fns[i]);

tbl <- lapply(fns, function(f) read.table(f, sep='\t', he=TRUE));
map <- lapply(tbl, function(x) {
  y <- x[, grep('entrez', colnames(x), ignore.case = TRUE)];
  a <- c(y[, 1], y[, 2]);
  b <- c(y[, 2], y[, 1]);
  c <- split(a, b);
  c <- lapply(c, unique); 
  d <- gns[names(c)]; 
  names(c) <- d;
  c[!is.na(d)]; 
});
names(map) <- names(url); 

nm <- unlist(lapply(map, names), use.names=FALSE);
mp <- do.call('c', map); 
mp0 <- split(mp, nm); 
mp0 <- lapply(mp0, function(x) unique(unlist(x, use.names=FALSE)));
map$Combined <- mp0; 
sapply(names(map), function(nm) 
  saveRDS(map[[nm]], paste(path, 'r', paste('interact_', nm, '.rds', sep=''), sep='/')))->x;

nm <- c(nm, names(mp0)); 
cl <- rep(names(map), sapply(map, length));
id <- paste(cl, nm, sep='_');
mp <- do.call('c', map);
url <- awsomics::UrlEntrezGene(sym[nm]);
nm <- paste(paste('Interact with', nm, 'protein'), paste(cl, 'set'), sep=', ');
df <- data.frame(Collection=cl, Name=nm, URL=url, stringsAsFactors = FALSE);
names(mp) <- rownames(df) <- id;

saveRDS(mp, paste(path, 'r', 'interact_all.rds', sep='/'));
saveRDS(df, paste(path, 'r', 'anno_all.rds', sep='/'));

##############################################################################################################
UpdateLog(out, paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/interactome', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/UpdateInteractome.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/log/', tm, '_UpdateInteractome.r' , sep='');
file.copy(fn0, fn1, overwrite = TRUE)

