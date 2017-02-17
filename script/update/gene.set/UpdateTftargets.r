# Download and parse the MSigDb database
path<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/tftargets', sep='/');
dir.create(path, showWarnings = FALSE); 
dir.create(paste(path, 'r', sep='/'), showWarnings = FALSE);
dir.create(paste(path, 'src', sep='/'), showWarnings = FALSE);

library(RCurl)
download.file(
  url = "https://raw.githubusercontent.com/slowkow/tftargets/master/data/tftargets.rda",
  destfile = paste(path, 'src', "tftargets.rda", sep='/')
);

nms <- load(paste(path, 'src', "tftargets.rda", sep='/'));
all <- list();
for (i in 1:length(nms)) {
  all[[i]] <- get(nms[i])
};
names(all) <- nms;

lst <- all$ENCODE;


##############################################################################################################
tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/UpdateGeneSet.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/log/', tm, '_UpdateGeneSet.r' , sep='');
file.copy(fn0, fn1, overwrite = TRUE);
