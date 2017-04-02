require(rchive); 

# Download and parse the MSigDb database
path<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/biogrid', sep='/');
dir.create(path, showWarnings = FALSE); 
dir.create(paste(path, 'r', sep='/'), showWarnings = FALSE);
dir.create(paste(path, 'src', sep='/'), showWarnings = FALSE);

fn <- paste(path, 'src', 'BIOGRID-ALL-3.4.146.tab2.zip', sep='/');
download.file('https://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.4.146/BIOGRID-ALL-3.4.146.tab2.zip',
              fn);
unzip(fn, exdir = paste(path, 'src', sep='/'))

tbl <- read.csv2(paste(path, 'src', 'BIOGRID-ALL-3.4.146.tab2.txt', sep='/'), sep='/' header = TRUE, comment.char = '');

##############################################################################################################
UpdateLog(out, paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/biogrid', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/UpdateBioGrid.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/log/', tm, '_UpdateBioGrid.r' , sep='');
file.copy(fn0, fn1, overwrite = TRUE)