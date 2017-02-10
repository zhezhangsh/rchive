library(devtools);
# source_url("https://raw.githubusercontent.com/zhezhangsh/rchive/master/load.r");
install_github("zhezhangsh/rchive");
library(rchive);

options(stringsAsFactors=FALSE);

# download source files from PubTator FTP server
ftp.files<-c(
  chemical="ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator//chemical2pubtator.gz",
  species="ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator//species2pubtator.gz",
  disease="ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator//disease2pubtator.gz",
  gene="ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator//gene2pubtator.gz",
  mutation="ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator//mutation2pubtator.gz"
)
for (i in 1:length(ftp.files)) {
  cat('Loading', ftp.files[i]); 
  ParsePubtator(ftp.files[i]);
}

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), 'source/script/update/literature/UpdatePubtator.r', sep='/');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/literature/log/', tm, '_UpdatePubtator.r' , sep='');
file.copy(fn0, fn1)

