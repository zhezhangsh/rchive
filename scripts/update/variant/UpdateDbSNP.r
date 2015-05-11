library(devtools);
#source_url("https://raw.githubusercontent.com/zhezhangsh/rchive/master/load.r");
install_github("zhezhangsh/rchive");
library(rchive);

vcfs<-c(
  'GRCh37' = "ftp://ftp.ncbi.nih.gov/snp//organisms/human_9606_b142_GRCh37p13/VCF/All_20150415.vcf.gz",
  'GRCh38' = "ftp://ftp.ncbi.nih.gov/snp//organisms/human_9606_b142_GRCh38/VCF/All_20150416.vcf.gz"
)

t<-sapply(names(vcfs), function(nm) ParseDbSNP(nm, vcfs[nm]));

# Save run
tm<-as.character(Sys.Date());
fn0<-paste(RCHIVE_HOME, 'source/update/variant/UpdateDbSNP.r', sep='/');
fn1<-paste(RCHIVE_HOME, '/source/update/variant/log/', tm, '_UpdateDbSNP.r' , sep='');
file.copy(fn0, fn1)
