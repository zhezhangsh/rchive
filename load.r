# Load functions through GitHub repo

RCHIVE_HOME<-"/zhangz/rchive";

library(devtools);

###
source_url("https://raw.githubusercontent.com/zhezhangsh/rchive/master/extra/HtmlUtilities.r");


###
source_url("https://raw.githubusercontent.com/zhezhangsh/rchive/master/assembly/MapChromosome.r");
source_url("https://raw.githubusercontent.com/zhezhangsh/rchive/master/gene/ParseGtf.r");
source_url("https://raw.githubusercontent.com/zhezhangsh/rchive/master/gene/MapGtf.r");
source_url("https://raw.githubusercontent.com/zhezhangsh/rchive/master/gene/ParseEntrez.r");
source_url("https://raw.githubusercontent.com/zhezhangsh/rchive/master/gene.set/MapKegg.r");
source_url("https://raw.githubusercontent.com/zhezhangsh/rchive/master/gene.set/MapBiosystems.r");
source_url("https://raw.githubusercontent.com/zhezhangsh/rchive/master/variant/ParseDbSNP.r");

