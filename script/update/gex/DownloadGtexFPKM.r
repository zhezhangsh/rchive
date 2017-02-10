devtools::install_github("zhezhangsh/rchive");

require(RoCA);
require(rchive);
require(awsomics);
require(DEGandMore);

require(grex); 

path    <- paste(Sys.getenv("RCHIVE_HOME"), 'data/gex/public/gtex', sep='/'); 
url.fkm <- "http://www.gtexportal.org/static/datasets/gtex_analysis_v6p/rna_seq_data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz"
url.smp <- "http://www.gtexportal.org/static/datasets/gtex_analysis_v6/annotations/GTEx_Data_V6_Annotations_SampleAttributesDS.txt"
url.sub <- "http://www.gtexportal.org/static/datasets/gtex_analysis_v6/annotations/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt"

if (!file.exists(path)) dir.create(path, recursive=TRUE);
if (!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'));
if (!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'));

fn.fkm <- paste(path, 'src', TrimPath(url.fkm), sep='/');
if (!file.exists(fn.fkm)) download.file(url.fkm, fn.fkm);

fn.smp <- paste(path, 'src', TrimPath(url.smp), sep='/');
if (!file.exists(fn.smp)) download.file(url.smp, fn.smp);

fn.sub <- paste(path, 'src', TrimPath(url.sub), sep='/');
if (!file.exists(fn.sub)) download.file(url.sub, fn.sub);

fkm <- as.matrix(read.table(fn.fkm, sep='\t', header=TRUE, row.names=1, skip=2)[, -1]);
smp <- read.table(fn.smp, sep='\t', header=TRUE, row.names=1, quote='', stringsAsFactors = FALSE);
sub <- read.table(fn.sub, sep='\t', header=TRUE, row.names=1, quote='', stringsAsFactors = FALSE);

cnm <- colnames(fkm);
cnm <- gsub('\\.', '-', cnm);
smp <- smp[cnm, ];
colnames(fkm) <- cnm; 

smp.nm <- rownames(smp);
sub.id <- paste('GTEX', sapply(strsplit(smp.nm, '-'), function(x) x[2]), sep='-');

age    <- sub[sub.id, 'AGE'];
age    <- paste(sapply(strsplit(age, '-'), function(x) x[1]), 's', sep='')
dth    <- sub[sub.id, 'DTHHRDY'];
gnd    <- sub[sub.id, 'GENDER'];
gnd    <- c('Male', 'Female')[gnd];
tis    <- smp$SMTS;
tis    <- gsub(' ', '_', tis);
tis    <- gsub('_Tissue$', '', tis);
spe    <- smp$SMTSD;
spe    <- sapply(strsplit(spe, ' - '), function(x) x[length(x)]);
spe    <- sapply(strsplit(spe, '\\('), function(x) x[1]);
spe    <- gsub(' ', '_', spe);
spe    <- gsub('_$', '', spe);
smp    <- data.frame(row.names=smp.nm, stringsAsFactors=FALSE, 
                     Tissue=tis, Tissue_Specific=spe, Donor=sub.id, Gender=gnd, Age=age, Death_Hardy_Code=dth, Site=smp$SMCENTER, RIN=smp$SMRIN, 
                     Ischemic_Time=smp$SMTSISCH, Isolation_Batch=smp$SMNABTCH, Isolation_Type=smp$SMNABTCHT, Isolation_Date=smp$SMNABTCHD,
                     Experiment_Batch=smp$SMGEBTCH, Experiment_Type=smp$SMGEBTCHT, Experiment_Date=smp$SMGEBTCHD);


# Map to entrez gene ID
eid <- cleanid(rownames(fkm)); 
ann <- grex(eid);
gid <- ann$entrez_id;
map <- split(rownames(fkm), gid); 
mp1 <- map[sapply(map, length)==1]; 
mp2 <- map[sapply(map, length)>1]; 
ex1 <- fkm[as.vector(unlist(mp1)), ]; 
ex2 <- t(sapply(mp2, function(mp) colSums(fkm[mp, ]))); 
rownames(ex1) <- names(mp1);
rownames(ex2) <- names(mp2); 
ex0 <- rbind(ex1, ex2); 
ord <- order(as.numeric(rownames(ex0))); 
ex0 <- ex0[ord, ]; 

saveRDS(smp, file=paste(path, 'src', 'sample.rds', sep='/'));
saveRDS(fkm, file=paste(path, 'src', 'fpkm.rds', sep='/'));
saveRDS(ex0, file=paste(path, 'src', 'fpkm_entrez.rds', sep='/'));

# Normalization
ex0 <- ex0[rowSums(ex0)>0, ]; 
ex0 <- log2(ex0+1); 
expr <- NormLoess(ex0); 
colnames(expr) <- sub('GTEX-', '', colnames(expr)); 
colnames(expr) <- gsub('-', '_', colnames(expr)); 
saveRDS(expr, file=paste(path, 'r', 'expr.rds', sep='/'));

##############################################################################################################
UpdateLog(smp, paste(Sys.getenv("RCHIVE_HOME"), 'data/gex/public/', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gex/DownloadGtexFPKM.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gex/log/', tm, '_DownloadGtexFPKM.r' , sep='');
file.copy(fn0, fn1)