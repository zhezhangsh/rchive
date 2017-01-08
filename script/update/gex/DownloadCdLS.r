devtools::install_github("zhezhangsh/rchive");

library(rchive);
library(affy);
library(devtools);

options(stringsAsFactors=FALSE);

path<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gex/public/cdls', sep='/');

path.gene<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/entrez/r', sep='/');

path.geodb<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gex/db', sep='/');
if (!file.exists(path.geodb)) dir.create(path.geodb);
fn.geodb<-paste(path.geodb, 'GEOmetadb.sqlite', sep='/');
GEOmetadb::getSQLiteFile(path.geodb);

if (!file.exists(path)) dir.create(path, recursive=TRUE);
if (!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'));
if (!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'));

gse.all<-c("GSE64034", "GSE12408", "GSE64031", "GSE27569", "GSE9613", "GSE54741", "GSE55316", "GSE30920", "GSE46316", "GSE59119");

# Data sets that provide raw CEL files from Affy chips
geo<-c("GSE64034", "GSE12408", "GSE64031", "GSE27569", "GSE54741", "GSE30920");
path1<-sapply(geo, function(id) ProcessGEO(id, paste(path, 'src', sep='/'))); 
names(path1)<-geo;

geo.nosupp<-c("GSE55316", "GSE46316");
path2<-sapply(geo.nosupp, function(id) ProcessGEO(id, paste(path, 'src', sep='/'), FALSE)); 

fn.GSE9613<-paste(path, 'src/GSE9613-GPL570_series_matrix.txt.gz', sep='/');
download.file('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE9nnn/GSE9613/matrix//GSE9613-GPL570_series_matrix.txt.gz', fn.GSE9613);
hd<-scan(fn.GSE9613, what='', sep='\n', flush=TRUE, comment='!', n=1);
sid<-strsplit(hd, '\t')[[1]][-1];
sid<-gsub('"', '', sid);
path3<-ProcessGEO("GSE9613", paste(path, 'src', sep='/'), TRUE, sid);
names(path3)<-"GSE9613";

gse<-getGEO('GSE59119', destdir=path, GSEMatrix = FALSE, AnnotGPL = FALSE, getGPL = FALSE);
gsm<-gse@gsms;
sid<-names(gsm);
sid<-sid[gsm[[1]]@header$library_source == 'transcriptomic'];
path4<-ProcessGEO("GSE59119", paste(path, 'src', sep='/'),  FALSE, sid);
names(path4)<-"GSE59119";

####################################################################################################################
# Non-CEL expression data sets
{
sid<-c("GSM1334027", "GSM1334029", "GSM1334031", "GSM1334033", "GSM1334034", "GSM1334035");
download.file("http://www.dtd.nlm.nih.gov/geo/download/?acc=GSE55316&format=file&file=GSE55316%5Frpkm%2Etxt%2Egz",
              paste(path, 'src', 'raw_GSE55316.txt', sep='/'));
lns<-scan(paste(path, 'src', 'raw_GSE55316.txt', sep='/'), what='', flush=TRUE, sep='\n');
hd<-strsplit(lns[1], '\t')[[1]];
rpkm<-t(sapply(lns[-1], function(ln) as.numeric(strsplit(ln, '\t')[[1]][3:8])));
rownames(rpkm)<-sapply(strsplit(lns[-1], '\t'), function(x) x[1])
rpkm<-rpkm[rowSums(rpkm)>0, , drop=FALSE];
gn<-readRDS(paste(path.gene, 'yeast_genes_full.rds', sep='/'));
gid<-rownames(gn);
names(gid)<-gn[, 2];
gid<-gid[names(gid) %in% rownames(rpkm)];
mp<-split(names(gid), gid);
mp<-sapply(mp, function(mp) if (length(mp)==1) mp[1] else mp[rowSums(rpkm[mp, ])==max(rowSums(rpkm[mp, ]))][1]);
rpkm<-rpkm[mp, ];
rownames(rpkm)<-names(mp); 
colnames(rpkm)<-hd[3:8];
expr<-log2(rpkm+1);
expr<-affy::normalize.loess(expr, log.it=FALSE);
colnames(expr)<-sid;
saveRDS(expr, file=paste(path, 'src', 'expr_GSE55316.rds', sep='/'));
}

{
sid<-c("GSM1128586", "GSM1128587", "GSM1128588", "GSM1128589", "GSM1128590", "GSM1128591");
url<-paste("http://www.ncbi.nlm.nih.gov/geo/download/?acc=", sid, "&format=file&file=", sid,
           "%5FL22%5F16", c(84, 83, 79, 80, 81, 82), "%5FmRNAonly%5Fread%5Fcount%2Etxt%2Egz", sep='');
fn<-paste(path, '/src/raw_', sid, '.txt', sep='')
for (i in 1:length(url)) download.file(url[i], fn[i]);
raw<-lapply(fn, function(f) {
  tbl<-read.table(f, sep='\t', header=FALSE);
  e<-tbl[, 8];
  names(e)<-sub('_mRNA$', '', tbl[, 4]);
  e;
});
rpkm<-sapply(1:6, function(i) raw[[i]][names(raw[[1]])]);
rpkm<-rpkm[rowSums(rpkm)>0, , drop=FALSE];
library(org.Mm.eg.db);
x <- org.Mm.egREFSEQ;
mapped_genes <- mappedkeys(x);
xx <- as.list(x[mapped_genes]);
xx[1:2]
ref<-unlist(xx, use.names=FALSE);
gid<-rep(names(xx), sapply(xx, length));
names(gid)<-ref;
gid<-gid[names(gid) %in% rownames(rpkm)];
mp<-split(names(gid), gid);
xx<-xx[names(xx) %in% names(mp)];
mp<-mp[names(xx)];
mp<-sapply(mp, function(mp) if (length(mp)==1) mp[1] else mp[rowSums(rpkm[mp, ])==max(rowSums(rpkm[mp, ]))][1]);
rpkm<-rpkm[mp, ];
rownames(rpkm)<-names(mp);
rpkm<-rpkm[order(as.numeric(rownames(rpkm))), ];
colnames(rpkm)<-sid;
rpkm<-rpkm[rowSums(rpkm)>0, ];
expr<-log2(rpkm+1);
expr<-affy::normalize.loess(expr, log.it=FALSE);
saveRDS(expr, file=paste(path, 'src', 'expr_GSE46316.rds', sep='/'));
}

{
  sid<-c("GSM1428903", "GSM1428904", "GSM1428905", "GSM1428906", "GSM1428907", "GSM1428908", "GSM1428909", "GSM1428910");
  url<-paste("http://www.ncbi.nlm.nih.gov/geo/download/?acc=", sid, "&format=file&file=", sid,
             "%5F", rep(c('Abrain', 'KO%5FEB', 'Pancreas', 'WT%5FEB'), each=2), 
             '%5Fpool', c(1,2,5,6,1,2,3,4), "%2Egenes%2Efpkm%5Ftracking%2Egz", sep='');
  fn<-paste(path, '/src/raw_', sid, '.txt', sep='')
  for (i in 1:length(url)) download.file(url[i], fn[i]);
  raw<-lapply(fn, function(f) {
    tbl<-read.table(f, sep='\t', header=TRUE);
    e<-tbl[, 'FPKM'];
    names(e)<-tbl[,1];
    e;
  });
  rpkm<-sapply(raw, function(x) x[names(raw[[1]])]);
  x <- org.Mm.egENSEMBL;
  mapped_genes <- mappedkeys(x);
  xx <- as.list(x[mapped_genes]);
  gid<-rep(names(xx), sapply(xx, length));
  names(gid)<-unlist(xx, use.names=FALSE);
  gid<-gid[names(gid) %in% rownames(rpkm)];
  mp<-split(names(gid), gid);
  xx<-xx[names(xx) %in% names(mp)];
  mp<-mp[names(xx)];
  mp<-sapply(mp, function(mp) if (length(mp)==1) mp[1] else mp[rowSums(rpkm[mp, ])==max(rowSums(rpkm[mp, ]))][1]);
  rpkm<-rpkm[mp, ];
  rownames(rpkm)<-names(mp);
  rpkm<-rpkm[order(as.numeric(rownames(rpkm))), ];
  colnames(rpkm)<-sid;
  rpkm<-rpkm[rowSums(rpkm)>0, ];
  expr<-log2(rpkm+1);
  expr<-affy::normalize.loess(expr, log.it=FALSE);
  saveRDS(expr, file=paste(path, 'src', 'expr_GSE59119.rds', sep='/'));
}

####################################################################################################################
paths<-c(path1, path3, path2, path4)[gse.all];
fn.phen<-paste(paths, 'phen.rdata', sep='/');
phen.all<-lapply(fn.phen, function(f) eval(parse(text=load(f))));
phen<-lapply(phen.all, function(x) x[, c(colnames(x)=='title', grep('characteristics_ch1', colnames(x)), grep('description', colnames(x)))]);
phen[[7]]<-phen[[7]][c(1,3,5,7,8,9), ];
names(phen)<-names(paths);

fn.expr<-paste(path, '/src/expr_', names(paths), '.rds', sep='');
expr<-lapply(fn.expr, readRDS);
expr<-lapply(expr, function(g) g[!is.na(as.numeric(rownames(g))), , drop=FALSE]);

smp<-list();
smp[[1]]<-data.frame(row.names=rownames(phen[[1]]), Name=phen[[1]]$description);
ph<-phen[[1]];
sb<-sapply(strsplit(ph[, 3], '[ ,]'), function(x) x[3]);
rp<-sapply(strsplit(ph$title, '[ ,]'), function(x) x[6]);
di<-rep(c('CdLS', 'CHOPS', 'Control'), c(4,4,8));
mt<-rep(c('NIPBL', 'AFF4', 'none'), c(4,4,8));
nm<-paste(di, sb, rp, sep='_');
age<-rep(rep(c(7, 10, 6, 12), each=2), 2);
gnd<-rep(rep(c('female', 'male'), each=2), 4);
smp[[1]]<-data.frame(row.names=rownames(ph), stringsAsFactors=FALSE, Name=paste(di, sb, rp, sep='_'), Disease=sub('CHOP', 'CHOPS', ph[,6]), 
                     Subject=sb, Replicate=rp, Age=as.numeric(age), Gender=gnd, Race="Caucasian", 
                     Gene=mt, Genotype=rep(c('mutation', 'wiletype'), each=8), Cell='fibroblast');

ph<-phen[[2]];
nm<-sapply(strsplit(ph[, 1], ' '), function(x) x[3]);
nm<-sub('\\(', '', nm);
di<-c(rep(c('Control', 'CdLS'), c(17, 16)), 'Control', 'CdLS', 'Roberts', 'Roberts', 'Alagille', 'Alagille');
gn<-rep(c('none', 'NIPBL', 'none', 'NIPBL', 'ESCO2', 'JAG1'), c(17, 16, 1, 1, 2, 2));
mt<-rep(c('wildtype', 'mutation', 'wildtype', 'mutation'), c(17, 16, 1, 5));
nm<-paste(di, c(nm[1:35], 1, 2, 1, 2), sep='_');
di[di=='Control']<-'Healthy control';
di[di!='Healthy control']<-paste(di[di!='Healthy control'], 'proband');
smp[[2]]<-data.frame(row.names=rownames(ph), stringsAsFactors=FALSE, Name=nm, Disease=di, Subject=nm, Age=as.numeric(sub("Age, year: ", '', ph[, 3])), 
                     Gender=tolower(sub("Gender: ", '', ph[, 2])), Race=ph[, 10], Gene=gn, Genotype=mt, Cell='lymphoblastoid');

ph<-phen[[3]];
nm<-sapply(strsplit(ph[, 1], '[ ,]'), function(x) paste(x[2], x[3], x[6], sep='_'));
smp[[3]]<-data.frame(row.names=rownames(ph), stringsAsFactors=FALSE, Name=nm, Disease=ph[, 6], 
                     Subject=sapply(strsplit(ph[, 3], '[ ,]'), function(x) x[3]), Replicate=rep(1:2, 5), 
                     Age=rep(c(6, 12, 11, 11, 8), each=2), Gender=rep(c('female', 'male', 'female', 'female', 'male'), each=2), 
                     Race='Caucasian', Gene=rep(c('AFF4', 'none'), c(4, 6)), Genotype=rep(c('mutation', 'wildtype'), c(4, 6)), Cell='fibroblast');
                     

ph<-phen[[4]];
nm<-sapply(strsplit(ph$title, '[ ,]'), function(x) paste(x[1], x[length(x)-4], x[length(x)], sep='_'));
smp[[4]]<-data.frame(row.names=rownames(ph), stringsAsFactors=FALSE, Name=nm, Age=sub('age: ', '', ph[,4]),
                     Gene=rep(c('none', 'esco2', 'none', 'esco2'), c(4,4,3,4)), 
                     Treatment=rep(c('none', 'knockdown', 'none', 'knockdown'), c(4,4,3,4)), Cell='whole embryos');

ph<-phen[[5]];
nm<-paste(sapply(strsplit(ph[, 2], '[ ,]'), function(x) paste(x[5], x[6], sep='_')), rep(1:2, 3), sep='_');
smp[[5]]<-data.frame(row.names=rownames(ph), stringsAsFactors=FALSE, Name=nm, Gene=rep(c('CTCF', 'hScc1', 'none'), each=2),
                     Treatment=rep(c('CTCF RNAi', 'hScc1 RNAi', 'Control RNAi'), each=2), Cell='Hela');

ph<-phen[[6]];
x<-apply(ph[, -1], 2, function(x) sapply(strsplit(x, ': '), function(x) x[length(x)]));
nm<-as.vector(gsub('\\.', '_', x[, 5]));
gn<-sapply(strsplit(x[, 2], ' mutant'), function(x) x[1]);
gn<-sub(' double', '', gn);
gn[gn=="wild type control"]<-'none';
gn<-sub(' ', '+', gn);
gt<-rep('mutant', length(gn));
gt[gn=='none']<-'wildtype';
x[1,]
smp[[6]]<-data.frame(row.names=rownames(ph), stringsAsFactors=FALSE, Name=nm, Strain=x[, 1], Batch=x[, 4], Gene=gn, Genotype=gt);

ph<-phen[[7]];
x<-apply(ph[, -1], 2, function(x) sapply(strsplit(x, ': '), function(x) x[length(x)]));
x[,2]<-gsub('[- ]', '_', x[, 2]);
smp[[7]]<-data.frame(row.names=rownames(ph), stringsAsFactors=FALSE, Name=paste(x[, 2], rep(c('wildtype', 'mutant'), each=3), x[,1], sep='_'), 
                     Gene=rep(c('none', 'scc2-4'), each=3), Genotype=rep(c('wildtype', 'mutant'), each=3), Treatment=x[,2]);

ph<-phen[[8]];
x<-apply(ph[, -1], 2, function(x) sapply(strsplit(x, ': '), function(x) x[length(x)]));
gn<-sapply(strsplit(ph$title, ' '), function(x) x[2]);
nm<-paste(gn, rep(1:2, 4), sep='_')
smp[[8]]<-data.frame(row.names=rownames(ph), stringsAsFactors=FALSE, Name=nm, Gene=gn, Treatment='knockdown', Cell='V6.5 mESC')

ph<-phen[[9]];
x<-apply(ph[, -1], 2, function(x) sapply(strsplit(x, ': '), function(x) x[length(x)]));
nm<-sub('^RNA[-_]Seq_', '', ph$title);
gn<-sapply(strsplit(nm, '_'), function(x) x[1]);
smp[[9]]<-data.frame(row.names=rownames(ph), stringsAsFactors=FALSE, Name=nm, Gene=gn, Treatment='shRNA', Cell='V6.5 mESC')
smp[[8]]

ph<-phen[[10]];
x<-apply(ph[, -1], 2, function(x) sapply(strsplit(x, ': '), function(x) x[length(x)]));
age<-sub(' old$', '', x[, 3]);
cll<-x[, 4];
trt<-sub(' ', '', sub('^SA1 ', '', x[, 2]));
nm<-paste(sub('^whole ', '', cll), trt, rep(1:2, 4), sep='_');
trt[trt=='wildtype']<-'none';
smp[[10]]<-data.frame(row.names=rownames(ph), stringsAsFactors=FALSE, Name=nm, Strain='C57BL/6', Age=age, Gene=rep(c('none', 'SA1', 'none', 'none'), each=2), Treatment=trt, Cell=cll);

names(smp)<-names(paths);
####################################################################################################################
cnm<-c("Name", "Disease", "Cell", "Strain", "Species", "Gene", "Genotype", "Treatment", "Time_Post", "Subject", "Replicate", "Age", "Gender", "Race", "Batch");      
smp.all<-matrix('', nrow=sum(sapply(smp, nrow)), ncol=length(cnm), dimnames=list(unlist(lapply(smp, rownames)), cnm))
for (i in 1:length(smp)) for (j in 1:ncol(smp[[i]])) smp.all[rownames(smp[[i]]), colnames(smp[[i]])[j]]<-smp[[i]][, j];
smp.all<-smp.all[, colnames(smp.all)!='Cell' & colnames(smp.all)!='Strain'];
s<-lapply(rownames(smp.all), function(s) getGEO(s, GSEMatrix = FALSE, AnnotGPL = FALSE, getGPL=FALSE));
names(s)<-rownames(smp.all);
smp.all<-data.frame(smp.all, stringsAsFactors=FALSE);
smp.all$Source<-sapply(s, function(s) s@header$source_name_ch1)
smp.all$Original_Title<-sapply(s, function(s) s@header$title)

sp.id<-sapply(s, function(s) s@header$taxid_ch1);
sp<-c('Human', 'Zebrafish', 'Yeast', 'Yeast', 'Mouse');
names(sp)<-unique(sp.id);
smp.all[, 'Species']<-sp[sp.id];

####################################################################################################################
grp.all<-list();
grp.all[[1]]<-awsomics::GroupSamples(smp[[1]], 'Disease')[3:1];
grp.all[[2]]<-awsomics::GroupSamples(smp[[2]], 'Disease');
grp.all[[3]]<-awsomics::GroupSamples(smp[[3]], 'Disease')[2:1];
grp.all[[4]]<-split(smp[[4]], sub('_[1-4]$', '', smp[[4]][, 1]));
grp.all[[5]]<-split(smp[[5]], sub('_RNAi_[1-4]$', '', smp[[5]][, 1]))[c(2, 1, 3)];
grp.all[[6]]<-awsomics::GroupSamples(smp[[6]], c('Gene'))[c(6, 1, 3, 4, 8, 2, 5, 7)];
grp.all[[7]]<-awsomics::GroupSamples(smp[[7]], c('Genotype'))[2:1];
grp.all[[8]]<-awsomics::GroupSamples(smp[[8]], c('Gene'))[4:1];
grp.all[[9]]<-awsomics::GroupSamples(smp[[9]], c('Gene'))[c(1, 2, 4, 3)];
grp.all[[10]]<-split(smp[[10]], sub('_[1-4]$', '', smp[[10]][, 1]))[4:1];

geo.db<-dplyr::src_sqlite(fn.geodb);
gse.tbl<-dplyr::tbl(geo.db, 'gse');
gse.tbl<-as.data.frame(dplyr::collect(gse.tbl), stringsAsFactors=FALSE);
rownames(gse.tbl)<-gse.tbl$gse;
gse.tbl<-gse.tbl[names(smp), -1];
ds.nm<-c("CHOPS and CdLS syndrome", "Defective Cohesin in CdLS Mediates Gene Expression", "Transcriptome characterization of CHOPS syndrome",
         "Expression data from zebrafish depleted of Esco2", "Cohesin mediates transcriptional insulation by CCCTC", "Analysis of a yeast eco1 mutant", 
         "Cohesin loader Scc2 regulates development", "Loading and Translocation of Condensin II in mESC", "Positioning of structural maintenance complexes",
         "Contribution of cohesin-SA1 to chromatin and expression");
ds<-data.frame(row.names=sub('^1', 'D', 10000+1:length(smp)), stringsAsFactors=FALSE, Name=ds.nm, GEO=gse.tbl$gse, PubMed=as.character(gse.tbl$pubmed_id),
                Num_Gene=sapply(expr, nrow), Num_Group=sapply(grp.all, length), Num_Sample=sapply(smp, nrow), Species=smp.all[sapply(grp.all, function(x) rownames(x[[1]])[1]), 'Species'],
                Type=sapply(strsplit(gse.tbl$type, ';'), function(x) x[1]), Original_Title=gse.tbl$title);
ds$PubMed[is.na(ds$PubMed)]<-'';
ds[c('D0001', 'D0003'), 'PubMed']<-'25730767'
ds.by.id<-apply(gse.tbl, 1, as.list);
names(ds.by.id)<-rownames(ds);

grp.nm<-unlist(lapply(grp.all, names), use.names=FALSE);
grp.nm<-gsub(' ', '_', grp.nm);
grp<-data.frame(row.names=sub('^1', 'G', 10000+1:length(grp.nm)), stringsAsFactors=FALSE, Name=grp.nm, Dataset=rep(rownames(ds), sapply(grp.all, length)), 
                Dataset_Name=ds[rep(rownames(ds), sapply(grp.all, length)), 'Name'], 
                Num_Sample=unlist(lapply(grp.all, function(g) sapply(g, nrow)), use.names=FALSE));

grp$Disease<-unlist(lapply(grp.all, function(x) lapply(x, function(x) {y<-unique(x$Disease); if(is.null(y)) y<-'none'; y;})));
grp$Gene<-unlist(lapply(grp.all, function(x) lapply(x, function(x) {y<-unique(x$Gene); if(is.null(y)) y<-'none'; y;})));
grp$Genotype<-unlist(lapply(grp.all, function(x) lapply(x, function(x) {y<-unique(x$Genotype); if(is.null(y)) y<-'none'; y;})));
grp$Treatment<-unlist(lapply(grp.all, function(x) lapply(x, function(x) {y<-unique(x$Treatment); if(is.null(y)) y<-'none'; y;})));
grp.by.id<-apply(grp, 1, as.list);
names(grp.by.id)<-rownames(grp);

sm.id<-as.vector(unlist(lapply(grp.all, function(g) lapply(g, rownames))));
sm<-data.frame(row.names=sub('^1', 'S', 10000+1:length(sm.id)), stringsAsFactors = FALSE, Name=smp.all[sm.id, 'Name'], GEO=sm.id);
sm$Group<-rep(rownames(grp), unlist(lapply(grp.all, function(g) sapply(g, nrow))));
sm$Group_Name<-grp[sm$Group, 'Name'];
sm$Dataset<-grp[sm$Group, 'Dataset'];
sm$Dataset_Name<-grp[sm$Group, 'Dataset_Name'];
sm<-data.frame(sm, smp.all[sm$GEO, -1], stringsAsFactors = FALSE);
sm.by.id<-s[sm$GEO];
names(sm.by.id)<-rownames(sm);

for (i in 1:length(expr)) colnames(expr[[i]])<-sapply(strsplit(colnames(expr[[i]]), '_'), function(x) x[1]);
sid<-rownames(sm);
names(sid)<-sm$GEO;
for (i in 1:length(expr)) colnames(expr[[i]])<-as.vector(sid[colnames(expr[[i]])]);

gex<-expr;
names(gex)<-rownames(ds);

meta<-list(Dataset=ds, Group=grp, Sample=sm);
meta.by.id<-list(Dataset=ds.by.id, Group=grp.by.id, Sample=sm.by.id);

brow<-meta;
brow[[1]][, 'GEO']<-awsomics::AddHref(ds[, 'GEO'], 
                                      paste("http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", ds[, 'GEO'], sep=''));
brow[[1]][, 'PubMed']<-awsomics::AddHref(ds[, 'PubMed'], 
                                         paste("http://www.ncbi.nlm.nih.gov/pubmed/?term=", ds[, 'PubMed'], sep=''))
brow[[2]][, 'Dataset']<-awsomics::AddHref(grp[, 'Dataset'], 
                                          paste("http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", ds[grp$Dataset, 'GEO'], sep=''));
brow[[3]][, 'GEO']<-awsomics::AddHref(sm[, 'GEO'], 
                                      paste("http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", sm[, 'GEO'], sep=''));
brow<-lapply(brow, function(t) data.frame(ID=rownames(t), t, stringsAsFactors = FALSE));

saveRDS(gex, file=paste(path, 'r', 'gex.rds', sep='/'));
saveRDS(meta, file=paste(path, 'r', 'metadata.rds', sep='/'));
saveRDS(meta.by.id, file=paste(path, 'r', 'metadata_by_id.rds', sep='/'));
saveRDS(brow, file=paste(path, 'r', 'browse_table.rds', sep='/'));

##############################################################################################################
UpdateLog(paths, paste(Sys.getenv("RCHIVE_HOME"), 'data/gex/public/', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gex/DownloadCdLS.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gex/log/', tm, '_DownloadCdLS.r' , sep='');
file.copy(fn0, fn1)
