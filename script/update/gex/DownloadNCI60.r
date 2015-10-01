devtools::install_github("zhezhangsh/rchive");

library(rchive);
library(affy);
library(devtools);

options(stringsAsFactors=FALSE);

path<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gex/public/nci60', sep='/');

path.gene<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/entrez/r', sep='/');


if (!file.exists(path)) dir.create(path, recursive=TRUE);
if (!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'));
if (!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'));

gse<-'GSE5949';
########################################################
path.out<-ParseGSE(gse, paste(path, 'src', sep='/'));
########################################################

# CEL files
fn.cel<-dir(paste(path.out, gse, sep='/'));
fn.cel<-fn.cel[grep('CEL.gz$', fn.cel)];
gsm<-sub('.CEL.gz$', '', fn.cel);
fn.cel<-paste(path.out, gse, fn.cel, sep='/');

# Sample metadata
smp<-lapply(gsm, function(s) getGEO(s, GSEMatrix = FALSE, AnnotGPL = FALSE, getGPL = FALSE))
meta<-lapply(smp, Meta);
saveRDS(meta, file=paste(path, 'r', 'original_metadata.rds', sep='/'));
ttl<-sapply(meta, function(x) x$title);
plt<-sapply(strsplit(ttl, ' '), function(x) x[length(x)]);
nms<-sapply(strsplit(ttl, ' '), function(x) x[1]);
names(nms)<-gsm;
names(meta)<-nms;
meta<-meta[!duplicated(nms)];
ch<-lapply(meta, function(x) x$characteristics_ch1);
ch<-lapply(ch, function(ch) {
  v<-sapply(strsplit(ch, ':'), function(x) x[2]);
  v<-gsub("\\s+", " ", v);
  v<-sub('^\\s+', '', v);
  names(v)<-sapply(strsplit(ch, ': '), function(x) x[1]);
  v;
});
cnm<-unique(unlist(lapply(ch, names), use.names=FALSE));
tbl<-matrix('', nr=length(ch), nc=length(cnm), dimnames=list(names(meta), cnm));
for (i in 1:length(ch)) tbl[i, names(ch[[i]])]<-as.vector(ch[[i]]);
fns<-split(fn.cel, plt);
gsm.id<-sub('.CEL.gz', '', sapply(strsplit(fns[[1]], '/'), function(x) x[length(x)]));
tbl<-cbind(Name=rownames(tbl), tbl);
tbl<-tbl[nms[gsm.id], ];
rownames(tbl)<-gsm.id;

# Process expression data matrix
cdf.old<-sapply(fns, function(f) affyio::read.celfile.header(f)[[1]]);
cdf.new<-sapply(cdf.old, InstallBrainarray);
dim.int<-lapply(fns, function(f) affyio::read.celfile.header(f)[[2]]);
names(fns)<-cdf.new;
raws<-lapply(1:5, function(i) {
  cat('Loading', cdf.old[i], '\n');
  exprs<-affyio::read_abatch(fns[[i]], FALSE, FALSE, FALSE, cdf.old[i], dim.int[[i]], TRUE);
  new("AffyBatch", exprs = exprs, cdfName = cdf.new[i], nrow = dim.int[[i]][2], ncol = dim.int[[i]][1],
      annotation = cleancdfname(cdf.new[i], addcdf = FALSE));
});
exprs<-lapply(raws, function(raw) exprs(rma(raw)));
exprs<-lapply(exprs, function(e) {
  cnm<-colnames(e); 
  cnm<-sapply(strsplit(cnm, '/'), function(x) x[length(x)]);
  cnm<-sub('.CEL.gz', '', cnm);
  colnames(e)<-cnm;
  e<-e[!grepl('^AFFX', rownames(e)), , drop=FALSE];
  rownames(e)<-sub('_at$', '', rownames(e)); 
  colnames(e)<-nms[colnames(e)];
  e<-e[, unique(nms)];
  e;
});
expr<-do.call('rbind', exprs[length(exprs):1]);
expr<-expr[!duplicated(rownames(expr)), ];
expr<-expr[order(as.numeric(rownames(expr))), ];

# Put together outputs
mp<-awsomics::GroupSamples(data.frame(tbl, stringsAsFactors = FALSE), 'OrganismPart');
n.smp<-sapply(mp, nrow);
mp<-lapply(mp, function(mp) awsomics::GroupSamples(mp, 'TargetedCellType'));
names(mp)<-sub("central nervous system", 'cns', names(mp));

ds<-data.frame(row.names=sub('^1', 'D', 10000+1:length(mp)), stringsAsFactors = FALSE, Name=names(mp), 
               Num_Gene=nrow(expr), Num_Group=sapply(mp, length), Num_Sample=n.smp, Species='Human');
grp.ds<-rep(rownames(ds), sapply(mp, length));
grp.nm<-gsub(' ', '_', unlist(lapply(mp, names), use.names=FALSE));
grp<-data.frame(row.names=sub('^1', 'G', 10000+1:length(grp.nm)), stringsAsFactors = FALSE, Name=grp.nm,
                Dataset=grp.ds, Organ=ds[grp.ds, 'Name'], Num_Sample=unlist(lapply(mp, function(mp) sapply(mp, nrow)), use.names=FALSE));

smp.grp<-rep(rownames(grp), unlist(lapply(mp, function(mp) sapply(mp, nrow))));
smp.ds<-rep(rownames(ds), n.smp);
smp<-data.frame(row.names=sub('^1', 'S', 10000+1:length(smp.ds)), stringsAsFactors = FALSE, Name=tbl[, 'Name'], 
                GEO=rownames(tbl), Group=smp.grp, Dataset=smp.ds, Organ=tbl[, 'OrganismPart'], Cell=tbl[, 'TargetedCellType'], 
                Disease=tbl[, 'DiseaseState'], Sex=tbl[, 'Sex'], Age=as.numeric(tbl[, 'Age']), Treatment=tbl[, 'Prior Treatment'],
                MDR_Function=tbl[, 'MDR Function'], p53_Status=tbl[, 'p53 Status']);
smp.brow<-smp;
smp.brow[, 'GEO']<-awsomics::AddHref(smp[, 'GEO'], paste("ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", smp[, 'GEO'], sep=''));

meta<-list(Dataset=ds, Group=grp, Sample=smp);

brow<-meta;
brow[[3]]<-smp.brow;
brow<-lapply(brow, function(x) data.frame(ID=rownames(x), x, stringsAsFactors = FALSE));

expr<-expr[, smp[, 'Name'], drop=FALSE];
colnames(expr)<-rownames(smp);
gex<-lapply(rownames(ds), function(id) expr[, rownames(smp)[smp$Dataset==id], drop=FALSE]);
names(gex)<-rownames(ds);

saveRDS(gex, file=paste(path, 'r', 'gex.rds', sep='/'));
saveRDS(meta, file=paste(path, 'r', 'metadata.rds', sep='/'));
saveRDS(brow, file=paste(path, 'r', 'browse_table.rds', sep='/'));

##############################################################################################################
UpdateLog(meta, paste(Sys.getenv("RCHIVE_HOME"), 'data/gex/public/', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gex/DownloadNCI60.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gex/log/', tm, '_DownloadNCI60.r' , sep='');
file.copy(fn0, fn1)
