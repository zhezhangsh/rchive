# Summarize all worm gene sets
path.out<-paste(Sys.getenv('RCHIVE_HOME'), 'data/gene.set/r', sep='/');

if (!file.exists(path.out)) dir.create(path.out, recursive=TRUE);

species<-'worm';
sp.mp<-readRDS(paste(Sys.getenv('RCHIVE_HOME'), 'data/gene/public/wormbase/r/human2worm_full.rds', sep='/'))$list;

meta<-readRDS(file=paste(path.out, 'metadata.rds', sep='/'));

##############################################################################################################
# MSigDB
path1<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/msigdb/r', sep='/');
f1<-dir(path1);
f1<-f1[grep('entrez', f1)];
f1<-paste(path1, f1, sep='/');
collections<-c(
  'h.all' = "C0_Hallmark",
  'c1.all' = "C1_Positional",
  "c2.CGP" = "C2_Chemical_and_genetic_perturbations",
  "c2.cp.BIOCARTA" = "C2_BioCarta_Pathways",
  "c3.MIR" = "C3_MicroRNA_targets",
  "c3.TFT" = "C3_TF_targets",
  "c4.CGN" = "C4_Cancer_gene_neighborhoods",
  "c4.CM" = "C4_Cancer_modules",
  'c6.all' = "C6_Oncogenic_signatures",
  'c7.all' = "C7_Immunologic_signatures"
);

f1<-unlist(sapply(names(collections), function(nm) f1[grep(nm, f1, ignore.case=TRUE)])); 
tp<-gsub('\\.', '_', names(f1));
set1<-lapply(f1, readRDS);
names(set1)<-names(tp)<-collections[names(f1)];
meta1<-do.call('rbind', lapply(names(set1), function(nm) {
  m<-set1[[nm]]$meta;
  id<-paste(tp[nm], 1:nrow(m), sep='_');
  data.frame(row.names=id, Collection=rep(nm, nrow(m)), Name=rownames(m), Species=species, Size=m$N, URL=m$URL, stringsAsFactors=FALSE);
}));
lst1<-do.call('c', lapply(set1, function(s) s$gene.sets));
if (length(lst1) != nrow(meta1)) stop("Error: MSigDB, un-matched length of gene lists and metadata table")
names(lst1)<-rownames(meta1);
meta1<-meta1[, c('Collection', 'Name', 'URL')];

if (species!='human' & exists('sp.mp')) {
  lst1<-lapply(lst1, function(l) unique(unlist(sp.mp[l], use.names=FALSE)));
  n<-sapply(lst1, length);
  lst1<-lst1[n>0]; 
  meta1<-meta1[n>0, , drop=FALSE]; 
}

##############################################################################################################
# KEGG
path2<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/kegg/r', sep='/');
f2<-dir(path2);
f2<-f2[grep('2gene.rds$', f2)];
sp<-sapply(strsplit(f2, '_'), function(x) x[1]);
tp<-sub('2gene.rds$', '', sapply(strsplit(f2, '_'), function(x) x[length(x)]));
f2<-paste(path2, f2, sep='/');
set2<-lapply(f2, readRDS);
nm<-unlist(lapply(set2, function(x) x$Name), use.names=FALSE);
sp<-rep(sp, sapply(set2, function(x) length(x$Name)));
tp<-rep(tp, sapply(set2, function(x) length(x$Name)));
id<-do.call('c', lapply(set2, function(x) names(x$Name)));
url<-paste("http://www.genome.jp/dbget-bin/www_bget?", tp, '+', id, sep='');
lst2<-do.call('c', lapply(set2, function(x) x$Map));
if (length(lst2) != length(id)) stop("Error: KEGG, un-matched length of gene lists and metadata table");
n<-sapply(lst2, length);
meta2<-data.frame(row.names=paste(sp, id, sep=':'), Collection=tp, Name=nm, Species=sp, Size=n, URL=url, stringsAsFactors=FALSE);
names(lst2)<-rownames(meta2);
meta2<-meta2[, c('Collection', 'Name', 'URL')];
meta2<-meta2[grep(paste('^', species, ':', sep=''), rownames(meta2)), , drop=FALSE];
meta2<-meta2[meta2$Collection %in% c('compound', 'dgroup', 'disease', 'drug', 'enzyme', 'module', 'pathway', 'reaction'), , drop=FALSE];
lst2<-lst2[rownames(meta2)];
names(lst2)<-rownames(meta2)<-sub(paste('^', species, ':', sep=''), '', rownames(meta2));
meta2$Collection<-paste('KEGG', meta2$Collection, sep='_');

##############################################################################################################
# BioSystems
path3<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/biosystems/r', sep='/');

f3<-dir(path3);
f3<-f3[grep(species, f3)];
f3<-f3[grep('biosystem2gene', f3)];
f3<-paste(path3, f3, sep='/');
f3<-f3[file.exists(f3)];
lst3<-lapply(f3, readRDS);
lst3<-do.call('c', lst3);

meta3<-meta$BioSystems;

# Metadata of GO terms
path.go<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/biosystems/r', sep='/');
fn.go<-dir(path.go);
fn.go<-fn.go[grep(paste("biosystem_", species, '_GO.rds', sep=''), fn.go)];
fn.go<-paste(path.go, fn.go, sep='/');
meta.go<-readRDS(fn.go);
nm<-paste(meta.go$Accession, meta.go$Name, sep='_');
url<-paste("http://amigo.geneontology.org/amigo/term", meta.go$Accession, sep='/');
ct<-c('functional_set'='GO_MF', 'pathway'='GO_BP', 'structural_complex'='GO_CC');
cl<-ct[meta.go$Type];
meta.go<-data.frame(row.names=rownames(meta.go), stringsAsFactors = FALSE, Collection=cl, Name=nm, URL=url);

meta3<-rbind(meta3[, c('Collection', 'Name', 'URL')], meta.go);
meta3<-meta3[rownames(meta3) %in% names(lst3), , drop=FALSE];
lst3<-lst3[rownames(meta3)];

##############################################################################################################
# OMIM diseases
path.omim<-paste(Sys.getenv("RCHIVE_HOME"), 'data/disease/public/omim/r', sep='/');

meta.omim<-readRDS(paste(path.omim, 'omim.rds', sep='/'));
lst.omim<-readRDS(paste(path.omim, 'map_omim2gene.rds', sep='/'));
meta.omim<-meta.omim[rownames(meta.omim) %in% names(lst.omim), , drop=FALSE];
lst.omim<-lst.omim[rownames(meta.omim)];
url<-paste("http://omim.org/entry", rownames(meta.omim), sep='/');
cll<-paste('OMIM', meta.omim$Type, sep='_');
meta.omim<-data.frame(row.names=rownames(meta.omim), stringsAsFactors = FALSE, Collection=cll, 
                      Name=meta.omim$Title, URL=url);
lst.omim<-lst.omim[rownames(meta.omim)];
names(lst.omim)<-rownames(meta.omim)<-paste('OMIM', rownames(meta.omim), sep='');

if (species!='human' & exists('sp.mp')) {
  lst.omim<-lapply(lst.omim, function(l) unique(unlist(sp.mp[l], use.names=FALSE)));
  n<-sapply(lst.omim, length);
  lst.omim<-lst.omim[n>0]; 
  meta.omim<-meta.omim[n>0, , drop=FALSE]; 
}

##############################################################################################################
# PubTator pubmed
path.pm<-paste(Sys.getenv("RCHIVE_HOME"), 'data/literature/public/pubtator/r/gene2pubmed.rds', sep='/');
lst.pm<-readRDS(path.pm);
sp.gn<-readRDS(paste(Sys.getenv("RCHIVE_HOME"), '/data/gene/public/entrez/r/', species, '_genes_full.rds', sep=''));
id<-rep(names(lst.pm), sapply(lst.pm, length));
gn<-unlist(lst.pm, use.names=FALSE);
id<-id[gn %in% rownames(sp.gn)];
gn<-gn[gn %in% rownames(sp.gn)];
lst.pm<-split(gn, id);
lst.pm<-lapply(lst.pm, unique);
meta.pm<-data.frame(row.names=names(lst.pm), stringsAsFactors = FALSE, Collection='PubMed', 
                    Name=names(lst.pm), URL=paste("http://www.ncbi.nlm.nih.gov/pubmed", names(lst.pm), sep='/'));
lst.pm<-lst.pm[rownames(meta.pm)];
names(lst.pm)<-rownames(meta.pm)<-paste('PMID', rownames(meta.pm), sep='');

###############################
# Combine
meta<-list(MSigDb=meta1, KEGG=meta2, BioSystems=meta3, OMIM=meta.omim, PubTator=meta.pm);
list<-list(lst1, lst2, lst3, lst.omim, lst.pm);
src<-rep(names(meta), sapply(meta, nrow));
list<-do.call('c', list);
n<-sapply(list, length);
meta<-do.call('rbind', meta);
meta<-data.frame(row.names=names(list), stringsAsFactors = FALSE, Source=src, Collection=meta$Collection, 
                 Name=meta$Name, Size=n, URL=meta$URL);

saveRDS(list(meta=meta, list=list), file=paste(path.out, '/default_set_', species, '_full.rds', sep=''));
saveRDS(list(meta=meta[n>=5&n<=1000, ], list=list[n>=5&n<=1000]), file=paste(path.out, '/default_set_', species, '_5-1000.rds', sep=''));

##############################################################################################################
tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/CollectDefaultGeneSet.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/log/', tm, '_CollectDefaultGeneSet.r' , sep='');
file.copy(fn0, fn1)
