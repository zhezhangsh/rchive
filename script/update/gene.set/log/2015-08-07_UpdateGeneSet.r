# Summarize all gene sets
path.out<-paste(Sys.getenv('RCHIVE_HOME'), 'data/gene.set/r', sep='/');

if (!file.exists(path.out)) dir.create(path.out, recursive=TRUE);

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
  "c2.cp.KEGG" = "C2_KEGG_pathways",
  "c2.cp.REACTOME" = "C2_Reactome_pathways",
  "c3.MIR" = "C3_MicroRNA_targets",
  "c3.TFT" = "C3_TF_targets",
  "c4.CGN" = "C4_Cancer_gene_neighborhoods",
  "c4.CM" = "C4_Cancer_modules",
  "c5.BP" = "C5_GO_biological_processes",
  "c5.CC" = "C5_GO_cellular_components",
  "c5.MF" = "C5_GO_molecular_functions",
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
  data.frame(row.names=id, Collection=rep(nm, nrow(m)), Name=rownames(m), Species='human', Size=m$N, URL=m$URL, stringsAsFactors=FALSE);
}));
list1<-do.call('c', lapply(set1, function(s) s$gene.sets));
if (length(list1) != nrow(meta1)) stop("Error: MSigDB, un-matched length of gene lists and metadata table")
names(list1)<-rownames(meta1);


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
list2<-do.call('c', lapply(set2, function(x) x$Map));
if (length(list2) != length(id)) stop("Error: KEGG, un-matched length of gene lists and metadata table");
n<-sapply(list2, length);
meta2<-data.frame(row.names=paste(sp, id, sep=':'), Collection=tp, Name=nm, Species=sp, Size=n, URL=url, stringsAsFactors=FALSE);
names(list2)<-rownames(meta2);

##############################################################################################################
# BioSystems
path3<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/biosystems/r', sep='/');

src<-c('KEGG', 'BIOCYC', 'REACTOME', 'WikiPathways', 'Pathway-Interaction-Database');
f3<-dir(path3);

# All BioSystems ID to gene lists
full.list<-readRDS(paste(path3, f3[grep('biosystem2gene_list.rds$', f3)][1], sep='/'));

f3<-unlist(lapply(src, function(x) f3[grep(paste('_', x, '.rds', sep=''), f3)]), use.names=FALSE)
f3<-f3[-grep("biosystem2", f3)];
f3<-f3[!(f3 %in% paste('biosystem_', src, '.rds', sep=''))];
f3<-paste(path3, f3, sep='/');

sp<-sapply(strsplit(f3, '_'), function(x) x[2]);
cl<-sub('.rds$', '', sapply(strsplit(f3, '_'), function(x) x[3]));

tbls<-lapply(f3, readRDS);
nm<-lapply(tbls, function(t) paste(t$Accession, t$Name, sep=': '));
meta3<-data.frame(Collection=rep(cl, sapply(tbls, nrow)), Name=unlist(nm, use.names=FALSE), Species=rep(sp, sapply(tbls, nrow)), stringsAsFactors=FALSE);
rownames(meta3)<-unlist(lapply(tbls, rownames), use.names=FALSE);
meta3$Size<-sapply(full.list[rownames(meta3)], length);
meta3$URL<-paste("http://www.ncbi.nlm.nih.gov/biosystems", rownames(meta3), sep='/');
meta3<-meta3[meta3$Size>0, , drop=FALSE];
list3<-full.list[rownames(meta3)];

################################################
# GO lists, without species-specific Biosystem IDs
f3<-dir(path3);
f3<-f3[grep('_GO.rds$', f3)];
f3.mp<-f3[grep('biosystem2gene_', f3)];
f3<-sub('biosystem2gene_', 'biosystem_', f3.mp);
flg<-file.exists(paste(path3, f3, sep='/')) & file.exists(paste(path3, f3.mp, sep='/'));
f3<-f3[flg];
f3.mp<-f3.mp[flg];
sp<-sub('_GO.rds', '', sub('biosystem_', '', f3));

go<-lapply(1:length(sp), function(i) {
  anno<-readRDS(paste(path3, f3[i], sep='/'));
  lst<-readRDS(paste(path3, f3.mp[i], sep='/'));
  anno<-anno[rownames(anno) %in% names(lst), , drop=FALSE];
  lst<-lst[rownames(anno)];
  id<-paste(anno$Accession, sp[i], sep='_');
  tbl<-data.frame(row.names=id, stringsAsFactors=FALSE, Collection='GO', Name=anno$Name, Species=sp[i], Size=sapply(lst, length), 
                  URL=paste("http://www.ebi.ac.uk/QuickGO/GTerm?id=", anno$Accession, sep=''));
  names(lst)<-id;
  list(tbl, lst);
});
go.tbl<-do.call('rbind', lapply(go, function(go) go[[1]]));
go.lst<-do.call('c', lapply(go, function(go) go[[2]]));

meta3<-rbind(meta3, go.tbl);
list3<-c(list3, go.lst);

##############################################################################################################
meta<-list(BioSystems=meta3, KEGG=meta2, MSigDB=meta1);
list<-list(BioSystems=list3, KEGG=list2, MSigDB=list1);
saveRDS(meta, file=paste(path.out, 'metadata.rds', sep='/'));
saveRDS(list, file=paste(path.out, 'all_list.rds', sep='/'));

meta.tree<-lapply(meta, function(m) lapply(split(m, m$Collection), function(x) split(x, x$Species)))
meta.tree<-lapply(meta.tree, function(m) lapply(m, function(m) lapply(m, function(m) m[, c('Name', 'Size', 'URL')])));
saveRDS(meta.tree, file=paste(path.out, 'metadata_as_tree.rds', sep='/'));

sapply(names(list), function(nm) saveRDS(list[[nm]], file=paste(path.out, '/', tolower(nm), '_list.rds', sep='')));


##############################################################################################################
tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/UpdateGeneSet.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/log/', tm, '_UpdateGeneSet.r' , sep='');
file.copy(fn0, fn1)
