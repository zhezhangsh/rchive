# Summarize all gene sets
path.out<-Sys.getenv('RCHIVE_HOME', 'data/gene.set/r', sep='/');

if (!file.exists(path.out)) dir.create(path.out, recursive=TRUE);

##############################################################################################################
# MSigDB
path1<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/msigdb/r', sep='/');
f1<-dir(path1);
f1<-f1[grep('entrez', f1)];
f1<-paste(path.in, f1, sep='/');
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
set1<-lapply(f1, readRDS);
names(set1)<-collections[names(f1)];
meta1<-do.call('rbind', lapply(names(set1), function(nm) {
  m<-set1[[nm]]$meta;
  data.frame(Collection=rep(nm, nrow(m)), Name=rownames(m), Species='human', Size=m$N, URL=m$URL, stringsAsFactors=FALSE);
}));
list1<-do.call('c', lapply(set1, function(s) s$gene.sets));
if (length(list1) != nrow(meta1)) stop("Error: MSigDB, un-matched length of gene lists and metadata table")
names(list1)<-rownames(meta1);


##############################################################################################################
# KEGG
path2<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/kegg/r', sep='/');
f2<-dir(path2);
f2<-f2[grep('2gene.rds$', f2)];
nm<-sapply(strsplit(f2, '_'), function(x) x[1]);
f2<-paste(path2, f2, sep='/');
set2<-lapply(f2, readRDS);
names(set2)<-nm;
meta2<-do.call('rbind', lapply(names(set2), function(nm) {
  m<-set2[[nm]]$Pathway;
  n<-as.vector(sapply(names(m), function(x) length(set2[[nm]]$Pathway2Gene[[x]])));
  url<-paste("http://www.genome.jp/dbget-bin/www_bget?pathway+", names(m), sep='');
  data.frame(Collection=rep(nm, length(m)), Name=paste(names(m), m, sep=': '), Species=rep(nm, length(m)), Size=n, URL=url, stringsAsFactors=FALSE);
}));
list2<-do.call('c', lapply(set2, function(s) s$Pathway2Gene));
if (length(list2) != nrow(meta2)) stop("Error: KEGG, un-matched length of gene lists and metadata table")
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
meta3$URL<-paste("http://www.ncbi.nlm.nih.gov/biosystems/", rownames(meta3), sep='/');
meta3<-meta3[meta3$Size>0, , drop=FALSE];
list3<-full.list[rownames(meta3)];
