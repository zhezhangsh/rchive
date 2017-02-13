# Summarize all gene sets
path.out<-paste(Sys.getenv('RCHIVE_HOME'), 'data/gene.set/r', sep='/');

if (!file.exists(path.out)) dir.create(path.out, recursive=TRUE);

##############################################################################################################
##############################################################################################################
{ # MSigDB
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
}; 

##############################################################################################################
##############################################################################################################
{ # KEGG
path2<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/kegg/r', sep='/');
f2 <- dir(path2);
f2 <- f2[grep('2gene.rds$', f2)];
sp <- sapply(strsplit(f2, '_'), function(x) x[1]);
tp <- sub('2gene.rds$', '', sapply(strsplit(f2, '_'), function(x) x[length(x)]));
f2 <- paste(path2, f2, sep='/');
set2 <- lapply(f2, readRDS);
all2 <- lapply(1:length(set2), function(i) {
  x <- set2[[i]]; 
  y <- sp[i];
  z <- tp[i];
  f <- paste(path2, paste('anno_', z, '.rds', sep=''), sep='/');
  a <- readRDS(f); 
  a <- a[rownames(a) %in% names(x), , drop=FALSE]; 
  x <- x[names(x) %in% rownames(a)]; 
  if (z!='pathway' & z!='module') rownames(a) <- names(x) <- paste(y, names(x), sep=':'); 
  cat(y, z, length(x), '\n'); 
  d <- data.frame(Collection=z, Name=a$Name, Species=y, Size=sapply(x, length), URL=a$URL, stringsAsFactors = FALSE);
  list(meta=d, list=x); 
});
meta2 <- lapply(all2, function(x) x[[1]]);
list2 <- lapply(all2, function(x) x[[2]]);
meta2 <- do.call('rbind', meta2); 
list2 <- do.call('c', list2); 
}

##############################################################################################################
##############################################################################################################
{ # BioSystems, species specific
path3<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/biosystems/r', sep='/');

src<-c('KEGG', 'BIOCYC', 'REACTOME', 'WikiPathways', 'Pathway-Interaction-Database');
f3<-dir(path3);

# All BioSystems ID to gene lists
full.list<-readRDS(paste(path3, f3[grep('biosystem2gene_list.rds$', f3)][1], sep='/'));

f3 <- unlist(lapply(src, function(x) f3[grep(paste('_', x, '.rds', sep=''), f3)]), use.names=FALSE)
f3 <- f3[-grep("biosystem2", f3)];
f3 <- f3[!(f3 %in% paste('biosystem_', src, '.rds', sep=''))];
f3 <- paste(path3, f3, sep='/');

sp <- sapply(strsplit(f3, '_'), function(x) x[2]);
cl <- sub('.rds$', '', sapply(strsplit(f3, '_'), function(x) x[3]));

tbls  <- lapply(f3, readRDS);
nm    <- lapply(tbls, function(t) paste(t$Accession, t$Name, sep=': '));
meta3 <- data.frame(Collection=rep(cl, sapply(tbls, nrow)), Name=unlist(nm, use.names=FALSE), Species=rep(sp, sapply(tbls, nrow)), stringsAsFactors=FALSE);
rownames(meta3) <- unlist(lapply(tbls, rownames), use.names=FALSE);
meta3$Size <- sapply(full.list[rownames(meta3)], length);
meta3$URL <- paste("http://www.ncbi.nlm.nih.gov/biosystems", rownames(meta3), sep='/');
meta3 <- meta3[meta3$Size>0, , drop=FALSE];
list3 <- full.list[rownames(meta3)];

x  <- strsplit(meta3[, 2], ': '); 
id <- sapply(x, function(x) x[1]);
nm <- sapply(x, function(x) x[2]); 
rownames(meta3) <- names(list3) <- paste(meta3[, 1], id, sep=':');
meta3[, 2] <- nm; 
}

###########################################
{# GO lists, without species-specific Biosystem IDs
f3 <- dir(path3);
f3 <- f3[grep('_GO.rds$', f3)];
f3.mp <- f3[grep('biosystem2gene_', f3)];
f3 <- sub('biosystem2gene_', 'biosystem_', f3.mp);
flg <- file.exists(paste(path3, f3, sep='/')) & file.exists(paste(path3, f3.mp, sep='/'));
f3 <- f3[flg];
f3.mp <- f3.mp[flg];
sp <- sub('_GO.rds', '', sub('biosystem_', '', f3));

go.root <- c('functional_set'='MF', 'pathway'='BP', 'structural_complex'='CC');

go<-lapply(1:length(sp), function(i) {
  anno<-readRDS(paste(path3, f3[i], sep='/'));
  lst<-readRDS(paste(path3, f3.mp[i], sep='/'));
  anno<-anno[rownames(anno) %in% names(lst), , drop=FALSE];
  lst<-lst[rownames(anno)];
  id<-paste(anno$Accession, sp[i], sep='_');
  tbl<-data.frame(row.names=id, stringsAsFactors=FALSE, Collection=paste('GO', go.root[anno$Type], sep='_'), Name=anno$Name, Species=sp[i], Size=sapply(lst, length), 
                  URL=paste("http://www.ebi.ac.uk/QuickGO/GTerm?id=", anno$Accession, sep=''));
  names(lst)<-id;
  list(tbl, lst);
});
go.tbl <- do.call('rbind', lapply(go, function(go) go[[1]]));
go.lst <- do.call('c', lapply(go, function(go) go[[2]]));
go.tbl <- go.tbl[rownames(go.tbl) %in% names(go.lst), , drop=FALSE];
go.lst <- go.lst[rownames(go.tbl)];
  
meta3<-rbind(meta3, go.tbl);
list3<-c(list3, go.lst);
}

##############################################################################################################
##############################################################################################################
# DisGeNet
path4 <- paste(Sys.getenv("RCHIVE_HOME"), 'data/disease/public/disgenet/r', sep='/');
l <- readRDS(paste(path4, 'disease2gene_full.rds', sep='/'));
d <- readRDS(paste(path4, 'disease_summary.rds', sep='/'));
l <- lapply(l, function(l) l[names(l) %in% rownames(d)]); 
i <- lapply(l, names);
n <- sapply(i, length); 
p <- names(i); 
a <- lapply(i, function(i) d[i, ]); 

meta4 <- data.frame(Collection=rep(p, n), Name=unlist(lapply(a, function(a) a$Name), use.names=FALSE),
                 Species='Human', Size=0, URL=unlist(lapply(a, function(a) a$URL), use.names=FALSE),
                 stringsAsFactors = FALSE);
list4 <- do.call('c', l); 
meta4$Size <- sapply(list4, length); 
id <- paste(unlist(i, use.names=FALSE), rep(names(l), n), sep='_');
id <- sub('^umls:', '', id); 
rownames(meta4) <- names(list4) <- id; 
##############################################################################################################

##############################################################################################################
##############################################################################################################
# Misc.
## OMIM
path5.1 <- paste(Sys.getenv("RCHIVE_HOME"), 'data/disease/public/omim/r', sep='/');
meta5.1 <- readRDS(paste(path5.1, 'phenotype_series.rds', sep='/'));
list5.1 <- readRDS(paste(path5.1, 'phenotype_series2gene.rds', sep='/'));
meta5.1 <- meta5.1[rownames(meta5.1) %in% names(list5.1), , drop=FALSE];
list5.1 <- list5.1[rownames(meta5.1)]; 
meta5.1 <- data.frame(Collection='OMIM_PS', Name=meta5.1$Name, Species='human', Size=sapply(list5.1, length), 
                      URL=meta5.1$URL, stringsAsFactors = FALSE);
rownames(meta5.1) <- names(list5.1) <- paste('OMIM', rownames(meta5.1), sep=':'); 

## 
path5.2 <- paste(Sys.getenv("RCHIVE_HOME"), 'data/literature/public/pubtator/r', sep='/');
list5.2 <- readRDS(paste(path5.2, 'pubmed2gene.rds', sep='/'));
# system.time(ttl <- MapPMID2Title(names(list5.2)));

##############################################################################################################


##############################################################################################################
meta<-list(BioSystems=meta3, KEGG=meta2, MSigDB=meta1, DisGeNET=meta4);
list<-list(BioSystems=list3, KEGG=list2, MSigDB=list1, DisGeNET=list4);
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
