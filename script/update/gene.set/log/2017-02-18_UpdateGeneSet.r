# Summarize all gene sets
path.out<-paste(Sys.getenv('RCHIVE_HOME'), 'data/gene.set/collection', sep='/');
if (!file.exists(path.out)) dir.create(path.out, recursive=TRUE);

species <- c("human", "mouse", "rat", "chimp", "pig", "chicken", "dog", "cow", "worm", "fly", "zebrafish", "yeast", "ecoli"); 

##############################################################################################################
##############################################################################################################
## PubTator
path <- paste(Sys.getenv("RCHIVE_HOME"), 'data/literature/public/pubtator/r', sep='/');
fns <- dir(path); 
fns <- fns[grep('^pubmed2gene', fns)]; 
fns <- fns[-grep('_all_', fns)]; 
fns <- lapply(species, function(s) {
  p <- paste('^pubmed2gene_', s, '_', sep='');
  f <- fns[grep(p, fns)];
  c <- sub(p, '', sub('.rds$', '', f));
  names(f) <- c;
  f;
});
names(fns) <- species;
cll <- names(fns[['human']]); 
fns <- lapply(fns, function(f) f[names(f) %in% cll]); 
fns <- fns[sapply(fns, length)>0]; 

ttl <- readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/literature/public/pubtator/r/pmid2title.rds', sep='/'));
all <- lapply(names(fns), function(nm1) {
  cat(nm1, '\n'); 
  fn <- fns[[nm1]];
  al <- lapply(names(fn), function(nm2) {
    cat(nm2, '\t'); 
    f <- paste(path, fn[nm2], sep='/'); 
    m <- readRDS(f); 
    t <- ttl[names(m)]; 
    i <- names(m); 
    u <- paste('https://www.ncbi.nlm.nih.gov/pubmed/?term=', i, sep='');
    t[t==''] <- i[t=='']; 
    i <- paste('pmid', i, '_', nm2, '_', nm1, sep='');
    b <- data.frame(Collection=nm2, Name=t, Species=nm1, Size=sapply(m, length), URL=u, stringsAsFactors = FALSE);
    rownames(b) <- names(m) <- i; 
    list(b, m); 
  });
  id <- unlist(lapply(al, function(a) names(a[[2]])), use.names=FALSE);
  mt <- do.call('rbind', lapply(al, function(a) a[[1]]));
  mp <- do.call('c', lapply(al, function(a) a[[2]]));
  names(mp) <- rownames(mt) <- id; 
  list(mt, mp); 
});
names(all) <- names(fns); 
ids <- unlist(lapply(all, function(a) names(a[[2]])), use.names=FALSE); 
meta <- do.call('rbind', lapply(all, function(a) a[[1]])); 
list <- do.call('c', lapply(all, function(a) a[[2]])); 
names(list) <- rownames(meta) <- ids; 
smm <- FormatGenesetCollection('PubTator', path.out, meta, list, species);
n <- meta$Size;
smm <- FormatGenesetCollection('PubTator_5-500', path.out, meta[n>=5, ], list[n>=5], species);

##############################################################################################################
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
smm1 <- FormatGenesetCollection('MSigDB', path.out, meta1, list1, species);

##############################################################################################################
##############################################################################################################
# KEGG
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
smm2 <- FormatGenesetCollection('KEGG', path.out, meta2, list2, species);

##############################################################################################################
##############################################################################################################
# BioSystems, species specific
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

###########################################
# GO lists, without species-specific Biosystem IDs
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
smm3 <- FormatGenesetCollection('BioSystems', path.out, meta3, list3, species);

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
meta4$Species <- tolower(meta4$Species); 
smm4 <- FormatGenesetCollection('DisGeNET', path.out, meta4, list4, species);
##############################################################################################################

##############################################################################################################
##############################################################################################################
# iProClass
path5 <- paste(Sys.getenv("RCHIVE_HOME"), 'data/protein/public/iproclass/r', sep='/');
fn1 <- dir(path5);
fn2 <- fn1[grep('2gene.rds', fn1)];
fn3 <- fn1[grep('^anno', fn1)];
names(fn3) <- sub('.rds', '', sapply(strsplit(fn3, '_'), function(x) x[2]));
fn3[1:length(fn3)] <- paste(path5, fn3, sep='/'); 
all <- lapply(names(fn3), function(nm1) { print(nm1); 
  ann <- readRDS(fn3[nm1]);
  fns <- fn2[grep(paste(nm1, '2gene.rds', sep=''), fn2)];
  spe <- sapply(strsplit(fns, '_'), function(x) x[1]); 
  rownames(ann) <- sub('^PIRSF', 'SF', rownames(ann)); 
  names(fns) <- spe;
  all <- lapply(spe, function(nm2) {
    f <- paste(path5, fns[nm2], sep='/'); print(nm2); 
    m <- readRDS(f);
    a <- ann[rownames(ann) %in% names(m), , drop=FALSE];
    m <- m[rownames(a)]; 
    list(a, m); 
  }); 
  l <- do.call('c', lapply(all, function(x) x[[2]])); 
  m <- lapply(all, function(x) x[[1]]);
  i <- unlist(lapply(m, rownames), use.names=FALSE);
  a <- unlist(lapply(m, function(x) x$Name), use.names=FALSE);
  u <- unlist(lapply(m, function(x) x$URL), use.names=FALSE);
  s <- rep(spe, sapply(m, nrow)); 
  n <- sapply(l, length);
  t <- data.frame(Collection=nm1, Name=a, Species=s, Size=n, URL=u, stringsAsFactors = FALSE);
  if (nm1 != 'omim') i <- paste(i, s, sep='_') else i <- paste('OMIM', i, sep=':'); 
  if (nm1 == 'pubmed') i <- paste('pmid', i, sep=''); 
  rownames(t) <- names(l) <- i;
  list(t, l)
});
meta5 <- lapply(all, function(x) x[[1]]);
list5 <- lapply(all, function(x) x[[2]]);
meta5 <- do.call('rbind', meta5);
list5 <- do.call('c', list5);
smm5 <- FormatGenesetCollection('iProClass', path.out, meta5, list5, species);

##############################################################################################################

##############################################################################################################
##############################################################################################################
# RegulatoryNetworks.org
hs1 <- readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/regulatorynetworks/r/anno_human.rds', sep='/'));
hs2 <- readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/regulatorynetworks/r/mapping_human.rds', sep='/'));
mm1 <- readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/regulatorynetworks/r/anno_mouse.rds', sep='/'));
mm2 <- readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/regulatorynetworks/r/mapping_mouse.rds', sep='/'));
meta7 <- rbind(hs1, mm1);
list7 <- c(hs2, mm2); 
meta7 <- data.frame(Collection=paste(meta7$Cell, 'cell'), Name=sub('targets', 'binding targets', meta7$Name), 
                    Species=rep(c('human', 'mouse'), c(nrow(hs1), nrow(mm1))), Size=sapply(list7, length), 
                    URL=meta7$URL, stringsAsFactors = FALSE);
smm7 <- FormatGenesetCollection('RegulatoryNetworks', path.out, meta7, list7, species);

##############################################################################################################

##############################################################################################################
##############################################################################################################
# interactome
meta8 <- readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/interactome/r/anno_all.rds', sep='/'));
list8 <- readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/interactome/r/interact_all.rds', sep='/'));
meta8 <- meta8[order(rownames(meta8)), ];
list8 <- list8[rownames(meta8)]; 
meta8 <- data.frame(Collection=meta8$Collection, Name=meta8$Name, Species='human', Size=sapply(list8, length),
                    URL=meta8$URL, stringsAsFactors = FALSE);
smm8 <- FormatGenesetCollection('CCSBInteractome', path.out, meta8, list8, species);
##############################################################################################################

##############################################################################################################
##############################################################################################################
# GeneSetDB
list9 <- list(
  human = readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/genesetdb/r/human.rds', sep='/')),
  mouse = readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/genesetdb/r/mouse.rds', sep='/')),
  rat = readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/genesetdb/r/rat.rds', sep='/'))
);
url <- readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/genesetdb/r/link.rds', sep='/'))
url <- lapply(list9, function(x) rep(url[names(x)], sapply(x, length)));
url <- as.vector(unlist(url));
nms <- lapply(list9, function(x) lapply(x, names));
cll <- lapply(nms, function(x) rep(names(x), sapply(x, length)));
ids <- lapply(list9, function(x) lapply(names(x), function(nm) paste(nm, 1:length(x[[nm]]), sep='_')));
spe <- rep(names(list9), sapply(nms, function(x) length(unlist(x, use.names=FALSE)))); 
nms <- unlist(nms, recursive = TRUE, use.names=FALSE);
ids <- unlist(ids, recursive = TRUE, use.names=FALSE);
cll <- unlist(cll, recursive = TRUE, use.names=FALSE);
ids <- paste(ids, spe, sep='_'); 
list9 <- lapply(list9, function(x) do.call('c', x));
list9 <- do.call('c', list9);
meta9 <- data.frame(Collection=cll, Name=nms, Species=spe, Size=sapply(list9, length), URL=url,
                    stringsAsFactors = FALSE);
rownames(meta9) <- names(list9) <- ids;
smm9 <- FormatGenesetCollection('GeneSetDB', path.out, meta9, list9, species);

##############################################################################################################

##############################################################################################################
##############################################################################################################
# Enrichment maps
lst <- list(
  human = readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/enrichmentmap/r/human_mapping_all.rds', sep='/')),
  mouse = readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/enrichmentmap/r/mouse_mapping_all.rds', sep='/')),
  rat = readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/enrichmentmap/r/rat_mapping_all.rds', sep='/'))
);
ann <- list(
  human = readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/enrichmentmap/r/human_anno_all.rds', sep='/')),
  mouse = readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/enrichmentmap/r/mouse_anno_all.rds', sep='/')),
  rat = readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/enrichmentmap/r/rat_anno_all.rds', sep='/'))
);
spe <- rep(names(ann), sapply(ann, nrow));
ids <- as.vector(unlist(lapply(ann, rownames)));
ids <- paste(ids, spe, sep='_'); 
ann <- do.call('rbind', ann); 
meta0 <- data.frame(Collection=ann$Collection, Name=ann$Name, Species=spe, Size=0, URL=ann$URL, 
                    stringsAsFactors = FALSE);
list0 <- do.call('c', lst); 
names(list0) <- rownames(meta0) <- ids;
smm0 <- FormatGenesetCollection('EnrichmentMap', path.out, meta0, list0, species);
##############################################################################################################

##############################################################################################################
##############################################################################################################
# Misc.
## OMIM
path6.1 <- paste(Sys.getenv("RCHIVE_HOME"), 'data/disease/public/omim/r', sep='/');
meta6.1 <- readRDS(paste(path6.1, 'phenotype_series.rds', sep='/'));
list6.1 <- readRDS(paste(path6.1, 'phenotype_series2gene.rds', sep='/'));q
meta6.1 <- meta6.1[rownames(meta6.1) %in% names(list6.1), , drop=FALSE];
list6.1 <- list6.1[rownames(meta6.1)]; 
meta6.1 <- data.frame(Collection='OMIM_PS', Name=meta6.1$Name, Species='human', Size=sapply(list6.1, length), 
                      URL=meta6.1$URL, stringsAsFactors = FALSE);
rownames(meta6.1) <- names(list6.1) <- paste('OMIM', rownames(meta6.1), sep=':'); 

## miRDB
path6.2 <- paste(Sys.getenv("RCHIVE_HOME"), 'data/mirna/public/mirdb/r', sep='/');
f <- dir(path6.2);
f <- f[grep('mir2target.rds$', f)]; 
s <- sapply(strsplit(f, '_'), function(x) x[1]);
f <- paste(path6.2, f, sep='/');
names(f) <- s;
list6.2 <- lapply(f, readRDS);
ids <- lapply(list6.2, names);
ann <- readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/mirna/public/mirbase/r/anno_mirna_mature.rds', sep='/'));
meta6.2 <- lapply(ids, function(i) ann[rownames(ann) %in% i, c('Name', 'URL'), drop=FALSE]);
list6.2 <- lapply(1:length(meta6.2), function(i) list6.2[[i]][rownames(meta6.2[[i]])])
ids <- paste('miRDB', unlist(ids, use.names=FALSE), sep=':'); 
meta6.2 <- data.frame(Collection='miRDB', Name=unlist(lapply(meta6.2, function(x) x$Name), use.names=FALSE), 
                      Species=rep(s, sapply(meta6.2, nrow)), Size=0, 
                      URL=unlist(lapply(meta6.2, function(x) x$URL), use.names=FALSE), stringsAsFactors = FALSE);
list6.2 <- do.call('c', list6.2); 
names(list6.2) <- rownames(meta6.2) <- ids; 
meta6.2$Size <- sapply(list6.2, length);
meta6.2 <- meta6.2[meta6.2$Size>0, , drop=FALSE];
list6.2 <- list6.2[rownames(meta6.2)]; 

## ENCODE TFBS
path6.3 <- paste(Sys.getenv("RCHIVE_HOME"), 'data/encode/public/tfbs/r', sep='/');
meta6.3 <- readRDS(paste(path6.3, 'tf_anno.rds', sep='/'));
list6.3 <- readRDS(paste(path6.3, 'tf_targets.rds', sep='/'));
meta6.3 <- meta6.3[rownames(meta6.3) %in% names(list6.3), , drop=FALSE];
list6.3 <- list6.3[rownames(meta6.3)];
meta6.3 <- data.frame(Collection='ENCODE_TFBS', Name=meta6.3$Name, Species='human', Size=sapply(list6.3, length), 
                      URL=meta6.3$URL, stringsAsFactors = FALSE);
rownames(meta6.3) <- names(list6.3) <- paste('ENCODE_TF', names(list6.3), sep=':'); 

##TRRUST
meta6.4 <- readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/trrust/r/tf_anno.rds', sep='/'));
list6.4 <- readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/trrust/r/tf2gene.rds', sep='/'));
meta6.4 <- data.frame(Collection='TRRUST', Name=meta6.4$Name, Species='human', Size=sapply(list6.4, length), 
                      URL=meta6.4$URL, stringsAsFactors = FALSE); 
rownames(meta6.4) <- names(list6.4) <- paste('TRRUST', names(list6.4), sep=':'); 

##TiGER
list6.5 <- readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/tiger/r/tissue2gene.rds', sep='/'));
meta6.5 <- data.frame(Collection='TiGER', Name=paste('Tissue-specific expression in', names(list6.5)), Species='human', 
                      Size=sapply(list6.5, length), URL='http://bioinfo.wilmer.jhu.edu/tiger/', stringsAsFactors = FALSE);
rownames(meta6.5) <- names(list6.5) <- paste('TiGER', names(list6.5), sep=':'); 

## GeneSigDB
list6.6 <- readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/genesigdb/r/mapping.rds', sep='/'));
meta6.6 <- readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/genesigdb/r/anno.rds', sep='/'));
nm <- meta6.6[, 1];
nm <- gsub(' ', '_', nm); 
sp <- tolower(sapply(strsplit(nm, '_'), function(x) x[1]));
nm <- sub('^Human_', '', nm); 
nm <- sub('^Mouse_', '', nm); 
nm <- sub('^Rat_', '', nm); 
meta6.6 <- data.frame(Collection='GeneSigDB', Name=nm, Species=sp, Size=sapply(list6.6, length),
                      URL=meta6.6$URL, stringsAsFactors = FALSE);
rownames(meta6.6) <- names(list6.6);

meta6 <- rbind(meta6.1, meta6.2, meta6.3, meta6.4, meta6.5, meta6.6);
list6 <- c(list6.1, list6.2, list6.3, list6.4, list6.5, list6.6); 
smm6 <- FormatGenesetCollection('Miscellaneous', path.out, meta6, list6, species);

##############################################################################################################

##############################################################################################################
sss <- c('BioSystems', 'KEGG', 'MSigDB', 'EnrichmentMap', 'GeneSetDB',  
         'iProClass', 'RegulatoryNetworks', 'DisGeNet', 'CCSBInteractome')

meta<-list(BioSystems=meta3, KEGG=meta2, MSigDB=meta1, DisGeNET=meta4, iProClass=meta5);
list<-list(BioSystems=list3, KEGG=list2, MSigDB=list1, DisGeNET=list4, iProClass=list5);
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
file.copy(fn0, fn1, overwrite = TRUE);
