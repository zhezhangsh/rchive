require(rchive); 

# Download and parse the MSigDb database
path<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/enrichmentmap', sep='/');
dir.create(path, showWarnings = FALSE); 
dir.create(paste(path, 'r', sep='/'), showWarnings = FALSE);
dir.create(paste(path, 'src', sep='/'), showWarnings = FALSE);

btch <- 'February_01_2017';
if (FALSE) system(paste('wget -r --no-parent', paste('http://download.baderlab.org/EM_Genesets', btch, sep='/'), 
                        '-P', paste(path, 'src', sep='/')));

pth <- paste(path, 'src', 'download.baderlab.org/EM_Genesets', btch, sep='/');
fns <- dir(pth, rec=TRUE);
fns <- fns[grep('.gmt$', fns)]; 
fns <- fns[grep('Entrezgene', fns)]; 
spe <- tolower(sapply(strsplit(fns, '/'), function(x) x[1]));

cll <- as.vector(sapply(strsplit(fns, '/'), function(x) rev(x)[1]));
cll <- sub(btch, '', cll);
cll <- sub('_entrezgene', '', cll, ignore.case = TRUE);
cll <- gsub('__', '_', cll);
cll <- gsub('_.gmt', '.gmt', cll);
cll <- sub('human_', '', cll, ignore.case = TRUE);
cll <- sub('mouse_', '', cll, ignore.case = TRUE);
cll <- sub('rat_', '', cll, ignore.case = TRUE);
cll <- sub('.gmt$', '', cll);
cll <- sub('_with_GO_iea', '_with_iea', cll);
cll <- sub('_no_GO_iea', '', cll); 
cll <- sub('HumanCyc', 'MetaCyc', cll); 
cll <- as.vector(cll); 

fns <- paste(pth, fns, sep='/'); 
names(spe) <- names(fns) <- cll;
names(cll) <- fns;
tbl <- cbind(spe, as.vector(cll), fns); 

cl0 <- cll[-grep('AllPathways', cll)];
cl0 <- cl0[-grep('GOALL', cl0)]; 
cl0 <- cl0[cl0!='DrugBank_smallmolecule']; 
cl0 <- cl0[cl0!='DrugBank_all']; 
tbl <- tbl[tbl[,2] %in% cl0, ]; 

all <- apply(tbl, 1, function(x) {
  print(x[2]); 
  lns <- readLines(x[3]); 
  lns <- strsplit(lns, '\t');
  ids <- sapply(lns, function(x) x[1]); 
  ids <- sapply(strsplit(ids, '%'), function(x) rev(x)[1]);
  ids <- sapply(strsplit(ids, ','), function(x) x[1]); 
  ids <- gsub(' ', '_', ids); 
  if (x[2] == 'Reactome') {
    pre <- c('human'='R-HSA-', 'mouse'='R-MMU-', 'rat'='R-RNO-')[x[1]];
    ids[-grep(pre, ids)] <- paste(pre, ids[-grep(pre, ids)], sep='');
    ids <- sapply(strsplit(ids, '\\.'), function(x) x[1]); 
  };
  nms <- sapply(lns, function(x) x[2]); 
  lst <- lapply(lns, function(x) if (length(x)>2) x[3:length(x)] else NULL); 
  names(lst) <- ids; 
  {
    if (x[2]=='DiseasePhenotypes') url <- paste('http://www.ontobee.org/ontology/HP?iri=http://purl.obolibrary.org/obo', ids, sep='/') else
      if (grepl('^DrugBank_', x[2])) url <- paste('https://www.drugbank.ca/drugs', ids, sep='/') else
        if (grepl('^GO_', x[2])) url <- paste('http://amigo.geneontology.org/amigo/term', ids, sep='/') else 
          if (grepl('MSigdb', x[2])) url <- 'http://software.broadinstitute.org/gsea/msigdb/genesets.jsp' else
            if (x[2] == 'KEGG') url <- paste('http://www.kegg.jp/kegg-bin/show_pathway?', tolower(ids), sep='') else
              if (grepl('Cyc$', x[2])) url <- paste('https://metacyc.org/META/NEW-IMAGE?type=NIL&object=', ids, sep='') else 
                if (x[2] == 'IOB') url <- 'http://www.ibioinformatics.org' else 
                  if (x[2] == 'NCI_Nature') url <- 'https://pid.nci.nih.gov/' else 
                    if (x[2] == 'NetPath') url <- 'http://www.netpath.org/browse' else 
                      if (x[2] == 'Panther') url <- paste('http://www.pantherdb.org/pathway/pathwayList.do?searchType=basic&fieldName=all&organism=all&listType=8&fieldValue=', ids, sep='') else 
                        if (x[2] == 'Reactome') url <- paste('http://www.reactome.org/content/detail', ids, sep='/') else url <- ''
  };
  ann <- data.frame(Name=nms, URL=url, stringsAsFactors=FALSE); 
  rownames(ann) <- ids;
  f <- paste(path, 'r', paste(x[1], '_', x[2], '.rds', sep=''), sep='/'); 
  saveRDS(lst, f); 
  list(anno=ann, mapping=lst);  
});
names(all) <- names(tbl[, 2]);

all <- split(all, tbl[, 1]); 
ids <- lapply(names(all), function(s) {
  ann <- lapply(all[[s]], function(x) x[[1]]); 
  lst <- lapply(all[[s]], function(x) x[[2]]); 
  lst <- do.call('c', lst); 
  ids <- unlist(lapply(ann, rownames), use.names=FALSE); 
  cll <- rep(names(ann), sapply(ann, nrow)); 
  ann <- do.call('rbind', ann); 
  ann <- cbind(Collection=cll, ann); 
  ids <- paste(cll, ids, sep=':'); 
  rownames(ann) <- names(lst) <- ids;
  f1 <- paste(path, 'r', paste(s, '_mapping_all.rds', sep=''), sep='/'); 
  f2 <- paste(path, 'r', paste(s, '_anno_all.rds', sep=''), sep='/'); 
  saveRDS(lst, file=f1);
  saveRDS(ann, file=f2); 
  ids;
}); 

##############################################################################################################
UpdateLog(ids, paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/enrichmentmap', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/UpdateEnrichmentMap.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/log/', tm, '_UpdateEnrichmentMap.r' , sep='');
file.copy(fn0, fn1, overwrite = TRUE);

