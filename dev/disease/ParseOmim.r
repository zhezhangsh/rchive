# Parse OMIM database
ParseOmim<-function(url="ftp://ftp.omim.org/OMIM", path=paste(Sys.getenv("RCHIVE_HOME"), 'data/disease/public/omim', sep='/'), download.all=FALSE) {

  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  if(!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'), recursive=TRUE);
  if(!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'), recursive=TRUE);
 
  # Download from OMIM ftp site
  fn.src<-c("genemap", "genemap.key", "genemap2.txt", "mim2gene.txt", "morbidmap", "omim.txt.Z");
  fn.ftp<-paste(url, fn.src, sep='/');
  fn.loc<-paste(path, 'src', fn.src, sep='/');
  names(fn.loc)<-fn.src;
  for (i in 1:length(fn.ftp)) if (download.all | !file.exists(fn.loc[i])) download.file(fn.ftp[i], fn.loc[i]);
  
  ############################################################################
  ## Process full record of OMIM entries
  # Load file
  fn.omim<-fn.loc['omim.txt.Z'];
  cd<-system(paste('uncompress -kf', fn.omim));
  lns<-scan(sub('.Z$', '', fn.omim), what='', flush=TRUE, sep='\n');
  
  # assign entry ID to lines
  ind1<-grep('^\\*RECORD\\*', lns); # every new OMIM entry
  lns[ind1]<-'';
  ent<-rep(1:length(ind1), c(ind1[-1]-ind1[-length(ind1)], length(lns)-ind1[length(ind1)]+1)); 
  ent[ind1]<-'';
  ent.id<-as.character(1:length(ind1));
  
  # split into individual fields
  ind2<-grep('^\\*FIELD\\*', lns); # every fields
  typ<-sub('^\\*FIELD\\* ', '', lns[ind2]); # field type
  fld<-c('', rep(1:length(ind2), c(ind2[-1]-ind2[-length(ind2)], length(lns)-ind2[length(ind2)]+1))); 
  fld[c(ind1, ind2)]<-''; # exclude separator lines
  grp<-split(lns, fld)[-1][as.character(1:length(typ))]; # group lines into individual fields
  names(grp)<-ent[ind2];
  grp<-split(grp, typ); # group by field type
  
  ##########################################################################
  ## Retrieve individual fields
  # Unique OMIM ID
  id<-sapply(grp[['NO']], function(x) x[[1]][1]);
  n<-sapply(grp[['NO']], length);
  if (max(n) > 1) warning("Some entries have more than one OMIM ID\n");
  if (length(id) != length(ind1)) warning("Not all entry have OMIM ID\n");
  id<-id[ent.id];
  id[is.na(id)]<-'';
  names(id)<-ent.id
  for (i in 1:length(grp)) names(grp[[i]])<-id[names(grp[[i]])];
  
  # OMIM title
  title<-sapply(grp[['TI']], function(x) paste(x, collapse=' ')); 
  title<-sapply(split(as.vector(title), names(title)), function(x) paste(x, collapse='; '));
  title<-title[id];
  names(title)<-id;
  title[is.na(title)]<-'';
  title<-gsub(' ;;', ';;', gsub(';; ', ';;', title));
  title<-strsplit(title, ';;');
  title.alt<-lapply(title, function(x) x[-1]);
  title<-sapply(title, function(x) x[1]);
  title.id<-sapply(strsplit(title, ' '), function(x) x[1]);
  title<-sapply(title, function(x) substr(x, gregexpr(' ', x)[[1]][1]+1, nchar(x)));
  
  # Get entry status by first character
  code<-c(
    "Other, mainly phenotypes with suspected mendelian basis",
    '%' = "Phenotype description or locus, molecular basis unknown", 
    '#' = "Phenotype description, molecular basis known",
    '^' = "Record moved",
    '*' = "Gene description",
    '+' = "Gene and phenotype, combined"
  );
  status<-code[substr(title.id, 1, 1)];
  status[is.na(status)]<-code[1];
  names(status)<-names(title);

  # text
  tx<-sapply(grp[['TX']], function(x) paste(x, collapse=' '));
  tx<-split(tx, names(tx));
  
  # clinical synopsis
  cs<-lapply(names(grp[['CS']]), function(nm) { 
    #x<-unlist(strsplit(grp[['CS']][[nm]], ':'), use.names=FALSE);
    x<-sub(' $', '', sub('^ ', '', gsub('[ ]+', ' ', gsub(';$', '', grp[['CS']][[nm]]))));
    x[grep(':', x)]<-paste('\n', x[grep(':', x)], sep='');
    x<-sub('[:;.]$', '', x);
    x<-sub(': ', '; ', x);
    x<-paste(x, collapse=';');
    x<-gsub('; ', ';', x);
    x<-gsub(';;', ';', x);
    x<-sub('^\n', '', x);
    x<-strsplit(x, '\n')[[1]];
    x<-strsplit(x, ';');
    feature<-lapply(x, function(x) x[-1]);
    names(feature)<-sapply(x, function(x) x[1]);
    feature;
  });
  names(cs)<-names(grp[['CS']]);
  
  # contributor
  cn<-lapply(grp[['CN']], as.vector);
  
  # created date
  cd<-lapply(split(grp[['CD']], names(grp[['CD']])), function(x) as.vector(unlist(x, use.names=FALSE))); 
            
  # edit history
  ed<-lapply(split(grp[['ED']], names(grp[['ED']])), function(x) as.vector(unlist(x, use.names=FALSE))); 
  
  #ALLELIC VARIANTS 
  av<-sapply(grp[['AV']], function(x) paste(x, collapse=' '));
  
  # reference
  rf<-lapply(grp[['RF']], function(x) {
    x[grep('^[0-9]+. ', x)]<-paste('\n', x[grep('^[0-9]+. ', x)], sep='');
    x<-paste(x, collapse=' ');
    x<-sub('^\n', '', x);
    strsplit(x, '\n');
  });
  
  # see also
  sa<-lapply(grp[['SA']], function(x) {
    x<-paste(x, collapse=' ');
    strsplit(x, '; ')[[1]];
  });
  
  ################################################################################################
  # Load OMIM to gene mapping
  omim2gene<-read.table(paste(path, 'src', 'mim2gene.txt', sep='/'), sep='\t', comment.char='', row.names=1, header=TRUE, stringsAsFactors=FALSE);
  if (!setequal(rownames(omim2gene), id)) warnings("Unmatched OMIM IDs: OMIM2gene table vs. Full OMIM records\n");
  omim2gene[omim2gene=='-']<-'';
  saveRDS(omim2gene, file=paste(path, 'r', 'omim2gene.rds', sep='/'));
    
  ################################################################################################
  # prepare OMIM full annotation table
  omim<-data.frame(Title=title[rownames(omim2gene)], row.names=rownames(omim2gene), stringsAsFactors=FALSE);
  omim$Type<-omim2gene[, 'Type']
  omim$Gene_ID<-omim2gene[, "Entrez.Gene.ID"];
  omim$Gene_Symbol<-omim2gene[, "Approved.Gene.Symbol"];
  omim$Status<-status[rownames(omim)];
  omim[is.na(omim)]<-""; 
  saveRDS(omim, file=paste(path, 'r/omim.rds', sep='/'));
  
  ################################################################################################
  # Prepare a full information list by OMIM IDs
  cat('Structure entries by ID');
  by.id<-lapply(rownames(omim), function(id) {
    lst<-list(
      ID = id,
      Title = title[id],
      Alternative_Title = title.alt[[id]],
      Type = omim[id, 'Type'],
      Status = omim[id, 'Status'],
      Gene_ID =  omim[id, 'Gene_ID'],
      Gene_Symbol = omim[id, 'Gene_Symbol'],
      Gene_Ensembl = omim2gene[id, 'Ensembl.Gene.ID'],
      Text = if (id %in% names(tx)) tx[id][[1]] else c(),
      Clinical_Synopsis = if (id %in% names(cs)) cs[[id]] else c(),
      Allelic_Variants = if (id %in% names(av)) av[[id]] else c(),
      References = if (id %in% names(rf)) rf[[id]] else c(),
      See_Also = if (id %in% names(sa)) sa[[id]] else c(),
      Contributors = if (id %in% names(cn)) cn[[id]] else c(),
      Created_Date = if (id %in% names(cd)) cd[[id]] else c(),
      Edit_History = if (id %in% names(ed)) ed[[id]] else c()
    );
    lst[sapply(lst, length)>0];
  }); 
  saveRDS(by.id, file=paste(path, 'r/omim_by_id.rds', sep='/'));
  
  omim;
}