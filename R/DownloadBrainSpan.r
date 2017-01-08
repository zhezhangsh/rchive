DownloadBrainSpan<-function(url="http://www.brainspan.org/api/v2/well_known_file_download/267666525",
                            path=paste(Sys.getenv("RCHIVE_HOME"), 'data/gex/public/adhb', sep='/'),
                            fn.anno=paste(Sys.getenv("RCHIVE_HOME"), "data/gene/public/entrez/r/human_genes_full.rds", sep='/'),
                            update.all=FALSE) {
  
  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  if (!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'));
  if (!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'));
  
  # alternative way to group samples
  #path.str<-paste(path, '_structure');
  #if (!file.exists(path.str)) dir.create(path.str, recursive=TRUE);
  #if (!file.exists(paste(path.str, 'r', sep='/'))) dir.create(paste(path.str, 'r', sep='/'));
  #if (!file.exists(paste(path.str, 'src', sep='/'))) dir.create(paste(path.str, 'src', sep='/'));
  
  ####################################################################################################################################
  # Download source data
  fn.rnaseq<-paste(path, 'src', 'rnaseq_gene.zip', sep='/');
  fd.rnaseq<-sub('.zip$', '', fn.rnaseq);
  if (!file.exists(fn.rnaseq) | update.all) download.file(url, fn.rnaseq);
  unzip(fn.rnaseq, exdir=fd.rnaseq);
  
  # Exon microarray data
  #fn.array<-paste(path, 'src', 'array_gene.zip', sep='/');
  #fd.array<-sub('.zip$', '', fn.array);
  #if (!file.exists(fn.array) | update.all) download.file("http://www.brainspan.org/api/v2/well_known_file_download/267666527", fn.array);
  #unzip(fn.array, exdir=fd.array);

  meta<-read.table(paste(fd.rnaseq, 'columns_metadata.csv', sep='/'), row=1, he=TRUE, sep=',', stringsAsFactors=FALSE);
  saveRDS(meta, file=paste(path, 'r/original_metadata.rds', sep='/'));
  age<-meta$age;
  age.nm<-as.integer(sapply(strsplit(age, ' '), function(x) x[1]));
  age.ut<-sapply(strsplit(age, ' '), function(x) x[2]);
  age.ind<-which(age.ut=='yrs' & age.nm<=5);
  age.nm[age.ind]<-12*age.nm[age.ind];
  age.ut[age.ind]<-'mos';
  names(age.nm)<-names(age.ut)<-rownames(meta);

  # Match age to developomental stage
  stage<-readRDS(paste(path, 'src/developmental_stage.rds', sep='/')); # pre-processed file based on source document
  s<-sapply(names(age.nm), function(id) {
    ut<-stage[stage$Unit==age.ut[id], , drop=FALSE];
    rownames(ut)[ut$Min<=age.nm[id] & ut$Max>=age.nm[id]];
  });
  smp<-data.frame(stage=s, meta, stringsAsFactors=FALSE);

  ############################################################################################################
  fn.expr<-paste(path, 'r/expr.rds', sep='/');
  if (!file.exists(fn.expr) | update.all) {
    rpkm<-read.table(paste(fd.rnaseq, 'expression_matrix.csv', sep='/'), he=FALSE, row=1, sep=',');
    rpkm<-as.matrix(rpkm);
    anno<-read.table(paste(fd.rnaseq, 'rows_metadata.csv', sep='/'), he=TRUE, row=1, sep=',');
    gn<-readRDS(fn.anno);    
    anno[, 4]<-as.character(anno[, 4]);
    anno.no<-anno[!(anno[, 4] %in% rownames(gn)), , drop=FALSE];
    gn.id<-rownames(gn);
    names(gn.id)<-as.vector(gn$Symbol);
    anno.no[, 4]<-gn.id[as.vector(anno.no[, 3])];
    anno[rownames(anno.no), ]<-anno.no;
    anno<-anno[rownames(rpkm), ];
    anno<-anno[anno[, 4] %in% rownames(gn), ];
    anno<-anno[anno[, 4] %in% rownames(gn), ];
    
    expr<-log2(rpkm+1);
    expr<-expr[rownames(anno), ];
    rownames(expr)<-as.vector(anno[, 4]); 
    dup.expr<-expr[duplicated(rownames(expr)), ];
    dup.id<-unique(rownames(dup.expr));
    dup.expr<-t(sapply(dup.id, function(id) colMeans(expr[rownames(dup.expr) == id, ])));
    expr<-rbind(dup.expr, expr[!(rownames(expr) %in% rownames(dup.expr)), ]);
    
    anno<-gn[rownames(gn) %in% rownames(expr), ];
    expr<-expr[rownames(anno), ];
    colnames(expr)<-1:ncol(expr);
    
    expr<-affy::normalize.loess(expr, log.it=FALSE);

    saveRDS(anno, file=paste(path, 'r', 'annotation.rds', sep='/'));
    saveRDS(expr, file=fn.expr);
  } else {
    expr<-readRDS(fn.expr);
    anno<-readRDS(paste(path, 'r', 'annotation.rds', sep='/'));
  }
    
  ### Prepare metadata
  stage<-readRDS(paste(path, 'src/developmental_stage.rds', sep='/'));
  stage<-stage[rownames(stage) %in% s, , drop=FALSE];
  ds.name<-stage[[1]];
  ds.id<-sub('^1', 'D', 10000+1:nrow(stage));
  names(ds.id)<-rownames(stage);
  names(ds.name)<-ds.id;
  
  grp<-paste(sapply(strsplit(smp$stage, '_'), function(s) s[1]), smp$structure_acronym, sep='_');
  names(grp)<-ds.id[smp$stage];
  ds2grp<-lapply(split(as.vector(grp), names(grp)), function(g) sort(unique(g)));
  grp.name<-unlist(ds2grp, use.names=FALSE);
  grp.id<-sub('^1', 'G', 10000+1:length(grp.name));
  names(grp.id)<-grp.name;
  names(grp.name)<-grp.id;
  ds2grp<-lapply(ds2grp, function(g) {names(g)<-grp.id[g]; g});
  grp2ds<-rep(names(ds2grp), sapply(ds2grp, length));
  names(grp2ds)<-grp.id[unlist(ds2grp, use.names=FALSE)];
  
  smp.name<-paste(smp$structure_acronym, smp$donor_id, sep='_');
  smp$Name<-smp.name;
  grp2smp<-split(smp.name, grp.id[grp]);
  grp2smp<-lapply(grp2smp, sort);
  gid<-rep(names(grp2smp), sapply(grp2smp, length));
  smp.name<-unlist(grp2smp, use.names=FALSE);
  smp.id<-sub('^1', 'S', 10000+1:length(smp.name));
  names(smp.id)<-smp.name;
  grp2smp<-lapply(grp2smp, function(s) as.vector(smp.id[s]));
  smp2grp<-rep(names(grp2smp), sapply(grp2smp, length));
  names(smp2grp)<-unlist(grp2smp, use.names=FALSE);
  
  # Metadata tables
  tbl.smp<-data.frame(row.names=smp.id[smp$Name], stringsAsFactors=FALSE, Name=smp$Name, Original_ID=rownames(smp), Group=smp2grp[smp.id[smp$Name]], Dataset=grp2ds[smp2grp[smp.id[smp$Name]]],
                      Donor=smp$donor_name, Age=smp$age, Gender=smp$gender);
  tbl.smp<-cbind(tbl.smp, Developmental_Stage=ds.name[tbl.smp$Dataset], Part_Code=smp$structure_acronym, Part_Name=smp$structure_name);
  tbl.smp<-tbl.smp[order(rownames(tbl.smp)), ];
  
  tbl.grp<-data.frame(row.names=grp.id, stringsAsFactors=FALSE, Name=grp.name[grp.id], Dataset=grp2ds[grp.id], Num_Sample=sapply(grp2smp, length)[grp.id], Developmental_Stage=ds.name[grp2ds]);
  tbl.grp$Part_Code<-sapply(split(tbl.smp$Part_Code, tbl.smp$Group), unique)[rownames(tbl.grp)]
  tbl.grp$Part_Name<-sapply(split(tbl.smp$Part_Name, tbl.smp$Group), unique)[rownames(tbl.grp)]
  
  tbl.ds<-data.frame(row.names=ds.id, stringsAsFactors=FALSE, Name=ds.name, Num_Gene=nrow(expr), Num_Group=sapply(ds2grp, length)[ds.id], Num_Sample=as.vector(table(tbl.smp$Dataset)[ds.id]), Species='Human');
  tbl.ds$Stage_Code<-sapply(ds.id, function(id) names(ds.id[ds.id==id]));
  tbl.ds$Stage_Name<-sapply(ds.id, function(id) sub(paste('^', tbl.ds[id, 'Stage_Code'], '_', sep=''), '', tbl.ds[id, 'Name']));
  age<-stage[tbl.ds$Stage_Code, -1];
  colnames(age)<-c('Age_Min', 'Age_Max', 'Age_Unit');
  tbl.ds<-cbind(tbl.ds, age);
  tbl.ds$Available_Part<-sapply(split(tbl.smp$Part_Code, tbl.smp$Dataset), function(x) paste(sort(unique(x)), collapse='; '))[rownames(tbl.ds)]
  metadata<-list(Dataset=tbl.ds, Group=tbl.grp, Sample=tbl.smp);
  
  gex<-expr[, as.vector(tbl.smp$Original_ID)];
  colnames(gex)<-rownames(tbl.smp);
  gex<-lapply(rownames(tbl.ds), function(id) gex[, rownames(tbl.smp)[tbl.smp$Dataset==id, drop=FALSE]]);
  names(gex)<-rownames(tbl.ds);
  
  brow.tbl<-lapply(metadata, function(x) data.frame(ID=rownames(x), x));
  names(brow.tbl)<-c('Data set', 'Group', 'Sample');
  gene<-data.frame(ID=rownames(anno), Species='human', Name=anno$Symbol, Num_Dataset=length(gex), Num_Sample=nrow(tbl.smp), Type=anno$type_of_gene, Title=anno$description);
  gene$ID<-awsomics::AddHref(gene$ID, paste("http://www.ncbi.nlm.nih.gov/gene/?term=", gene$ID, sep=''));
  brow.tbl$Gene<-gene;
  
  saveRDS(metadata, file=paste(path, 'r', 'metadata.rds', sep='/'));
  saveRDS(brow.tbl, file=paste(path, 'r', 'browse_table.rds', sep='/'));
  saveRDS(gex, file=paste(path, 'r', 'gex.rds', sep='/'));
  
  metadata;
}