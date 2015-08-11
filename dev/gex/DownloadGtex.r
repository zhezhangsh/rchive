DownloadGtex<-function(path=paste(Sys.getenv("RCHIVE_HOME"), 'data/gex/public/gtex', sep='/'),
                       url="http://www.gtexportal.org/static/datasets/gtex_analysis_v4",
                       url.fn=c(fpkm = "rna_seq_data/GTEx_Analysis_V4_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz",
                                smp.dd = "annotations/GTEx_Data_V4_Annotations_SampleAttributesDD.xlsx",
                                smp.ds = "annotations/GTEx_Data_V4_Annotations_SampleAttributesDS.txt", 
                                sub.dd = "annotations/GTEx_Data_V4_Annotations_SubjectPhenotypes_DD.xlsx",
                                sub.ds = "annotations/GTEx_Data_V4_Annotations_SubjectPhenotypes_DS.txt"),
                       fn.anno=paste(Sys.getenv("RCHIVE_HOME"), "data/gene/public/entrez/r/human_genes_full.rds", sep='/'),
                       download.all = FALSE) {  

  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  if (!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'));
  if (!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'));
  
  nm<-names(url.fn);
  url.fn<-paste(url, url.fn, sep='/');
  fn<-paste(path, 'src', sapply(strsplit(url.fn, '/'), function(f) f[length(f)]), sep='/');
  names(fn)<-url.fn;
  if (!download.all) {
    fn<-fn[!file.exists(fn)];
    if (length(fn) > 0) f<-sapply(names(fn), function(url) download.file(url, fn[url]));
  }
  fn<-paste(path, 'src', sapply(strsplit(url.fn, '/'), function(f) f[length(f)]), sep='/');
  names(fn)<-nm;

  # Load FPKM data
  fn.rpkm<-paste(path, 'r', 'rpkm.rds', sep='/');
  if (download.all | !file.exists(fn.rpkm)) {
    gex<-read.table(fn[1], sep='\t', skip=2, header=TRUE, row=1, stringsAsFactors=FALSE);
    colnames(gex)<-gsub('\\.', '-', colnames(gex));
    anno.full<-readRDS(fn.anno);
    gex<-gex[gex[, 1] %in% anno.full$Symbol, , drop=FALSE];
    gex<-gex[rowSums(gex[, -1]) > 0, , drop=FALSE];
    sym2id<-lapply(split(rownames(gex), gex[, 1]), unique);
    id<-sapply(sym2id, function(id) if (length(id)==1) id else {sm<-rowSums(gex[id, -1]); id[sm==max(sm)][1]});
    gex<-gex[id, , drop=FALSE];
    rownames(gex)<-gex[, 1];
    gex<-as.matrix(gex[, -1]);
    id<-rownames(anno.full);
    names(id)<-as.vector(anno.full$Symbol);
    id<-id[rownames(gex)];
    rownames(gex)<-id;
    anno<-anno.full[rownames(anno.full) %in% id, , drop=FALSE];
    gex<-gex[rownames(anno), , drop=FALSE];
    saveRDS(gex, fn.rpkm);
  } else gex<-readRDS(fn.rpkm);
  
  # Process annotation info.
  meta<-lapply(names(fn[-1]), function(nm) {
    f<-fn[nm];
    if (grepl('.xlsx$', f)) {
      m<-xlsx::read.xlsx(f, 1, header=TRUE);
      m<-m[, !grepl('^NA.', names(m))];
      rownames(m)<-m[, 1];
      m<-m[, -1];
    } else {
      m<-read.csv2(f, sep='\t', header=TRUE, row=1);
    }
    m;
  });
  names(meta)<-names(fn)[-1];
  saveRDS(meta, file=paste(path, 'r', 'original_metadata.rds', sep='/'));
  
  # Define groups by tissue type
  smp0<-meta$smp.ds;
  smp0<-smp0[rownames(smp0) %in% colnames(gex), , drop=FALSE];
  sub0<-meta$sub.ds;
  tis<-smp0$SMTS;
  tis<-gsub(' ', '_', tis);
  tis<-gsub('_Tissue$', '', tis);
  tis.spe<-smp0$SMTSD;
  tis.spe<-sapply(strsplit(tis.spe, ' - '), function(x) x[length(x)]);
  tis.spe<-sapply(strsplit(tis.spe, '\\('), function(x) x[1]);
  tis.spe<-gsub(' ', '_', tis.spe);
  tis.spe<-gsub('_$', '', tis.spe);
  
  smp.nm<-rownames(smp0);
  smp.id<-sapply(strsplit(smp.nm, '-'), function(x) x[length(x)])
  sub.id<-paste('GTEX', sapply(strsplit(smp.nm, '-'), function(x) x[2]), sep='-');
  age<-sub0[sub.id, 'AGE'];
  age<-paste(sapply(strsplit(age, '-'), function(x) x[1]), 's', sep='')
  dth<-sub0[sub.id, 'DTHHRDY'];
  gnd<-sub0[sub.id, 'GENDER'];
  gnd<-c('Male', 'Female')[gnd];
  
  smp<-data.frame(row.names=smp.nm, stringsAsFactors=FALSE,
                  Name=smp.nm, GTEX_ID=smp.nm, Group='', Dataset='', Donor=sub.id, Gender=gnd, Age=age, Death_Hardy_Code=dth,
                  Tissue=tis, Tissue_Specific=tis.spe, Site=smp0$SMCENTER, RIN=smp0$SMRIN, Ischemic_Time=smp0$SMTSISCH,
                  Isolation_Batch=smp0$SMNABTCH, Isolation_Type=smp0$SMNABTCHT, Isolation_Date=smp0$SMNABTCHD,
                  Experiment_Batch=smp0$SMGEBTCH, Experiment_Type=smp0$SMGEBTCHT, Experiment_Date=smp0$SMGEBTCHD);

  # data set table
  ds.nm<-paste(tis, tis.spe, sep='-');
  ds.nm[tis==tis.spe]<-tis[tis==tis.spe];
  ds<-data.frame(Name=unique(ds.nm), Num_Gene=nrow(gex), Num_Group=0, Num_Sample=0, Species='Human', stringsAsFactors=FALSE);
  ds<-ds[order(ds[,1]), , drop=FALSE];
  rownames(ds)<-ds.id<-sub('^1', 'D', 10000+1:nrow(ds));  
  ds$Num_Sample<-sapply(ds$Name, function(nm) length(ds.nm[ds.nm==nm]));
  
  names(ds.id)<-ds$Name;
  smp$Dataset<-ds.id[ds.nm];
  
  # group table
  grp.nm<-paste(age, ds.nm, sep='-');
  grp2smp<-split(rownames(smp), grp.nm);
  grp2ds<-sapply(split(ds.nm, grp.nm), unique);
  grp<-data.frame(Name=names(grp2smp), Dataset=ds.id[grp2ds], Num_Sample=sapply(grp2smp, length), Age=sapply(strsplit(names(grp2smp), '-'), function(x) x[1]), stringsAsFactors=FALSE);
  grp<-grp[order(grp$Dataset, grp$Name), , drop=FALSE];
  grp$Tissue<-ds[grp$Dataset, 'Name'];
  rownames(grp)<-grp.id<-sub('^1', 'G', 10000+1:nrow(grp));
  ds$Num_Group<-sapply(split(rownames(grp), grp$Dataset)[rownames(ds)], length);
  
  names(grp.id)<-grp$Name;
  smp$Group<-grp.id[grp.nm];
  
  smp$Name<-paste(sub('^GTEX-', '', smp$Donor), ds[smp$Dataset, 'Name'], sep='_');
  smp<-smp[order(smp$Dataset, smp$Group, smp$Name), ];
  rownames(smp)<-sub('^1', 'S', 10000+1:nrow(smp));
  
  meta<-list(Dataset=ds, Group=grp, Sample=smp);
  saveRDS(meta, file=paste(path, 'r', 'metadata.rds', sep='/'));
  
  # Finalize FPKM data
  gex<-gex[, smp$GTEX_ID, drop=FALSE];
  colnames(gex)<-rownames(smp);
  gex<-log2(gex+1);
  
  gex<-lapply(rownames(ds), function(id) gex[, rownames(smp)[smp$Dataset==id], drop=FALSE]);
  names(gex)<-rownames(ds);
  gex<-lapply(gex, function(gex) affy::normalize.loess(gex, log.it=FALSE)); 
  
  saveRDS(gex, file=paste(path, 'r', 'gex.rds', sep='/'));
  
  bro<-lapply(meta, function(meta) data.frame(ID=rownames(meta), meta, stringsAsFactors=FALSE));
  names(bro)[1]<-'Data set';
  bro$Gene<-data.frame(row.names=rownames(anno), stringsAsFactors=FALSE,
                       ID=awsomics::AddHref(rownames(anno), awsomics::UrlEntrezGene(rownames(anno))), 
                       Species='human', Name=as.vector(anno$Symbol), Num_Dataset=nrow(ds), Num_Sample=nrow(smp), 
                       Type=anno$type_of_gene, Title=anno$description);
  saveRDS(bro, file=paste(path, 'r', 'browse_table.rds', sep='/'));
  
  meta;
}