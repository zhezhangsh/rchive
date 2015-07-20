# Download and process Mouse Atlas of Gene Expression data set from mouse development SAGE data

DownloadMouseAtlas<-function(url.table="http://www.mouseatlas.org/data/mouse/project_libraries_all_view", 
                             path=paste(Sys.getenv("RCHIVE_HOME"), 'data/gex/public/mage', sep='/'), 
                             fn.anno=paste(Sys.getenv("RCHIVE_HOME"), "data/gene/public/entrez/r/mouse_genes_full.rds", sep='/'),
                             update.all=FALSE, tag.length=17) {
  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  if (!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'));
  if (!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'));
  
  library(RCurl);
  library(XML);
  
  # get sample IDs
  ln<-strsplit(getURL(url.table), '\n')[[1]];
  ids<-sort(sapply(strsplit(CleanHtmlTags(ln[grep('/libraries/SM', ln)]), ' '), function(l) l[1]));

  # get sample metadata
  url.smp<-paste("http://www.mouseatlas.org/data/mouse/libraries", ids, sep='/');
  fn<-paste(path, 'r', 'original_meta.rds', sep='/');
  if (!file.exists(fn) | update.all) {
    meta<-lapply(url.smp, function(u) {
      cat(u, '\n');
      ln<-strsplit(getURL(u), '\n')[[1]];
      ind<-grep('<th>', ln);
      k<-CleanHtmlTags(ln[ind], FALSE);
      v<-CleanHtmlTags(ln[ind+1], FALSE);
      names(v)<-k;
      v;
    }); 
    names(meta)<-ids;
    saveRDS(meta, file=fn);
  } else meta<-readRDS(fn);

  # Download individual libraries
  fn<-sapply(ids, function(id) {
    fn<-paste(path, '/src/', id, '.txt', sep='');
    url<-paste("http://www.mouseatlas.org/data/mouse/libraries/", id, "/plain_text/tagcounts_", id, "_Q99.txt", sep='');
    if (!file.exists(fn) | update.all) download.file(url, fn);
    fn;
  });

  # Read in tag count
  ct<-lapply(fn, function(fn) { 
    c<-read.table(fn, sep='\t'); 
    ct<-c[,2]; 
    names(ct)<-as.vector(c[, 1]); 
    ct; 
  })
  # combine tags
  seq<-sort(unique(unlist(lapply(ct, names)), use.names=FALSE));
  seq<-seq[nchar(seq)==tag.length];

  # Align all tags to chromosomes one by one
  library(BSgenome.Mmusculus.UCSC.mm10);
  chr<-names(Mmusculus)[1:22];
  fn<-sapply(chr, function(chr) {
    cat('Mapping to', chr, '\n');
    fn<-paste(path, '/r/tag2genome_', chr, '.rds', sep='');
    if (!file.exists(fn) | update.all) {
      gr<-MatchSage(seq, Mmusculus, chr);
      saveRDS(gr, file=fn);
    }
     fn;
  });
  
  # Alignment results
  fn<-paste(path, '/r/tag2genome', '.rds', sep='');
  if (file.exists(fn)) gr<-readRDS(fn) else {
    grs<-lapply(fn, readRDS);
    gr<-grs[[1]];
    for (i in 2:length(grs)) gr<-c(gr, grs[[i]]); 
    saveRDS(gr, fn);
  }
  
  # Map to exon/gene
  library(TxDb.Mmusculus.UCSC.mm10.knownGene);
  txdb<-TxDb.Mmusculus.UCSC.mm10.knownGene;
  ex<-exons(txdb, columns=c('exon_id', 'gene_id')); # all known exons
  gr.hit<-gr[countOverlaps(gr, ex, type='within')>0]; # tag mapped to at least one exon
  olap<-as.matrix(findOverlaps(gr.hit, ex, type='within'));
  gr.olap<-gr.hit$Seq_ID[olap[,1]];
  gn.olap<-ex$gene_id[olap[,2]];
  gr2gn<-split(unlist(gn.olap), rep(gr.olap, elementLengths(gn.olap))); # tag to gene
  gr2gn<-lapply(gr2gn, unique);
  gr2gn<-lapply(gr2gn, function(gg) gg[!is.na(gg)]);
  gr2gn.uni<-unlist(gr2gn[sapply(gr2gn, length)==1]); # tags mapped to one and only one gene
  gr.mapped<-gr[gr$Seq_ID %in% as.integer(names(gr2gn.uni))];
  occur<-as.vector(table(gr.mapped$Seq_ID)[names(gr2gn.uni)]);
  tag2gn<-data.frame(row.names=as.character(seq[as.integer(names(gr2gn.uni))]), Gene_ID=as.vector(gr2gn.uni), Hits_Total=occur, stringsAsFactors=FALSE); # Tag to gene final results  
  
  # Tag count matrix
  cnt<-matrix(0, nr=nrow(tag2gn), nc=length(ct), dimnames=list(rownames(tag2gn), names(ct)));
  for (i in 1:length(ct)) {
    c<-ct[[i]];
    c<-c[names(c) %in% rownames(tag2gn)];
    if (length(c)>0) cnt[names(c), i]<-as.integer(c);
  }
  saveRDS(tag2gn, file=paste(path, 'r', 'tag2gene_all.rdse', sep='/'));
  saveRDS(cnt, file=paste(path, 'r', 'tagcount_all.rds', sep='/'));
  
  # Gene Count matix
  cnt<-cnt[tag2gn[,2]<=10, colSums(cnt)>10000];
  gn<-tag2gn[rownames(cnt), 1];
  c<-apply(cnt, 2, function(cnt) sapply(split(cnt, gn), sum));
  c<-c[, colSums(c)>10000];
  saveRDS(c, file=paste(path, 'r', 'gene_count.rds', sep='/'));
  
  fn.expr<-paste(path, 'r', 'expr.rds', sep='/');
  if (!file.exists(fn.expr) | update.all==TRUE) {
    library(affy);
    expr<-log2(c+1);
    expr<-normalize.loess(expr, which(rowMeans(expr)>median(rowMeans(expr))), log.it=FALSE);
    saveRDS(expr, file=paste(path, 'r', 'expr.rds', sep='/'));
  } else expr<-readRDS(fn.expr);
 
  
  # Clean up sample metadata
  meta<-meta[colnames(expr)];
  cnm<-Reduce('intersect', lapply(meta, names));
  smp<-t(sapply(meta, function(meta) meta[cnm]));
  saveRDS(smp, file=paste(path, 'r', 'meta_full.rds', sep='/'));
  l<-strsplit(smp[, 'Location'], ' [-–] ');
  loc<-sapply(l, function(l) l[1]);
  tis<-sapply(l, function(l) l[2]);
  smp<-data.frame(row.names=smp[, "Library Name"], Organ=loc, Part=tis, Age=smp[, 'Age'], Stage=smp[, "Developmental Stage"], Sex=smp[, 'Sex'],
                  Strain=smp[, 'Strain'], Condition=smp[, "Organsim Condition"], Protocol=smp[, "SAGE Protocol"], 
                  Lab=smp[, 'Lab'], stringsAsFactors=FALSE);  
  saveRDS(smp, file=paste(path, 'r', 'sample.rds', sep='/'));
  
  ##########################################################################################################################
  ### Prepare files required by GeEx as a data collection
  
  ## Final data matrix
  anno=readRDS(fn.anno);
  gex<-expr[rownames(anno[rownames(anno) %in% rownames(expr), ]), ];
  anno<-anno[rownames(gex), ];
  
  # re-order samples
  smp.e<-smp[grep('^E', smp$Age), ];
  smp.e<-smp.e[order(as.numeric(sub('^E', '', smp.e$Age))), ];
  smp.p<-smp[-grep('^E', smp$Age), ];
  smp.p<-smp.p[order(smp.p$Age), ];
  smp<-rbind(smp.e, smp.p);
  smp<-smp[order(smp$Organ), ];
  smp[is.na(smp)]<-'';
  
  # Prepare metadata IDs
  ds.name<-unique(smp$Organ);
  ds.id<-sub('1', 'D', 10000+1:length(ds.name));
  names(ds.id)<-ds.name;
  grp.name<-paste(smp$Organ, smp$Age, sep=':');
  names(grp.name)<-smp$Organ;
  grp.name<-grp.name[!duplicated(grp.name)];
  grp.id<-sub('1', 'G', 10000+1:length(grp.name));
  names(grp.id)<-grp.name;
  smp.id<-sub('1', 'S', 10000+1:nrow(smp));
  smp.name<-paste(smp$Organ, smp$Age, sep=':');
  names(smp.name)<-smp.id;
  grp2smp<-split(smp.name, smp.name);
  smp.name<-paste(unlist(grp2smp), unlist(lapply(grp2smp, function(x) 1:length(x))), sep='_');
  names(smp.name)<-unlist(lapply(grp2smp, names));
  smp.name<-smp.name[smp.id];
  
  # metadata tables
  tbl.smp<-data.frame(row.names=smp.id, stringsAsFactors=FALSE, Name=smp.name, Original_ID=rownames(smp), 
                      Group=grp.id[paste(smp$Organ, smp$Age, sep=':')], Dataset=ds.id[smp$Organ], smp);
  tbl.grp<-data.frame(row.names=grp.id, stringsAsFactors=FALSE, Name=grp.name, Dataset=ds.id[names(grp.name)], Num_Sample=sapply(grp2smp, length)[grp.name], 
                      Organ=sapply(strsplit(grp.name, ':'), function(x) x[1]), Age=sapply(strsplit(grp.name, ':'), function(x) x[2]));
  tbl.ds<-data.frame(row.names=ds.id, stringsAsFactors=FALSE, Name=ds.name, Num_Gene=nrow(gex), 
                     Num_Group=as.integer(table(tbl.grp$Dataset)[ds.id]), Num_Sample=as.integer(table(tbl.smp$Dataset)[ds.id]), Species='Mouse');
  metadata<-list(Dataset=tbl.ds, Group=tbl.grp, Sample=tbl.smp);
  
  # Browse tables
  brow.gene<-data.frame(row.names=rownames(anno), stringsAsFactors=FALSE, ID=awsomics::AddHref(rownames(anno), paste('http://www.ncbi.nlm.nih.gov/gene/?term=', rownames(anno), sep='')),
                        Species='mouse', Name=as.vector(anno$Symbol), Num_Dataset=nrow(tbl.ds), Num_Sample=nrow(tbl.smp), Type=as.vector(anno$type_of_gene), Title=anno$description);
  brow.smp<-data.frame(ID=rownames(tbl.smp), tbl.smp, stringsAsFactors=FALSE);
  brow.smp$Original_ID<-awsomics::AddHref(brow.smp$Original_ID, paste('http://www.mouseatlas.org/data/mouse/libraries/', brow.smp$Original_ID, sep=''));
  brow.grp<-data.frame(ID=rownames(tbl.grp), tbl.grp, stringsAsFactors=FALSE);
  brow.ds<-data.frame(ID=rownames(tbl.ds), tbl.ds, stringsAsFactors=FALSE);
  browse_table<-list(Dataset=brow.ds, Group=brow.grp, Sample=brow.smp, Gene=brow.gene);
  
  saveRDS(gex, file=paste(path, 'r', 'gex.rds', sep='/'));
  saveRDS(metadata, file=paste(path, 'r', 'metadata.rds', sep='/'));
  saveRDS(browse_table, file=paste(path, 'r', 'browse_table.rds', sep='/'));
  
  metadata;
  ##########################################################################################################################  
}

# Match SAGE tags to genome
MatchSage<-function(seq, genome, chromosome=names(genome)) {
  # seq         Sequences of SAGE tags, all should have the same length
  # genome      Reference genome as an BSgenome object
  # chromosome  Name of the chromosomes to map to
  seq<-DNAStringSet(seq);
  chromosome<-chromosome[chromosome %in% names(genome)];
  if (length(chromosome) == 0) NA else {
    gr<-lapply(chromosome, function(chr) {
      cat('Aligning to', chr, '\n');
      mt1<-matchPDict(PDict(seq), genome[[chr]]);
      cb1<-do.call('c', mt1);
      gr1<-GRanges(rep(chr, length(cb1)), IRanges(start(cb1), end(cb1)), strand='+');
      gr1$Seq_ID<-rep(1:length(seq), sapply(mt1, length));
      
      mt2<-matchPDict(PDict(reverseComplement(seq)), genome[[chr]]);
      cb2<-do.call('c', mt2);
      gr2<-GRanges(rep(chr, length(cb2)), IRanges(start(cb2), end(cb2)), strand='-');
      gr2$Seq_ID<-rep(1:length(seq), sapply(mt2, length));
      
      c(gr1, gr2);
    });
    do.call('c', gr);
  }
}