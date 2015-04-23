# Map and update human chromosome names
# The mapping uses the first-come-first-use rule to avoid ambiguity, the final data might miss some rare contig
path<-paste(RCHIVE_HOME, 'data/assembly/public/chromosome', sep='/');

if (!file.exists(path)) dir.create(path, recursive=TRUE);
if (!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'), recursive=TRUE);
if (!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'), recursive=TRUE);

# file names to chromosome mapping and sets
fn.map<-paste(path, 'r', 'human_chromosome_mapping.rds', sep='/');
fn.set<-paste(path, 'r', 'human_chromosome_sets.rds', sep='/');

if (file.exists(fn.map)) mp<-readRDS(fn.map) else mp<-list();
if (file.exists(fn.set)) st<-readRDS(fn.set) else st<-list();

# NCBI versions; name a version if want to create a set for it
# version 18 to 20 doesn't exist
url<-paste("ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.", c(29:21, 17, 14:8), ".assembly.txt", sep='');
names(url)<-c('GRCh38.p3', 'GRCh38.p2', 'GRCh38.p1', '', 'GRCh37.p13', 'GRCh37.p12', 'GRCh37.p11', 'GRCh37.p10');
names(url)[is.na(names(url))]<-'';

# Loading all tables
tbls<-lapply(url, function(url) {
  cat(url, '\n'); 
  read.table(url, sep='\t', stringsAsFactors=FALSE);
})
tbls<-lapply(tbls, as.matrix);

##############################################################################################
# merging new mapping list to existing one
mergeList<-function(mp0, mp1) {
  if (length(mp0) == 0) mp<-mp1 else { # if old list is empty, simply return the new list
    for (i in 1:length(mp1)) {
      nm1<-mp1[[i]];
      nm1<-nm1[nm1!='' & !is.na(nm1)];
      ind<-rep(1:length(mp), sapply(mp, length))[unlist(mp) %in% nm1];
      if (length(ind) == 0) mp0<-append(mp0, list(nm1)) else mp0[[ind[1]]]<-unique(c(mp0[[ind[1]]], setdiff(nm1, unlist(mp0[[-ind[1]]]))));
      #ct<-sapply(mp0, function(nm0) length(nm1[toupper(nm1) %in% toupper(nm0)]));
      #for (i in 1:length(ct)) if (ct[i]==0) mp0<-append(mp0, nm1[[i]]) else mp0[[i]]<-c(mp0[[i]], setdiff(nm1[[i]], unlist(mp0, use.names=FALSE)));
      #if (max(ct) == 0) mp0<-append(mp0, list(nm1)) else mp0[which(ct>0)]<-lapply(mp0[which(ct>0)], function(nm0) c(nm0, nm1));
    }
    mp<-mp0;
  }
  #mp<-lapply(mp, function(mp) c(mp, sub('_CTG[0-9]+', '', mp)));
  lapply(mp, function(mp) sort(unique(mp)));
}
##############################################################################################

##############################################################################################
# Merge new mapping lists to existing one
mp<-lapply(mp, function(x) x[x!='' & !is.na(x)]);

# UCSC extra chromosomes
mp1<-as.list(strsplit("chr6_apd_hap1;chr6_cox_hap2;chr6_dbb_hap3;chr6_mcf_hap5;chr6_qbl_hap6;chr4_ctg9_hap1;chr6_mann_hap4;chr6_ssto_hap7;chrUn_gl000211;chrUn_gl000212;chrUn_gl000213;chrUn_gl000215;chrUn_gl000218;chrUn_gl000219;chrUn_gl000220;chrUn_gl000222;chrUn_gl000223;chrUn_gl000227;chrUn_gl000228;chr17_ctg5_hap1;chr1_gl000191_random;chr1_gl000192_random;chr4_gl000193_random;chr4_gl000194_random;chr7_gl000195_random;chr17_gl000205_random;chr19_gl000209_random", ';')[[1]]);
mp<-mergeList(mp, mp1);

# Downloaded from NCBI
for (i in 1:length(tbls)) {
  print(i); 
  tb<-tbls[[i]]; 
  tb[tb=='na']<-'';
  cind<-c(1, 5, 7, 10);
  mp1<-lapply(1:nrow(tb), function(i) as.vector(tb[i, cind[cind<ncol(tb)]]));
  mp1<-lapply(mp1, function(x) x[x!='' & !is.na(x)]);
  mp<-mergeList(mp, mp1);
}

mp<-lapply(mp, sort);
saveRDS(mp, file=fn.map);
##############################################################################################

# save chromosome set of a specific version
ts<-tbls[names(tbls) != ''];
for (i in 1:length(ts)) {
  if (ncol(ts[[i]])>8) len<-as.numeric(ts[[i]][, 9]) else len<-rep(NA, nrow(ts[[i]]));
  names(len)<-ts[[i]][,1];
  st[as.character(names(ts)[i])][[1]]<-len;
}
saveRDS(st, file=fn.set);

# Mapping list indexed by chromosome name
mp.ind<-lapply(mp, function(mp) {
  ind<-lapply(1:length(mp), function(i) mp[-i]);
  names(ind)<-mp;
  ind;
});
mp.ind<-do.call('c', mp.ind);
mp.ind<-split(mp.ind, names(mp.ind));
mp.ind<-lapply(mp.ind, function(x) sort(unique(unlist(x, use.names=FALSE))));
saveRDS(mp.ind, file=sub('.rds', '_indexed.rds', fn.map));


##############################################################################################################
tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(RCHIVE_HOME, 'source/update/assembly/UpdateHumanChromosome.r', sep='/');
fn1<-paste(RCHIVE_HOME, '/source/update/assembly/log/', tm, '_UpdateHumanChromosome.r' , sep='');
file.copy(fn0, fn1)

