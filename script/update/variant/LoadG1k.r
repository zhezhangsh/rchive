library(VariantAnnotation);
library(BSgenome.Hsapiens.UCSC.hg19);

library(devtools);
install_github("zhezhangsh/GtUtility");
library(GtUtility);

path<-paste(Sys.getenv("RCHIVE_HOME"), 'data/variant/public/g1k', sep='/');

if (!file.exists(path)) dir.create(path, recursive=TRUE);
if(!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'), recursive=TRUE);
if(!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'), recursive=TRUE);


# This is a VCF file based on the ones downloaded from 1000 Genomes FTP site
# --- All chromosomes are combined
# --- The problem with Ns in ALT was fixed
fn.vcf<-"/nas/is1/reference/human/1kg_hg19_GTcalls/ALL.chr1_to_chrX_phase3_shapeit2_mvncall_integrated.20130502.genotypes.fixed.vcf.gz";
if (!file.exists(fn.vcf)) stop("Error: VCF file not found");

# split genome into 10 Mbp sections for uploading
chr<-c(1:22, 'X', 'Y');
len<-seqlengths(Hsapiens)[paste('chr', chr, sep='')];
names(len)<-chr;
rng<-lapply(chr, function(chr) data.frame(chr=rep(chr, ceiling(len[chr]/10^7)), start=seq(1, len[chr], 10000000), end=c(seq(10000000, len[chr], 10000000), len[chr]), stringsAsFactors=FALSE));
rng<-do.call('rbind', rng);
gr<-GRanges(rng[,1], IRanges(rng[,2], rng[,3]));

# File names of loaded VCF 
fn<-paste(rng[, 1], ceiling(rng[,2]/10^7), sep='_');
names(gr)<-fn;
gr$file<-paste(path, '/r/', fn, '.rds', sep='');

sapply(1:length(gr), function(i) {
  g<-gr[i];
  if (!file.exists(g$file)) { # load from one 10-Mbp section 
    cat('\nLoading region: ', names(g));
    vcf<-readVcf(fn.vcf, 'hg19', param=ScanVcfParam(which=g));
    cat('\nLoaded', nrow(info(vcf)), 'variants.');
    saveRDS(vcf, file=g$file);
  } else vcf<-readRDS(g$file); # previously loaded
  
  # Parse VCF
  cat('\nParsing regions: ', names(g));
  gt<-parseVcf(vcf);
  saveRDS(gt, file=paste('genotype', g$file, sep='_'));
});




