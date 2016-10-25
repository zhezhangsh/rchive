library(VariantAnnotation);
library(BSgenome.Hsapiens.UCSC.hg19);


# VCF file with snpEff annotation
fn.vcf<-"/nas/is1/rchive/data/variant/public/g1k/src/ALL.chr1_to_chrX_phase3_shapeit2_mvncall_integrated.20130502.genotypes.fixed.vcf.bgz";
# Output file with snpEff annotation
fn.out<-"/nas/is1/rchive/data/variant/public/g1k/src/snpeff.txt";

# Run parser
source('~/R/source/GtUtility/R/getEff.r');
parseEFF(fn.vcf, fn.out, -1, 1000000, 10);

eff<-read.table(fn.out, he=FALSE, row=1, skip=601)


