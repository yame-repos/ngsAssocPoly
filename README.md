# ngsAssocPoly
This repository contains phenotype data, genotype data and an R package ("ngsAssocPoly") corresponding to a manuscript titled "Genetic Mapping in Autohexaploid Sweet Potato with Low-coverage NGS-based Genotyping Data".
The main objective of the study was to performe genetic mapping in polyploids using low-coverage NGS data.

The following files are included with this distribution.

  - For Windows:   		ngsAssocPoly_1.0.0.zip
  - For MacOS and Linux:	ngsAssocPoly_1.0.0.tgz
  - Reference manual:	ngsAssocPoly_manual.pdf
  - ddRAD-seq data:	KX-F1.vcf.gz, X18-S1.vcf.gz
  - Phenotype data: KX-F1_Pheno.csv, X18-S1_Pheno.csv
  - Genetic map information: GeneticMap_Shirasawa.et.al_2016.csv (https://www.nature.com/articles/srep44207)

<!-- end list -->

## Example

Following libraries are required to perform this analysis.
``` r
require(qqman)
require(vcfR)
require(RColorBrewer)
require(ngsAssocPoly)
```

Simple (but low accuracy) allele dosage calculation is performed like
``` r 
alleleDosageEstimation("KX-F1.vcf.gz")
```

Use other libraries if you would like to perform allele dosage estimation wiht higher accuracy.
Following script is an example that uses "updog" (Gerard et al. 2018, https://www.genetics.org/content/210/3/789)
``` r 
vcf.file.name = "KX-F1.vcf.gz"
project.id = strsplit(vcf.file.name, ".vcf")[[1]][1]
Map.name = paste(project.id, "_Map.rda", sep = "")
Geno.name = paste(project.id, "_Geno.rda", sep = "")

count.miss = function(vec) length(vec[is.na(vec)])

vcf = vcfR::read.vcfR(vcf.file.name)
chrom = vcfR::getCHROM(vcf)
pos = vcfR::getPOS(vcf)
marker = paste(chrom, pos, sep = "_")
Map = data.frame(marker = marker, chrom = chrom, pos = pos)
rm(chrom, pos, marker) ; gc()

dp = vcfR::extract.gt(vcf, element="DP")
mode(dp) = "numeric"
rd = vcfR::extract.gt(vcf, element="RD")
mode(rd) = "numeric"
rm(vcf) ; gc()

dp = dp[ , !is.element(colnames(dp), rm.col)]
rd = rd[ , !is.element(colnames(rd), rm.col)]

dp[dp > max.dp] = NA
dp[dp < min.dp] = NA

miss = apply(dp, 1, count.miss) / ncol(dp)
use = miss < max.miss
Map = Map[use, , drop = FALSE]
dp = dp[use, , drop = FALSE]
rd = rd[use, , drop = FALSE]
rm(max.dp, min.dp, count.miss, miss)
gc()

mout = updog::multidog(refmat = rd, sizemat = dp, 
                       ploidy = ploidy, model = model, 
                       p1_id = p1_id, p2_id = p2_id, 
                       nc = n.cores)
indmat = as.matrix(mout$inddf[,-c(1:7)])
rownames(indmat) = mout$inddf$ind
snpix = mout$inddf$snp

markers = mout$snpdf$snp
n.markers = length(markers)
sample.names = unique(mout$inddf$ind)
n.samples = length(sample.names)

Geno = vector(mode="list", length=n.markers)
names(Geno) = markers
for (p in 1:n.markers) {
  print(paste(p, "/", n.markers))
  Geno[[p]] = indmat[is.element(snpix, markers[p]), ]
}#p

freq = numeric(n.markers)
not.miss = numeric(n.markers)
for (i in 1:n.markers) {
  m = Geno[[i]]
  m = m[!is.na(m[,1]),,drop=FALSE]
  not.miss[i] = nrow(m) / nrow(Geno[[i]])
  freq[i] = max(apply(m, 2, sum) / sum(m))
}
good = (not.miss > (1-max.miss)) & (freq < max.freq)

Geno = Geno[good]
Map = Map[good,]

save(Geno, file=Geno.name)
save(Map, file=Map.name)
```

Following two files will be created.

  - KX-F1_Geno.rda (Allele dosage matrices)
  - KX-F1_Map.rda (Map position of SNP markers in KX-F1_Geno.rda)

<!-- end list -->

The reference genome used in the tutorial dataset is not chromosome-level genome.
Therefore, we allocated the SNP markers to linkage groups identified in Shirasawa et al. (2017, https://www.nature.com/articles/srep44207)
``` r
load("KX-F1_Map.rda")
genetic.map = read.csv("GeneticMap_Shirasawa.et.al_2016.csv")
LG = as.character(unique(genetic.map$tempLG))
Scaffold = as.character(unique(Map$chrom))
M = matrix(0, nrow=length(Scaffold), ncol=length(LG))
rownames(M) = Scaffold
colnames(M) = LG
for (i in 1:nrow(M)) {
  ix = rownames(M)[i]
  ix = unique(genetic.map$tempLG[genetic.map$Scaffold==ix])
  if (length(ix)!=1) next
  M[i,is.element(colnames(M), ix)] = 1
}
for (i in 1:ncol(M)) {
  ix = rownames(M)[M[,i]==1]
  Map$tempLG[is.element(Map$chrom, ix)] = colnames(M)[i]
}
Map$tempLG[is.na(Map$tempLG)] = "Unanc"
LG = c(sort(unique(Map$tempLG)), "Unanc")
for (i in LG) {
  ix = nrow(Map[Map$tempLG==i,,drop=FALSE])
  Map$order[Map$tempLG==i] = 1:ix
}
Map$chrom = Map$tempLG
Map$pos = Map$order
head(Map)
Map = Map[,is.element(colnames(Map), c("marker","chrom","pos"))]
save(Map, file = "KX-F1_Map.rda")
```

The association test is performed like
``` r
alleleDosageGLM(Geno.name = "KX-F1_Geno.rda",
                Map.name = "KX-F1_Map.rda",
                pheno.file.name = "KX-F1_Pheno.csv")
```

The products are

  - KX-F1_scores.csv (Scores of the association test)
  - KX-F1_color__manhattan.jpeg, KX-F1_internode_length__manhattan.jpeg (Manhattan plots)
  - KX-F1_color__qq.jpeg, KX-F1_internode_length__qq.jpeg (QQ plots)

<!-- end list -->
