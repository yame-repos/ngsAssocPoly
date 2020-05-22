\name{alleleDosageEstimation}
\alias{alleleDosageEstimation}
\title{
Allele dosage estimation using NGS-based genotyping data
}
\description{
Performs allele dosage estimation from NGS-based genotyping data.
}
\usage{
alleleDosageEstimation(vcf.file.name,
                        ploidy = 6,
                        min.dp = 10,
                        max.dp = 1000,
                        max.miss = 0.5,
                        max.freq = 0.95,
                        round.up = 1,
                        cut.off = 0,
                        read.err.prob = 0.001)
}
\arguments{
  \item{vcf.file.name}{
Character string of VCF file name.
}
  \item{ploidy}{
Ploidy level of species analyzed. Default is 6.
}
  \item{min.dp}{
Minimum read depth. If a genotype has a read depth less than min.dp, it is assigned NA. Default is 10.
}
  \item{max.dp}{
Maximum read depth. If a genotype has a read depth greater than max.dp, it is assigned NA. Default is 1000.
}
  \item{max.miss}{
Maximum proportion of missing genotype data. 0 allows sites that are completely missing. 1 indicates no missing data allowed.  Default is 0.5.
}
  \item{max.freq}{
Maximum proportion of major genotype (dosage) frequency. 1 indicates monomorphic markers allowed. Default is 0.95.
}
  \item{round.up}{
In allele dosage estimation, a dosage greater than this value is assigned 1, and the others are assigned 0. Default is 1.
}
  \item{cut.off}{
In allele dosage estimation, dosages less than this value is assigned 0. Default is 0.
}
  \item{read.err.prob}{
This value is set to handle unknown experimental errors that breaks calculations within this function. Default is 0.001.
}

}
\details{
Read depth (DP field in VCF) and reference allele depth (RD field) are used as input values for simple binomial probability function to calculate prior probabilities of allele dosage. Then, relative probabilities for all possible allele dosage (i.e. 0, 1/N, 2/N, ..., N/N, where N is the ploidy level) are returned as the result. It should be noted that the method used in this naive and therefore, other elaborated methods such as updog and polyRAD are recommended for more accurate allele dosage estimation. It should also be noted that accuracy of allele dosage estimation has little effect on subsequent association analysis (alleleDosageGLM).
}
\value{
Returns R binary files for estimated allele dosage information (XXX_Geno.Rda) and the chromosomal location information (XXX_Map.Rda). The estimated allele dosage information is a list vector whose each element is n x p matrix of probabilities for allele dosages. n is number of individuals in the population and p is ploidy + 1. The chromosomal location information is a matrix whose row element is SNP marker and the column elements are SNP marker name,  chromosome and position.
}
\author{
Eiji Yamamoto
}
