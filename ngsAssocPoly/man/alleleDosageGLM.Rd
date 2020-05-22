\name{alleleDosageGLM}
\alias{alleleDosageGLM}
\title{
Association analysis based on estimated allele dosage
}
\description{
Performs association analysis  based on estimated allele dosage with generalized linear model (GLM).
}
\usage{
alleleDosageGLM(Geno.name,
                Map.name,
                Pheno.file.name,
                method = "dogmat",
                Plot = TRUE,
                verbose = TRUE)
}
\arguments{
  \item{Geno.name}{
Character string of file name of estimated allele dosage that is created by function alleleDosageEstimation.
}
  \item{Map.name}{
Character string of file name of map information for Geno.name that is created by function alleleDosageEstimation.
}
  \item{Pheno.file.name}{
Character string of CSV file name that include phenotype data. The first colmn must be the line (or individual) name that corresponds to Geno.name.
}
  \item{method}{
Specifies mode of input genotype. "dogmat" is for matricies of estimated allele dosage probabilities. "continuous" is for vectors of continuous allele dosage values. "diplodized" is for genotypes in GT field of VCF files.
}
  \item{Plot}{
If TRUE, qq and Manhattan plots are generated, default TRUE
}
  \item{verbose}{
If TRUE, the state of progress is printed, default TRUE
}
}
\details{
Association analysis based on estimated allele dosage. If the trait is binary, binomial function is used. Otherwise, gaussian is used.
}
\value{
Returns a data csv file where the first three columns are the marker name, chromosome, and position, and subsequent columns are the marker scores (-log_{10}p) for the traits, and qq and Manhattan plots are generated in newly created directly.
}
\author{
Eiji Yamamoto
}
