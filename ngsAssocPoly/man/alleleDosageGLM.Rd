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
Specifies mode of input genotype. "dosage" is for matricies of estimated allele dosage probabilities. "diploidized" is for genotypes in GT field of VCF files. "continuous" is for continuous dosage genotypes calculated from estimated allele dosage probabilities.
}
  \item{Plot}{
IF TRUE, qq and Manhattan plots are generated.
}
  \item{verbose}{
If TRUE, print the run information to the console.
}
}
\details{
Association analysis based on estimated allele dosage. If the trait is binary, binomial function is used. Otherwise, gaussian is used.
}
\value{
Returns a data csv file where the first three columns are the marker name, chromosome, and position, and subsequent columns are the marker scores (-log_{10}p) for the traits, and qq and Manhattan plots are generated in newly created directly.
}
\references{
Yamamoto et al. 2020. Genetic Mapping in Autohexaploid Sweet Potato with Low-Coverage NGS-based Genotyping Data.  G3 (Bethesda) doi: 10.1534/g3.120.401433.
}
\author{
Eiji Yamamoto
}
