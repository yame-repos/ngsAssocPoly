alleleDosageEstimation <- function(vcf.file.name,
                                  ploidy = 6,
                                  min.dp = 10,
                                  max.dp = 1000,
                                  max.miss = 0.5,
                                  max.freq = 0.95,
                                  round.up = 1,
                                  cut.off = 0,
                                  read.err.prob = 0.001
)
{
  if (!file.exists(vcf.file.name)) stop(paste(vcf.file.name, "does not exist."))

  make.prob.vec = function(ploidy=ploidy, err=read.err.prob) {
    v = numeric(ploidy+1)
    names(v) = c(0, 1:ploidy)
    v[1] = err
    for (i in 1:ploidy) v[(i+1)] = i/ploidy
    v[(ploidy+1)] = 1 - err
    return(v)
  }
  count.miss = function(vec) {
    length(vec[is.na(vec)])
  }

  print(paste("Input VCF file :", vcf.file.name))

  project.id <- strsplit(vcf.file.name, ".vcf")[[1]][1]
  Map.name = paste(project.id, "_Map.rda", sep="")
  Geno.name = paste(project.id, "_Geno.rda", sep="")

  vcf = vcfR::read.vcfR(vcf.file.name)
  chrom = vcfR::getCHROM(vcf)
  pos = vcfR::getPOS(vcf)
  marker = paste(chrom, pos, sep="_")
  Map = data.frame(marker=marker, chrom=chrom, pos=pos)
  rm(chrom, pos, marker)
  gc()

  prob.vec = make.prob.vec(ploidy, read.err.prob)
  rm(make.prob.vec, read.err.prob)
  gc()

  dp = vcfR::extract.gt(vcf, element="DP")
  mode(dp) = "numeric"
  rd = vcfR::extract.gt(vcf, element="RD")
  mode(rd) = "numeric"
  rm(vcf)
  gc()

  dp[dp > max.dp] = NA
  dp[dp < min.dp] = NA

  miss = apply(dp, 1, count.miss) / ncol(dp)
  use = miss < max.miss
  Map = Map[use,,drop=FALSE]
  dp = dp[use,,drop=FALSE]
  rd = rd[use,,drop=FALSE]
  rm(max.dp, min.dp, count.miss, miss)
  gc()

  n.samples = ncol(dp)
  n.markers = nrow(dp)
  sample.names = colnames(dp)

  Geno = vector(mode="list", length=n.markers)
  names(Geno) = Map$marker
  M.geno = matrix(NA, nrow=n.samples, ncol=ploidy+1)
  rownames(M.geno) = sample.names
  colnames(M.geno) = c(0:ploidy)
  geno.prob = numeric(ploidy+1)

  for (p in 1:n.markers) {
    print(paste("Calculating allelic dosage probability:  ", p, "/", n.markers, "completed"))
    M = M.geno
    for (i in 1:n.samples) {
      n = dp[p,i]
      if (is.na(n)) next
      r = rd[p,i]
      if (is.na(r)) next
      if (r > n) next
      prob = geno.prob
      for (j in 0:ploidy) prob[j+1] = dbinom(r, n, prob.vec[j+1])
      if (sum(prob)==0) next
      prob = prob/sum(prob)
      if (max(prob) >= round.up) prob[prob!=max(prob)] = 0
      prob[prob < cut.off] = 0
      M[i,] = prob/sum(prob)
    }
    Geno[[p]] = M
  }

  rm(rd, dp, M, M.geno, geno.prob, prob, prob.vec,
     n.markers, n.samples, sample.names,
     round.up, cut.off, n, r, i, j, p)
  gc()

  n.markers = length(Geno)
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

  rm(n.markers, freq, good,not.miss, max.miss, max.freq, i, m)
  gc()

  save(Geno, file=Geno.name)
  save(Map, file=Map.name)

  print(paste("Estimated allele dosage information :", Geno.name))
  print(paste("Map information :", Map.name))
  print(paste("The allele dosage estimation for [", vcf.file.name, "] has completed.", sep=""))

}
