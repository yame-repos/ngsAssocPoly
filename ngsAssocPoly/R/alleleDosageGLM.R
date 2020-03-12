alleleDosageGLM <- function(Geno.name,
                            Map.name,
                            pheno.file.name,
                            method = "dogmat")
{
  glmPoly = function(pheno, geno, method = method) {
    F.ratio = function(S0, S1, ndf, ddf) {
      if (S0 < S1) {
        score = 0
      } else {
        Fr = ((S0 - S1) / ndf) / (S1 / ddf)
        logp = pf(Fr, ndf, ddf, lower.tail=FALSE, log.p=TRUE)
        score = -1 * logp / log(10)
      }
      return(score)
    }
    make.design.mat <- function(row.id, col.id) {
      Z <- matrix(0, nrow=length(row.id), ncol=length(col.id))
      rownames(Z) <- row.id
      colnames(Z) <- col.id
      Z[cbind(1:nrow(Z),match(row.id,col.id))] <- 1
      return(Z)
    }
    n.markers <- length(geno)
    phenos = colnames(pheno)[-1]
    results = as.data.frame(matrix(NA, nrow=n.markers, ncol=length(phenos)))
    colnames(results) = phenos
    for (phenoi in phenos) {
      y = pheno[[phenoi]]
      names(y) <- pheno$gid
      y = y[!is.na(y)]
      if (length(unique(y))==2) glm.exe = function(y, x) glm(y~x, family=binomial)
      if (length(unique(y))==2) glm.exe = function(y, x) glm(y~x, family=gaussian)

      scores <- rep(NA, n.markers)
      for (p in 1:n.markers) {
        print(paste(phenoi, "::", p,"/",n.markers))
        X = geno[[p]]
        X = X[!is.na(X[,1]),]
        use = intersect(names(y), rownames(X))
        if (length(use) < 50) next
        yp = y[is.element(names(y), use)]
        if (method == "dogmat") {
          Xp = X[is.element(rownames(X), use),,drop=FALSE]
          svd.xp = svd(Xp)
          r = max(which(svd.xp$d > 1e-08))
          xp = as.matrix(svd.xp$u[,1:r])
          rownames(xp) = rownames(Xp)

        } else if (method == "continuous") {
          X = X %*% as.numeric(colnames(X))
          xp = X[is.element(rownames(X), use), , drop=FALSE]
        } else if (method == "diplodized") {
          X = X %*% as.numeric(colnames(X))
          xp = X[is.element(rownames(X), use), , drop=FALSE]
          xp[xp == 0] = 0
          xp[(xp != 0) & (xp < 6)] = 1
          xp[xp == 6] = 2
        } else {
          stop("method must be dogmat, continuous or diplodized.")
        }
        Z = make.design.mat(names(yp), rownames(xp))
        Zxp = Z %*% xp
        res = glm.exe(yp, Zxp)
        res1 = yp - res$fitted
        S1 = as.double(res1 %*% res1)
        S0 = var(yp) * (length(yp) - 1)
        ndf = ncol(Zxp)
        ddf = length(yp) - ncol(Zxp) - 1
        scores[p] = F.ratio(S0, S1, ndf, ddf)
      }
      results[[phenoi]] = scores
    }
    return(results)
  }

  if (!file.exists(Map.name)) stop(paste(Map.name, "does not exist."))
  if (!file.exists(Geno.name)) stop(paste(Geno.name, "does not exist."))
  if (!file.exists(pheno.file.name)) stop(paste(pheno.file.name, "does not exist."))

  project.id <- strsplit(Geno.name, "_Geno.rda")[[1]]

  load(Map.name)
  load(Geno.name)
  phenos = read.csv(pheno.file.name)
  pheno.names = colnames(phenos)[-1]

  scores.name = paste(project.id, "_scores.csv", sep="")
  scores = glmPoly(phenos, Geno, method)
  scores = cbind(Map, scores)
  write.csv(scores, scores.name, row.names = FALSE)

  thr.line = -log10(0.05/length(Geno))
  n.LG = length(unique(Map$chrom))
  col.list = rep(RColorBrewer::brewer.pal(8, "Dark2"), floor(n.LG/8))[1:n.LG]
  col.list[n.LG] = "black"
  rm(n.LG); gc()

  #new.dir.name = paste(project.id, "_Plots", sep="")
  #system(paste("mkdir", new.dir.name))
  #setwd(new.dir.name)
  pheno.names = colnames(scores)[-c(1:3)]
  for (phenoi in pheno.names) {
    p = scores[[phenoi]]
    input = data.frame(SNP=Map$marker,
                       CHR=Map$chrom,
                       BP=Map$pos,
                       P=10^(-1*p))
    input = input[order(input$BP),]
    input = input[order(input$CHR),]
    input$CHR = as.numeric(input$CHR)
    input = input[!is.na(input$P),,drop=FALSE]
    manhattan.name = paste(project.id, phenoi, "_manhattan.jpeg",sep="_")
    jpeg(filename=manhattan.name,width=960,height=640,pointsize=28,quality=100)
    qqman::manhattan(input,
                     col = col.list,
                     main = phenoi,
                     suggestiveline = Inf,
                     genomewideline = thr.line)
    dev.off()
    qq.name = paste(project.id, phenoi, "_qq.jpeg",sep="_")
    jpeg(filename=qq.name,width=640,height=640,pointsize=28,quality=100)
    qqman::qq(input$P, main=phenoi)
    dev.off()
  }
  #setwd("..")

}
