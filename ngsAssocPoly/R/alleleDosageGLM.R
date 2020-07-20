alleleDosageGLM = function (Geno.name,
                            Map.name,
                            Pheno.file.name,
                            method = "dogmat",
                            Plot = TRUE,
                            verbose = TRUE)
{
  make.design.mat = function(row.id, col.id) {
    Z = matrix(0, nrow = length(row.id), ncol = length(col.id))
    rownames(Z) = row.id
    colnames(Z) = col.id
    Z[cbind(1:nrow(Z), match(row.id, col.id))] = 1
    return(Z)
  }
  glmPoly = function(pheno, Geno, method = method) {
    n.markers = length(Geno)
    y = pheno[, 2]
    names(y) = pheno$gid
    if (length(unique(y)) == 2) glm.exe = function(y, x) glm(y ~ x, family = binomial)
    if (length(unique(y)) != 2) glm.exe = function(y, x) glm(scale(y) ~ x, family = gaussian)
    scores = rep(0, n.markers)
    for (p in 1:n.markers) {
      if (verbose) print(paste(colnames(pheno)[2], "::", p, "/", n.markers))
      X = Geno[[p]]
      X = X[!is.na(X[, 1]), ]
      use = intersect(names(y), rownames(X))
      if (length(use) < 50) next
      yp = y[is.element(names(y), use)]
      if (method == "dogmat") {
        Xp = X[is.element(rownames(X), use), , drop = FALSE]
        svd.xp = svd(Xp)
        r = max(which(svd.xp$d > 1))
        xp = as.matrix(svd.xp$u[, 1:r])
        rownames(xp) = rownames(Xp)
      }
      else if (method == "continuous") {
        X = X %*% as.numeric(colnames(X))
        xp = X[is.element(rownames(X), use), , drop = FALSE]
      }
      else if (method == "diploidized") {
        X = X %*% as.numeric(colnames(X))
        xp = X[is.element(rownames(X), use), , drop = FALSE]
        xp[xp < 1] = 0
        xp[(xp != 0) & (xp <= (ncol(Geno[[1]]) - 1))] = 1
        xp[xp > (ncol(Geno[[1]]) - 1)] = 2
      }
      else {
        stop("method must be dogmat, continuous or diplodized.")
      }
      Z = make.design.mat(names(yp), rownames(xp))
      Zxp = Z %*% xp
      soln = glm.exe(yp, Zxp)
      scores[p] = -log10(pchisq(soln$null.deviance - soln$deviance, df = soln$df.null - soln$df.residual, lower.tail = FALSE))
    }
    return(scores)
  }
  if (!file.exists(Map.name)) stop(paste(Map.name, "does not exist."))
  if (!file.exists(Geno.name)) stop(paste(Geno.name, "does not exist."))
  if (!file.exists(Pheno.file.name)) stop(paste(Pheno.file.name, "does not exist."))
  project.id = strsplit(Geno.name, "_Geno.rda")[[1]]
  load(Map.name)
  load(Geno.name)
  Pheno = read.csv(Pheno.file.name)
  pheno.names = colnames(Pheno)[-1]
  n.phenos = length(pheno.names)
  n.markers = length(Geno)
  results.name = paste(project.id, "_scores.csv", sep = "")
  Results = matrix(NA, nrow = n.markers, ncol = n.phenos)
  colnames(Results) = pheno.names
  for (i in 1:n.phenos) {
    pheno = Pheno[, is.element(colnames(Pheno), c("gid", pheno.names[i]))]
    pheno = pheno[!is.na(pheno[, 2]), ]
    Results[, i] = glmPoly(pheno, Geno, method)
  }
  Results = cbind(Map, Results)
  write.csv(Results, results.name, row.names = FALSE)
  if (Plot) {
    n.LG = length(unique(Map$chrom))
    col.list = rep(RColorBrewer::brewer.pal(8, "Dark2"), floor(n.LG/8))[1:n.LG]
    rm(n.LG)
    gc()
    for (phenoi in pheno.names) {
      p = Results[[phenoi]]
      input = data.frame(SNP = Map$marker, CHR = Map$chrom, BP = Map$pos, P = 10^(-1 * p))
      input = input[order(input$BP), ]
      input = input[order(input$CHR), ]
      input$CHR = as.numeric(as.factor(input$CHR))
      input = input[!is.na(input$P), , drop = FALSE]
      input = input[input$P != 0, , drop = FALSE]
      manhattan.name = paste(project.id, phenoi, "_manhattan.jpeg", sep = "_")
      jpeg(filename = manhattan.name, width = 960, height = 480, pointsize = 28, quality = 100)
      qqman::manhattan(input, col = col.list, main = phenoi, suggestiveline = Inf, genomewideline = Inf)
      dev.off()
      qq.name = paste(project.id, phenoi, "_qq.jpeg", sep = "_")
      jpeg(filename = qq.name, width = 480, height = 480, pointsize = 28, quality = 100)
      qqman::qq(input$P, main = phenoi)
      dev.off()
    }
  }
}
