#!/usr/bin/env Rscript

suppressMessages(library(WGCNA))
suppressMessages(library(multicore))
suppressMessages(library(getopt))

allowWGCNAThreads()

# read command line arguments
spec = matrix(c(
  'ratios','r',1,'character',
  'outdir', 'o', 1, 'character',
  'help', 'h', 0, 'logical'
), byrow=TRUE, ncol=4)

opt <- getopt(spec)

if (is.null(opt$ratios) || is.null(opt$outdir) || !is.null(opt$help)) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

residual <- function( rats ) {
  d.rows <- rowMeans( rats, na.rm=T )
  d.cols <- colMeans( rats, na.rm=T )
  d.all <- mean( d.rows, na.rm=T )
  rij <- rats + d.all
  rij <- rij - matrix( d.cols, nrow=nrow( rij ), ncol=ncol( rij ), byrow=T )
  rij <- rij - matrix( d.rows, nrow=nrow( rij ), ncol=ncol( rij ), byrow=F )
  average.r <- mean( abs( rij ), na.rm = TRUE )
  average.r
}

residual.norm <- function( rats, maxRowVar ) {
  d.rows <- rowMeans( rats, na.rm=T )
  d.cols <- colMeans( rats, na.rm=T )
  d.all <- mean( d.rows, na.rm=T )
  rij <- rats + d.all
  rij <- rij - matrix( d.cols, nrow=nrow( rij ), ncol=ncol( rij ), byrow=T )
  rij <- rij - matrix( d.rows, nrow=nrow( rij ), ncol=ncol( rij ), byrow=F )
  average.r <- mean( abs( rij ), na.rm = TRUE )
  row.var <- mean( apply( rats, 1, var, use = "pairwise.complete.obs" ), na.rm=T )
  if ( is.na( row.var ) || row.var > maxRowVar ) row.var <- maxRowVar
  average.r <- average.r / row.var
  average.r
}
getEigengene <- function (expr, rows, impute = TRUE, nPC = 1, align = "along average", 
    excludeGrey = FALSE, grey = ifelse(is.numeric(rows), 0, 
        "grey"), subHubs = FALSE, trapErrors = FALSE, returnValidOnly = trapErrors, 
    softPower = 1, scale = TRUE, verbose = 0, indent = 0) 
{
    nVarExplained = nPC
    modlevels = 1:length(rows)
    PrinComps = data.frame(matrix(NA, nrow = dim(expr)[[1]], 
        ncol = length(modlevels)))
    averExpr = data.frame(matrix(NA, nrow = dim(expr)[[1]], ncol = length(modlevels)))
    varExpl = data.frame(matrix(NA, nrow = nVarExplained, ncol = length(modlevels)))
    validMEs = rep(TRUE, length(modlevels))
    validAEs = rep(FALSE, length(modlevels))
    isPC = rep(TRUE, length(modlevels))
    isHub = rep(FALSE, length(modlevels))
    validColors = rows
    names(PrinComps) = paste(moduleColor.getMEprefix(), modlevels, 
        sep = "")
    names(averExpr) = paste("AE", modlevels, sep = "")
    for (i in c(1:length(modlevels))) {
        modulename = modlevels[i]
        restrict1 = rows[[modulename]]
        datModule = as.matrix(t(expr[ ,restrict1]))
        n = dim(datModule)[1]
        p = dim(datModule)[2]
        pc = try({
            if (scale) 
                datModule = t(scale(t(datModule)))
            svd1 = svd(datModule, nu = min(n, p, nPC), nv = min(n, 
                p, nPC))
            veMat = cor(svd1$v[, c(1:min(n, p, nVarExplained))], 
                t(datModule), use = "p")
            varExpl[c(1:min(n, p, nVarExplained)), i] = apply(veMat^2, 
                1, mean, na.rm = TRUE)
            svd1$v[, 1]
        }, silent = TRUE)
        PrinComps[, i] = pc
        ae = try({
            if (isPC[i]) 
                scaledExpr = scale(t(datModule))
            averExpr[, i] = apply(scaledExpr, 1, mean, na.rm = TRUE)
            if (align == "along average") {
                if (cor(averExpr[, i], PrinComps[, i], use = "p") < 
                0) 
                PrinComps[, i] = -PrinComps[, i]
            }
            0
        }, silent = TRUE)
    }
    allOK = (sum(!validMEs) == 0)
    if (returnValidOnly && sum(!validMEs) > 0) {
        PrinComps = PrinComps[, validMEs]
        averExpr = averExpr[, validMEs]
        varExpl = varExpl[, validMEs]
        validMEs = rep(TRUE, times = ncol(PrinComps))
        isPC = isPC[validMEs]
        isHub = isHub[validMEs]
        validAEs = validAEs[validMEs]
    }
    allPC = (sum(!isPC) == 0)
    allAEOK = (sum(!validAEs) == 0)
    list(eigengenes = PrinComps, averageExpr = averExpr, varExplained = varExpl, 
        nPC = nPC, validMEs = validMEs, validColors = validColors, 
        allOK = allOK, allPC = allPC, isPC = isPC, isHub = isHub, 
        validAEs = validAEs, allAEOK = allAEOK)
}

# Read in genes for each cluster
infile.gene <- paste(opt$outdir, 'cluster.members.genes.txt', sep='/')
d1 = read.csv(infile.gene,header=F)
biclustMembership.gene = list()
allGenes = c()
for(j in 1:length(d1[,1])) {
    biclustMembership.gene[[j]] = strsplit(as.character(d1[j,]),split=' ')[[1]][-1]
}

# Read in conditions for each cluster
infile.cond <- paste(opt$outdir, 'cluster.members.conditions.txt', sep='/')
d1 = read.csv(infile.cond,header=F)
biclustMembership.cond = list()
allGenes = c()
for(j in 1:length(d1[,1])) {
    biclustMembership.cond[[j]] = strsplit(as.character(d1[j,]),split=' ')[[1]][-1]
}

# Read in expression ratios file
ratios <- read.delim(file=opt$ratios, sep="\t", as.is=T, header=T,row.names=1 )
#rownames(ratios) <- toupper(rownames(ratios))
maxRowVar = mean( apply( ratios, 1, var, use="pair" ), na.rm=T )

# Calculate the residuals for all clusters in the second dataset
ks = length(biclustMembership.gene)
outNames = c('n.rows','n.cols','orig.resid','avg.perm.resid','perm.p','orig.resid.norm','avg.norm.perm.resid','norm.perm.p','pc1.var.exp','avg.pc1.var.exp','pc1.perm.p')
m1 = matrix(ncol=length(outNames),nrow=ks,dimnames=list(1:ks,outNames))
#permutations = 1/0.05*ks
permutations = 1000
print(paste('Running ',permutations,' permutations...',sep=''))
m1 = do.call(rbind, mclapply(1:ks, function(k) {
        r1 = rep(NA,12)
        r1[1] = k
        # Get and add number of rows and columns
        k.rows = biclustMembership.gene[[k]]
        k.cols = biclustMembership.cond[[k]]
        if(length(k.rows)>1 && length(k.cols)>1) {
            r1[2] = length(k.rows)
            r1[3] = length(k.cols)
            r1[4] = residual(as.matrix(ratios[k.rows,k.cols]))
            sub = sapply(1:permutations, function(i) { residual(as.matrix(ratios[sample(rownames(ratios),r1[2]), sample(colnames(ratios),r1[3])])) })
            r1[5] = mean(sub)
            r1[6] = length(which(sub <= r1[4]))/permutations
            r1[7] = residual.norm(as.matrix(ratios[k.rows,k.cols]),maxRowVar)
            sub = sapply(1:permutations, function(i) { residual.norm(as.matrix(ratios[sample(rownames(ratios),r1[2]), sample(colnames(ratios),r1[3])]), maxRowVar) })
            r1[8] = mean(sub)
            r1[9] = length(which(sub <= r1[7]))/permutations
            testEm.rows = list()
            testEm.rows[[1]] = k.rows
            for( i in 2:(permutations+1)) {
                testEm.rows[[i]] = sample(rownames(ratios),r1[2])
            }
            eg1 = getEigengene(t(ratios),testEm.rows)
            var.exp = t(eg1$varExplained)[,1]
            r1[10] = var.exp[1]
            r1[11] = mean(var.exp[2:length(var.exp)],na.rm=TRUE)
            r1[12] = length(which(na.omit(var.exp[2:length(var.exp)]) >= r1[10]))/length(na.omit(var.exp[2:length(var.exp)]))
        }
        print(c(k,as.numeric(r1)))
    }))
outNames = c('bicluster','bicluster','n.rows','n.cols','orig.resid','avg.perm.resid','perm.p','orig.resid.norm','avg.norm.perm.resid','norm.perm.p','pc1.var.exp','avg.pc1.var.exp','pc1.perm.p')
colnames(m1) = outNames
outfile <- paste(opt$outdir, 'residualPermutedPvalues_permAll.csv', sep='/')
write.csv(m1,file=outfile)

