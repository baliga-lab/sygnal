#! /usr/bin/env Rscript
suppressMessages(library(methods))
suppressMessages(library(MASS)) # standard, no need to install
suppressMessages(library(class))	# standard, no need to install
suppressMessages(library(cluster))
suppressMessages(library(impute)) # install it for imputing missing value
suppressMessages(library(getopt))

spec = matrix(c(
  'outdir', 'o', 1, 'character',
  'ratios', 'r', 1, 'character',
  'mirnas', 'm', 1, 'character',
  'som_muts', 's', 1, 'character',
  #'pat_seq', 'p', 1, 'character',
  'eigengene', 'e', 1, 'character',
  'tumor', 't', 1, 'character',
  'cores', 'c', 1, 'integer',
  'help', 'h', 0, 'logical'
  ), byrow=TRUE, ncol=4)

opt <- getopt(spec)

if (is.null(opt$outdir) || is.null(opt$tumor) || is.null(opt$ratios) || is.null(opt$mirnas) || is.null(opt$som_muts) || is.null(opt$eigengene) || is.null(opt$cores) || !is.null(opt$help)) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

cat('Running NEO...')

########################
## NEO for TCGA Tumor ##
########################

# Load gene expression and extract TFs
cat('\nLoading ratios...')
ratios <- read.csv(file=opt$ratios, as.is=T, header=T, row.names=1 )
ratios = ratios[,which(sapply(colnames(ratios), function(x) { sum(is.na(ratios[,x])) })!=length(rownames(ratios)))]
ratios = ratios[which(rowSums(ratios)!=0),]
tf1 = read.csv('NEO/humanTFsFINAL_ENTREZ_GO_0003700.csv')[,1]
tf_genes = rownames(ratios)[which(rownames(ratios) %in% tf1)]
tfExp = as.matrix(ratios[tf_genes,])
rownames(tfExp) = paste('X',tf_genes,sep='')

# Load miRNA expression
if (opt$mirnas != 'NA') {
    cat('\nLoading miRatios...')
    miRExp <- read.csv(file=opt$mirnas, as.is=T, header=T, row.names=1 )
    miRExp = miRExp[,which(sapply(colnames(miRExp), function(x) { sum(is.na(miRExp[,x])) })!=length(rownames(miRExp)))]
    miRExp = miRExp[which(rowSums(miRExp)!=0),]

    # Merge regulator expression
    cat('\nLoading regExp...')
    comSamp = intersect(colnames(tfExp),colnames(miRExp))
    regExp = rbind(tfExp[,comSamp], miRExp[,comSamp])
} else {
    regExp = tfExp
}

# Load somatic mutations
cat('\nLoading somatic mutations...')
mutations1 <- read.csv(file=opt$som_muts, as.is=T, header=T, row.names=1)
maf1 = rowSums(mutations1)
mutations2 = mutations1[names(maf1)[which(maf1>=3)],]

# Load bicluster eigengenes
cat('\nLoading bicluster eigengene...')
be1 = read.csv(file=opt$eigengene, row.names=1, header=T)
rownames(be1) = paste('bic',rownames(be1),sep='_')
ol1 = intersect(intersect(colnames(be1),colnames(mutations2)),colnames(regExp))
cat(paste('\nsom_muts = ',nrow(mutations2),'; regExp = ',nrow(regExp),'; be1 = ',nrow(be1),sep=''))
d2 = rbind(as.matrix(mutations2[,ol1]), as.matrix(regExp[,ol1]), as.matrix(be1[,ol1]))
d3 = t(na.omit(t(d2)))

# Error handling t.test function
t.test_err <- function(...) {
    obj = try(t.test(...), silent=TRUE)
    if (is(obj, 'try-error')) return(1) else return(obj$p.value)
}

## Use parallel processing to calculate t-test p-values and fold-changes faster
cat('\nCalculating associations...')
suppressMessages(library(doParallel))
registerDoParallel(cores=opt$cores)
if(!file.exists(paste('output/mut_reg_t_test_p_values_',opt$tumor,'.csv',sep=''))) {
    m1 = foreach(mut1=rownames(mutations2), .combine=rbind) %dopar% sapply(rownames(regExp), function(reg1) { t.test_err(as.numeric(d3[reg1,]) ~ as.numeric(d3[mut1,])) } )
    rownames(m1) = rownames(mutations2)
    fc1 = foreach(mut1=rownames(mutations2), .combine=rbind) %dopar% sapply(rownames(regExp), function(reg1) { median(2^as.numeric(d3[reg1,which(d3[mut1,]==0)]),na.rm=T)/median(2^as.numeric(d3[reg1,which(d3[mut1,]==1)]),na.rm=T) } )
    rownames(fc1) = rownames(mutations2)
    write.csv(m1,paste('output/mut_reg_t_test_p_values_',opt$tumor,'.csv',sep=''))
    write.csv(fc1,paste('output/mut_reg_fold_changes_',opt$tumor,'.csv',sep=''))
} else {
    m1 = read.csv(paste('output/mut_reg_t_test_p_values_',opt$tumor,'.csv',sep=''),header=T,row.names=1)
    fc1 = read.csv(paste('output/mut_reg_fold_changes_',opt$tumor,'.csv',sep=''),header=T,row.names=1)
}

## Select which mutations are associated with which regualtors
# Use Student's T-test p-value cutoff of 0.05 and fold change of FC>=1.2 or FC<=0.8 as combined cutoffs
cat('\nSelecting which associations to test...')
sigRegFC = sapply(rownames(m1), function(x) { colnames(m1)[intersect(which(m1[x,]<=0.05),union(which(fc1[x,]>=1.25),which(fc1[x,]<=0.8)))] })


########################
## Causality analysis ##
########################
## Use for filtering:
#  1. Signficant differntial expression of regulator between wt and mutant (FC <= 0.8 or FC >= 1.25, and T-test p-value <= 0.05)
cat('\nRunning NEO...')
source('NEO/neoDecember2015.R')
registerDoParallel(cores=opt$cores)
dir.create(paste('output/causality_',opt$tumor,sep=''), showWarnings=F)
foreach(mut1=names(sigRegFC)) %dopar% {
    # Make a place to store out the data from the analysis
    mut2 = mut1
    if(nchar(mut2)>75) {
        mut2 = substr(mut2,1,75)
    }
    dir.create(paste('output/causality_',opt$tumor,'/causal_', mut2, sep=''))

    # Change the names to be compatible with NEO
    print(paste('Starting ',mut1,'...',sep=''))
    for(reg1 in sigRegFC[[mut1]]) {
        # Make the data matrix will all genes strsplit(mut1)[[1]][1]
        d3 = t(na.omit(t(d2[c(mut1,reg1,rownames(be1)),])))
        dMut1 = matrix(data=as.numeric(d3),nrow=dim(d3)[1],ncol=dim(d3)[2],byrow=F,dimnames=dimnames(d3))
        print(paste('  Starting ',mut1,' vs. ', reg1,' testing ', length(rownames(be1)), ' biclusters...', sep=''))
        sm1 = try(single.marker.analysis(t(dMut1),1,2,3:length(rownames(dMut1))),silent=TRUE)
        if (!(class(sm1)=='try-error')) {
            write.csv(sm1[order(sm1[,6],decreasing=T),],paste('output/causality_',opt$tumor,'/causal_', mut2, '/sm.nonsilent_somatic.',mut2,'_',reg1,'.csv',sep=''))
            print(paste('Finished ',reg1,'.',sep=''))
        } else {
            print(paste('  Error ',mut1,'.',sep=''))
        }
    }
    print(paste('Finished ',mut1,'.',sep=''))
}
cat('\nDone!')
