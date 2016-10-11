# Read in the 
d1 = read.csv('clpMerge_21nov.csv',header=T,row.names=1)

# Identifier mapping
library(org.Hs.eg.db)
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])

library(MASS) # standard, no need to install
library(class)	# standard, no need to install
library(cluster)	
library(impute)# install it for imputing missing value

suppressMessages(library(getopt))

spec = matrix(c(
  'basedir', 'b', 1, 'character',
  'help', 'h', 0, 'logical'
), byrow=TRUE, ncol=4)

opt <- getopt(spec)

if (is.null(opt$basedir) || !is.null(opt$help)) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}
loc1 = opt$basedir

#########################
## GBM TCGA Biclusters ##
#########################

# Get gene identifiers for expression data
geneExp = d1[which(c(sapply(rownames(d1), function(x){ (strsplit(x,split=':')[[1]])[2] }))=='GEXP'),]
genes = c(sapply(rownames(geneExp), function(x){ (strsplit(x,split=':')[[1]])[3] }))

## Get the mutations and CNVs
mutations = as.matrix(d1[which(c(sapply(rownames(d1), function(x){ (strsplit(x,split=':')[[1]])[2] }))=='GNAB'),])
# Coding potential
mutations_code_potential = as.matrix(d1[intersect(which(c(sapply(rownames(d1), function(x){ (strsplit(x,split=':')[[1]])[2] }))=='GNAB'),which(c(sapply(rownames(d1), function(x){ (strsplit(x,split=':')[[1]])[8] }))=='code_potential_somatic')),])
mut2_code_potential = sapply(rownames(mutations_code_potential), function(mut1) { mean(as.numeric(mutations_code_potential[mut1,]), na.rm=T) })
gt0_05maf_code_potential = rownames(mutations_code_potential)[which(mut2_code_potential>=0.05)]
# Non-silent
mutations_nonsilent = as.matrix(d1[intersect(which(c(sapply(rownames(d1), function(x){ (strsplit(x,split=':')[[1]])[2] }))=='GNAB'),which(c(sapply(rownames(d1), function(x){ (strsplit(x,split=':')[[1]])[8] }))=='nonsilent_somatic')),])
mut2_nonsilent = sapply(rownames(mutations_nonsilent), function(mut1) { mean(as.numeric(mutations_nonsilent[mut1,]), na.rm=T) })
gt0_05maf_nonsilent = rownames(mutations_nonsilent)[which(mut2_nonsilent>=0.05)]
# Take all 35 of the code_potential and the 88 of the nonsilent for the NCI thingies
mutations2 = mutations[c(gt0_05maf_code_potential,gt0_05maf_nonsilent[36:123]),]
rownames(mutations2) = as.character(sapply(rownames(mutations2), function(x) { paste(gsub('-','_',strsplit(x,":")[[1]][3]), strsplit(x,":")[[1]][8],sep='_') }))
mutations2 = matrix(ncol=ncol(mutations2),nrow=nrow(mutations2),data=as.numeric(mutations2),dimnames=dimnames(mutations2))

## Get expression data for all regulators
tf1 = read.csv('humanTFsFINAL_GO_0003700.csv')[,1]
tf_genes = names(genes)[which(genes %in% tf1)]
tfExp = as.matrix(d1[tf_genes,])
rownames(tfExp) = as.character(sapply(rownames(tfExp), function(x) { gsub('-','_',strsplit(x,":")[[1]][3]) }))
miRNAExp = as.matrix(d1[which(c(sapply(rownames(d1), function(x){ (strsplit(x,split=':')[[1]])[2] }))=='MIRN'),])
rownames(miRNAExp) = as.character(sapply(rownames(miRNAExp), function(x) { gsub('-','_',strsplit(x,":")[[1]][3]) }))
regExp = rbind(tfExp, miRNAExp)
regExp = matrix(ncol=ncol(regExp), nrow=nrow(regExp), data=as.numeric(regExp), dimnames=dimnames(regExp))

## Load bicluster eigengenes
be1 = read.csv(paste(loc1,'sygnal/output/biclusterFirstPrincComponents.csv',sep=''),row.names=1, header=T)
rownames(be1) = paste('bic',rownames(be1),sep='_')
ol1 = intersect(colnames(be1),colnames(mutations2))
d2 = rbind(as.matrix(mutations2[,ol1]), as.matrix(regExp[,ol1]), as.matrix(be1[,ol1]))
d3 = t(na.omit(t(d2)))

## Use parallel processing to calculate t-test p-values and fold-changes faster
library(doParallel)
registerDoParallel(cores=4)
if(!file.exists(paste(loc1,'sygnal/output/mut_reg_t_test_p_values.csv',sep=''))) {
    m1 = foreach(mut1=rownames(mutations2), .combine=rbind) %dopar% sapply(rownames(regExp), function(reg1) { t.test(as.numeric(d3[reg1,]) ~ as.numeric(d3[mut1,]))$p.value } )
    rownames(m1) = rownames(mutations2)
    fc1 = foreach(mut1=rownames(mutations2), .combine=rbind) %dopar% sapply(rownames(regExp), function(reg1) { median(2^as.numeric(d3[reg1,which(d3[mut1,]==0)]))/median(2^as.numeric(d3[reg1,which(d3[mut1,]==1)])) } )
    rownames(fc1) = rownames(mutations2)
    write.csv(m1,paste(loc1,'sygnal/output/mut_reg_t_test_p_values.csv',sep=''))
    write.csv(fc1,paste(loc1,'sygnal/output/mut_reg_fold_changes.csv',sep=''))
} else {
    m1 = read.csv(paste(loc1,'sygnal/output/mut_reg_t_test_p_values.csv',sep=''),header=T,row.names=1)
    fc1 = read.csv(paste(loc1,'sygnal/output/mut_reg_fold_changes.csv',sep=''),header=T,row.names=1)
}

## Select which mutations are associated with which regualtors
# Use Student's T-test p-value cutoff of 0.05 and fold change of FC>=1.2 or FC<=0.8 as combined cutoffs
sigRegFC = sapply(rownames(m1), function(x) { colnames(m1)[intersect(which(m1[x,]<=0.05),union(which(fc1[x,]>=1.25),which(fc1[x,]<=0.8)))] })


########################
## Causality analysis ##
########################
## Use for filtering:
#  1. Signficant differntial expression of regulator between wt and mutant (FC <= 0.8 or FC >= 1.25, and T-test p-value <= 0.05)
source('neoDecember2015.R')
registerDoParallel(12)
dir.create(paste(loc1, 'sygnal/output/causality', sep=''))
for(mut1 in names(sigRegFC)) {
    # Make a place to store out the data from the analysis
    mut2 = mut1
    if(nchar(mut2)>75) {
        mut2 = substr(mut2,1,75)
    }
    dir.create(paste(loc1, 'sygnal/output/causality/causal_', mut2, sep=''))

    # Change the names to be compatible with NEO
    print(paste('Starting ',mut1,'...',sep=''))
    foreach(reg1=sigRegFC[[mut1]]) %dopar% {
        # Make the data matrix will all genes strsplit(mut1)[[1]][1]
        d3 = t(na.omit(t(d2[c(mut1,reg1,rownames(be1)),])))
        dMut1 = matrix(data=as.numeric(d3),nrow=dim(d3)[1],ncol=dim(d3)[2],byrow=F,dimnames=dimnames(d3))
        print(paste('  Starting ',mut1,' vs. ', reg1,' testing ', length(rownames(be1)), ' biclusters...', sep=''))
        sm1 = try(single.marker.analysis(t(dMut1),1,2,3:length(rownames(dMut1))),silent=TRUE)
        if (!(class(sm1)=='try-error')) {
            write.csv(sm1[order(sm1[,6],decreasing=T),],paste(loc1, 'sygnal/output/causality/causal_', mut2, '/sm.nonsilent_somatic.',mut2,'_',reg1,'.csv',sep=''))
            print(paste('Finished ',reg1,'.',sep=''))
        } else {
            print(paste('  Error ',mut1,'.',sep=''))
        }
    }
    print(paste('Finished ',mut1,'.',sep=''))
}

