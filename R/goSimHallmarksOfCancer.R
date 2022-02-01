########## Hallmarks of Cancer Enrichment from GO Terms #####
#                                                           #
# Copyright 2019 Baliga Lab, Institute for Systems Biology  #
# Original script by Chris Plaisier                         #
# Modified by Serdar Turkarslan October 2019                #
#                                                           #
#############################################################

##### Install required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GOSemSim")
BiocManager::install("org.Hs.eg.db")

##### Load libraries
#library('GOSim')
library('GOSemSim')
library('org.Hs.eg.db')

#####  List of GO Terms assigned to Hallmarks of Cancer. Original list from Plasier et al. 2016 
#####     with update by using Knijnenburg et al., 2015
l1 = list()
l1$SelfSufficiencyInGrowthSignals = unique(
    c(
        'GO:0009967',
        'GO:0030307',
        'GO:0008284',
        'GO:0045787',
        'GO:0007165',
        'GO:0008283',
        'GO:0016049',
        'GO:0007049',
        'GO:0051301',
        'GO:0051781'
    )
)
l1$InsensitivityToAntigrowthSignals = unique(
    c(
        'GO:0045786',
        'GO:0007165',
        'GO:0009968',
        'GO:0008285',
        'GO:0030308',
        'GO:0051782'
    )
)

l1$EvadingApoptosis = unique(
    c(
    'GO:0043069',
    'GO:0043066',
    'GO:0045768',
    'GO:0012501',
    'GO:0043067'
    )
)

l1$LimitlessReplicativePotential = unique(
    c(
        'GO:0001302',
        'GO:0032206',
        'GO:0090398',
        'GO:0032200',
        'GO:0000723',
        'GO:0032204',
        'GO:1900062',
        'GO:2000772'
    )
)
l1$SustainedAngiogenesis = unique(
    c(
        'GO:0045765',
        'GO:0045766',
        'GO:0030949',
        'GO:0001570',
        'GO:0001525',
        'GO:2001212',
        'GO:0008015'
    )
)
l1$TissueInvasionAndMetastasis = unique(
    c(
        'GO:0042060',
        'GO:0007162',
        'GO:0033631',
        'GO:0044331',
        'GO:0016477',
        'GO:0048870',
        'GO:0007155',
        'GO:0001837',
        'GO:0030155',
        'GO:0030030',
        'GO:0030036',
        'GO:0034330',
        'GO:0042330',
        'GO:0007163'
    )
)
l1$GenomeInstabilityAndMutation = unique(
    c(
        'GO:0051276',
        'GO:0045005',
        'GO:0006281',
        'GO:0031570',
        'GO:0045005',
        'GO:0006282'
    )
)
l1$TumorPromotingInflammation = unique(
    c(
        'GO:0002419',
        'GO:0002420',
        'GO:0002857',
        'GO:0002842',
        'GO:0050776',
        'GO:0006954',
        'GO:0002367',
        'GO:0002718',
        'GO:0042060',
        'GO:0061041',
        'GO:0050727',
        'GO:0042533'
    )
)
l1$ReprogrammingEnergyMetabolism = unique(
    c(
        'GO:0006096',
        'GO:0071456',
        'GO:0006006',
        'GO:0046323'
    )
)
l1$EvadingImmuneDetection = unique(
    c(
        'GO:0002837',
        'GO:0002418',
        'GO:0002367',
        'GO:0050776',
        'GO:0006955',
        'GO:0020012',
        'GO:0006897'
    )
)

##### Read GO Term enrichment results
d1 = read.csv("biclusterEnrichment_GOBP.csv",
              header = T,
              row.names = 1)

##### Split GO Terms and store in the list for the loop
l2 = list()
for (cluster in rownames(d1)) {
    l2[[cluster]] = strsplit(as.character(d1[cluster, 2]), ';')[[1]]
}

##### Semantic similarity search for the GO terms in the list
# Build Annotation data
hsGO <- godata('org.Hs.eg.db', ont="BP", keytype = "ENSEMBL")

# Create hallmarks matrix to fill in
hallmarks = matrix(
    ncol = length(names(l1)),
    nrow = length(names(l2)),
    dimnames = list(names(l2), names(l1))
)

# Loop through each cluster and each hallmark to compare their GO term similarities.
for (cluster in names(l2)) {
    if (!(length(l2[[cluster]]) == 0)) {
        for (hallmark in names(l1)) {
            d2 = mgoSim(l2[[cluster]], l1[[hallmark]], semData = hsGO, measure = "Lin", combine = NULL)
            if(length(d2) == 1 & is.na(d2)==T){
                hallmarks[cluster, hallmark] = 0
            }else{
                hallmarks[cluster, hallmark] = max(d2, na.rm = T)
            }
        }
    }
}

##### Write final results into a file.
write.csv(hallmarks, 'Lin_hallmarks.csv')
