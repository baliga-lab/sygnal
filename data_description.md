# Data Description of the SYGNAL output

## Introduction

SYGNAL stores its output in mostly 2 formats:

  1. serialized objects in Python pickle format
  2. csv files

In this document we will describe the structure of each file

## 1. The pickle files

### The cMonkeyWrapper pickle files

### c1.pkl

This is a serialized object of the class cMonkeyWrapper which contains the basic cmonkey2
information

### c1_all.pkl

Based on the c1.pkl object, adds some processing information. Also used to checkpoint the steps
in the SYGNAL computation, by checking the existence of member variables

### c1_postInterim.pkl

Based on the c1_all.pkl object, adds information from the post interim step (in perform_postprocessing()).

### c1_postProc.pkl

Based on the c1_postInterim.pkl object, adds information from the complete postprocessing step
(in perform_postprocessing())

### m2m.pkl

This is a serialized object of the class mirVestigator.

### meme_upstream.pkl

This object is a serialized list of lists of PSSM objects. Each entry in the list of the top level
list is from a separate MEME run.

### miRvestigatorResults.pkl

A serialized dictionary of dictionaries in the format

{ <motif_name>: { 'miRNA': <miRNAname(s)>, 'model': <k-mer>, 'mature_seq_ids': <individual miRNA names> } }


### pita_3pUTR.pkl

Serialized weeder results, in the format of a dictionary.

{ <run id>: [list of PSSM] }

### postProcessesd.pkl

Generated by parallelized invokations of the function post_process().
A serialized list of dictionaries
[ { 'k': <cluster>, 'k.rows': <number of rows>, 'k.cols': <number of conditions> } ]

### targetscan_3pUTR.pkl

Results of the set enrichment computation (targetscan). A list of results from running the function cluster_hypergeo()
on every bicluster.

### tfbs_db.pkl

Results of the set enrichment computation (TFBS DB). A list of results from running the function cluster_hypergeo()
on every bicluster.

### upstreamJasparTransfacComparison.pkl

Serialized TOMTOM results.

### weeder_3pUTR.pkl

Serialized results of Weeder on 3p UTR

### weeder_upstream.pkl

Serialized results of Weeder on upstream

### The miRNA directory (miRvestigator generated)

Is processed by the miRNA step to generate other output

## 2. The csv files

### biclusterEigengenes.csv

Results of the getEigengene.R script

### biclusterEnrichment_GOBP.csv

Results of the enrichment.R script


### biclusterFirstPrincComponents.csv

First principal components, exported from the main SYGNAL script

### biclusterVarianceExplained.csv

Result from getEigengene.R script


### causalitySummary.csv

Results from NEO, written in write_neo_summary() step in sygnal.py


### mut_reg_fold_changes_*.csv

Results from runNEO.R

### mut_reg_t_test_p_values_*.csv

Results from runNEO.R

### replicationPvalues_mesoTCGA.csv

TODO

### residualPermutedPvalues_permAll.csv

Result from permutedResidualPvalues_permAll_mc.R

### upstreamComparison_jaspar_transfac.csv

Result written by sygnal.py in the __tomtom_upstream_motifs() function

### upstreamMotifPermutedPValues.csv

Result written by sygnal.py in the __get_permuted_pvalues_for_upstream_meme_motifs() function

### The causality_* directory

This contains the results of the tool NEO (network edge orienting)

### postProcessed_example_pita.csv

result file written from sygnal.py and specified through the --res switch