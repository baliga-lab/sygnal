Preparing for SYGNAL run
========================

Overview of SYGNAL method
-------------------------

The SYGNAL pipeline is designed to be cloned into a completed cMonkey\ :sub:`2` run directory. It will then run the SYGNAL pipeline using the cMonkey\ :sub:`2` result database and put all the results into an output directory. SYGNAL performs all the post-hoc analysis of the cMonkey\ :sub:`2` biclusters that is required for subsequent regulatory network construction.

SYGNAL is designed for use on mammalian species and the cMonkey\ :sub:`2` runs have been tailored to use the more computationally efficient set-enrichment methods for bicluster training on TF and miRNA regulatory inputs. This is used in place of the MEME or WEEDER based promoter or 3' UTR sequence searching which very computationally expensive if done iteratively. The SYGNAL pipeline was developed to streamline post-hoc running of these methods. The primary input for the SYGNAL pipeline is the results from a cMonkey\ :sub:`2` run which is comprised of a set of biclusters typically stored as a SQLite database (https://github.com/baliga-lab/cmonkey2/wiki/Database-schema). SYGNAL then conducts the following analyses:

  1. Retrieve cMonkey\ :sub:`2` biclusters from database
  2. Discover putative TF and miRNA regulators
  3. Determine significance of bicluster co-expression
  4. Statistical tests comparing study phenotypes versus bicluster eigengene (e.g. survival analysis, T-test, or correlation)
  5. Replicate in addition studies if provided
  6. Functional enrichment of GO biological processes
  7. Semantic similarity of enriched GO biological processes to hallmarks of disease
  8. Causal analysis linking mutations to regulator to bicluster eigengene

Currently the SYGNAL pipeline requires some code manipulation to turn each of these analyses on or off. We are working towards making each of these optional and will allow users to turn them on or off through the configuration file. The causality analysis is a working example of this proposed functionality as it can be turned on or off through a flag in the ``sygnal_config.json`` configuration file.

Directory structure for running cMonkey2 and SYGNAL
---------------------------------------------------

.. highlight:: none

::

   +-- cMonkey2
   |  +-- (clone of https://github.com/baliga-lab/cmonkey2.git)
   +-- commonFiles
   |  +-- (common files such as synonymThesaurus.csv)
   +-- config_pita.ini (configuration file for PITA trained run of cMonkey2)
   +-- config_targetscan.ini (configuration file for TargetScan trained run of cMonkey2)
   +-- config_tfbs_db.ini (configuration file for TFBS_DB trained run of cMonkey2)
   +-- exprs
   |  +-- (user defined expression matrix formatted and standardized for cMonkey2, typically filtered to less than 8,000 genes).tsv
   +-- exprsAll
   |  +-- (all genes expression matrix not standardized).csv
   +-- miRNA_seq
   |  +-- (all miRNAs expression matrix not standardized).csv
   +-- pipeline.json (scoring configuration file for cMonkey2)
   +-- runs
   |  +-- run.pita
   |  |  +-- (link all cMonkey2 files here using 'ln -s ../../cMonkey2/* .')
   |  |  +-- cache
   |  |  +-- out
   |  |  |  +-- cmonkey_run.db (cMonkey2 SQLite database file used as input for SYGNAL)
   |  |  |  +-- (etc.)
   |  |  +-- start_cmonkey.sh (script to start cMonkey run)
   |  |  +-- sygnal
   |  |  |  +-- (link all sygnal files here using 'ln -s ../../../sygnal/* .')
   |  |  |  +-- output (directory holding output for SYGNAL)
   |  |  |  |  +-- biclusterEigengenes.csv (matrix of bicluster eigenegenes from WGCNA)
   |  |  |  |  +-- biclusterEnrichment_GOBP.csv (GO biological process enrichments per bicluster)
   |  |  |  |  +-- biclusterFirstPrincComponents.csv (matrix of bicluster first principal components)
   |  |  |  |  +-- biclusterVarianceExplained.csv (variance explained by first principal component for each bicluster)
   |  |  |  |  +-- c1_all.pkl (2nd SYGNAL checkpoint)
   |  |  |  |  +-- c1.pkl (initial SYGNAL checkpoint)
   |  |  |  |  +-- c1_postInterim.pkl (3rd SYGNAL checkpoint)
   |  |  |  |  +-- c1_postProc.pkl (final SYGNAL checkpoint)
   |  |  |  |  +-- causal (directory holding causality analyses)
   |  |  |  |  +-- causalitySummary.csv (summary of causality analyses)
   |  |  |  |  +-- cluster.members.conditions.txt (list of bicluster condition membership)
   |  |  |  |  +-- cluster.members.genes.txt (list of bicluster gene membership)
   |  |  |  |  +-- jiangConrath_hallmarks.csv (semantic similarity for each hallmark of disease)
   |  |  |  |  +-- m2m.pkl (miRvestigator checkpoint)
   |  |  |  |  +-- meme_upstream.pkl (checkpoint for MEME promoter motif discovery)
   |  |  |  |  +-- miRNA (directory for miRNA pre-computed files)
   |  |  |  |  +-- miRvestigatorResults.pkl (results from miRvestigator runs)
   |  |  |  |  +-- pita_3pUTR.pkl (enrichment results from PITA analyses)
   |  |  |  |  +-- postProcessed_meso_pita.csv (final result file)
   |  |  |  |  +-- postProcessed.pkl (checkpoint for phenotype analysis)
   |  |  |  |  +-- replicationPvalues_mesoTCGA.csv (results form replication study)
   |  |  |  |  +-- residualPermutedPvalues_permAll.csv (permuted p-values for co-expression significance)
   |  |  |  |  +-- targetscan_3pUTR.pkl (enrichment results from TargetScan analyses)
   |  |  |  |  +-- tfbs_db.pkl (enrichment results from TFBS_DB analyses)
   |  |  |  |  +-- upstreamComparison_jaspar_transfac.csv (results of comparing motifs versus motif databases)
   |  |  |  |  +-- upstreamJasparTransfacComparison.pkl (checkpoint of comparing motifs versus motif databases)
   |  |  |  |  +-- upstreamMotifPermutedPValues.csv (significance of discovered promoter motifs)
   |  |  |  |  +-- weeder_3pUTR.pkl (checkpoint for WEEDER 3' UTR motif discovery)
   |  |  |  |  +-- weeder_upstream.pkl (checkpoint for WEEDER promoter motif discovery)
   |  |  |  +-- start_sygnal.sh (script to start SYGNAL run)
   |  |  |  +-- tmp (directory to store temporary files)
   |  +-- run.targetscan
   |  |  +-- (repeat of run.pita with paths modified to feed off of run.targetscan)
   |  +-- run.tfbs_db
   |  |  +-- (repeat of run.pita with paths modified to feed off of run.tfbs_db)
   +-- somMuts
   |  +-- (mutations matrix with values as [0,1]).csv
   +-- sygnal
   |  +-- (clone of https://github.com/cplaisier/sygnal.git)


Running cMonkey\ :sub:`2`
-------------------------
The training configuration for cMonkey\ :sub:`2` should include co-expression, GeneMania gene-gene interaction network (http://genemania.org/data/current/), and enrichment of either TF or miRNA target genes using the set-enrichment module (Reiss et al., 2015). cMonkey\ :sub:`2` should be run three times. The first run used the TF-target gene interaction database as input to the set-enrichment module to discover TF mediated regulation. The second and third runs used PITA (Kertesz et al., 2007) and TargetScan (Friedman et al., 2009) as input to the set-enrichment module to discover miRNA mediated regulation.

``start_cMonkey.sh``
--------------------

(with cmonkey\ :sub:`2` installed through its pip package on pypi.org)

.. highlight:: none

::

   cmonkey2 --organism hsa --pipeline ../../pipeline.json --config ../../config_pita.ini --rsat_organism Homo_sapiens --rsat_dir ../../commonFiles/data/Homo_sapiens --string ../../commonFiles/data/STRING/string_conv.tsv --synonym_file ../../commonFiles/synonymsThesaurus.csv --verbose --case_sensitive --numclusters 534 ../../exprs/mesothelioma_norm_cv_z.tsv

``pipeline.json``
-----------------

.. highlight:: javascript

::

   {
     "row-scoring": {
        "id": "combiner",
        "function": { "module": "cmonkey.scoring", "class": "ScoringFunctionCombiner" },
        "args": {
          "functions": [
            { "id": "Rows",
              "function": { "module": "cmonkey.microarray", "class": "RowScoringFunction" }
            },
            { "id": "Networks",
              "function": { "module": "cmonkey.network", "class": "ScoringFunction" }
            },
            { "id": "SetEnrichment",
              "function": { "module": "cmonkey.set_enrichment", "class": "ScoringFunction" }
            }
          ]
        }
      },
      "column-scoring": { "id": "Columns",
                          "function": { "module": "cmonkey.scoring",
                                        "class": "ColumnScoringFunction"} }
    }

``pita_config.json``
--------------------

.. highlight:: none

::

   [General]
   normalize_ratios = True
   num_iterations = 2000

   [SetEnrichment]
   schedule = 1,7
   scaling_rvec=seq(1e-5, 0.5, length=num_iterations/2)
   map_to_ratio_genes = True
   set_types = pita

   [SetEnrichment-pita]
   set_file = ../../commonFiles/pita_miRNA_sets_entrez_hsa.json
   weight = 1.0


Overlaying SYGNAL on alternative clustering methods
---------------------------------------------------

It is possible to modify SYGNAL to sit on top of other clustering methods. The focus of such efforts should be to modify the `cMonkeyWrapper.py <https://github.com/cplaisier/sygnal/blob/master/cMonkeyWrapper.py>`_ such that it pulls the information from the new source of cluster information.
