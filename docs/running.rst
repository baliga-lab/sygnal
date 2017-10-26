Running SYGNAL
==============

Command line parameters
-----------------------

SYGNAL is written in python and takes command line parameters to adjust run behavior:

.. highlight:: none

::

   usage: sygnal.py [-h] [--config CONFIG] [--out OUT] [--res RES]
                 [--cores CORES] [--tmp TMP]
                 cm_rundb

**--config** = SYGNAL configuration file formatted as JSON

**--res** = overall result file name

**--out** = output directory

**--tmp** = temporary file storage directory 

**--cores** = number of cores to be used for 

**cm_rundb** = cMonkey run output SQLite database

Here is an example run training on PITA miRNA-target gene database that follows the directory structure described `previously <https://github.com/cplaisier/sygnal/wiki/Preparing-for-SYGNAL-run#directory-structure-for-running-cmonkey2-and-sygnal>`_:

.. highlight:: none

::

   python sygnal.py --config sygnal_config.json --res postProcessed_example_pita.csv --out output --tmp tmp --cores 19 ../out/cmonkey_run.db


Input files
-----------

The SYGNAL pipeline requires user supplied data to run. This information is used to configure the run.

``sygnal_config.json``
~~~~~~~~~~~~~~~~~~~~~~

Example ``sygnal_config.json``:

.. highlight:: javascript

::

   {
      "max-evalue": 10.0,
      "rho-cut": 0.3,
      "pvalue-cut": 0.05,
      "perc-targets": 0.1,
      "leo-nb-atob": 0.5,
      "mlogp-m-atob": 0.05,

      "run-neo": true,

      "synonyms-file": "../../../commonFiles/synonymsThesaurus.csv.gz",
      "phenotypes-file": "../../../commonFiles/phenotypes_meso.csv",
      "ratios-file": "../../../exprs/<example_cMonkey2_standardized_expression_data>.tsv",
      "all-ratios-file": "../../../exprs_all/<example_all_expression_data>.csv",
      "som-muts-file": "../../../som_muts/<somatic_mutations_matrix.csv",
      "mirna-fasta-file": "miRNA/hsa.mature.fa",
      "rand_pssms_dir": "randPSSMs",
      "promoterSeq": "seqs/promoterSeqs_Homo_sapiens.csv.gz",
      "p3utrSeq": "seqs/p3utrSeqs_Homo_sapiens.csv.gz",
      "gene_conv": "extras/entrez2ucsc.csv",
      "replication-dataset-names": ["<exampleReplication>"],

      "pita": "../../../commonFiles/pita_miRNA_sets_entrez_hsa.json",
      "targetscan": "../../../commonFiles/targetScan_miRNA_sets_entrez_hsa.json",
      "tfbsdb": "../../../commonFiles/tfbsDb_plus_and_minus_5000_entrez.json",

      "meme": {
        "upstream": {
          "nmotifs": 2,
          "motif-widths": [6, 12],
          "revcomp": "True",
          "bgfile": "seqs/bgFile.meme"
        }
      },
      "weeder": {
        "upstream": {
          "bgfile": "HS",
          "size": "small",
          "enriched": "T50",
          "revcomp": "True"
        },
        "3pUTR": {
          "bgfile": "HS3P",
          "size": "small",
          "enriched": "T50",
          "revcomp": "False"
        }
      },
      "tomtom": {
        "upstream": {
          "motif-files": [
          "../../../commonFiles/jasparCoreVertebrata_redundant.json",
          "../../../commonFiles/transfac_2012.1_PSSMs_vertabrate.json",
          "../../../commonFiles/uniprobePSSMsNonRedundant.json",
          "../../../commonFiles/selexPSSMsNonRedundant.json"]
        }
      },
      "tfExpansion": {
        "tfs": "TF/humanTFs_All.csv",
        "tfFamilies": "TF/tfFamilies.csv"
      },
      "mirvestigator": {
        "seedModel": [6,7,8],
        "minor": "True",
        "minor": "True",
        "p5": "True",
        "p3": "True",
        "wobble": "False",
        "wobbleCut": 0.25,
        "species": "hsa"
      },
      "first-principal-comps-result-file": "biclusterFirstPrincComponents.csv"
    }


cMonkey\ :sub:`2` output database
---------------------------------

cmonkey-python writes the results of its computation to an SQLite database. This choice was made, because SQLite is a free, open source and portable data store which is available on many systems and has programming interfaces to a large number of programming languages. Another important aspect is that the entire database is stored in a single file, which can be easily copied, archived and analyzed. The database structure and its function is explained in further detail in the `cMonkey\ :sub:`2` Wiki <https://github.com/baliga-lab/cmonkey2/wiki/Database-schema>`_.

User supplied cMonkey\ :sub:`2` expression matrix
--------------------------------------------------

The cMonkey\ :sub:`2` expression data matrix is usually filtered based on differential expression, most variant genes, or expression above a threshold. The reason for this step is that genes with little variance or low expression are not likely to yield much information.

User supplied expression data (microarray or RNAseq) as a tab separated value (TSV) file that has been properly normalized. The expression data should be standardized (mean subtracted and divided by the standard deviation to turn each expression value into a Z score). The header for the expression data is expected to be without the leading tab.

====== ====== ======== ======
Gene   Cond_1 ...      Cond_N
====== ====== ======== ======
Gene_1 Z[1,1] ...      Z[1,N]
...    ...    ...      ...
Gene_M Z[M,1] Z[M,N-1] Z[M,N]
====== ====== ======== ======

The expression data matrix is specified in the ``sygnal_config.json`` using the ``"ratios-file"`` option.

Complete user supplied expression matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the pre-filtered version of the expression data matrix with all genes. This allows correlation of all putative TF regulators versus bicluster eigengenes even if they were filtered out from the expression data matrix used for cMonkey<sub>2</sub> biclustering.

User supplied expression data (microarray or RNAseq) as a comma separated value (CSV) file that has been properly normalized. The expression data is not standardized. The header for the expression data is not offset .

====== ====== ======== ======
Gene   Cond_1 ...      Cond_N
====== ====== ======== ======
Gene_1 Z[1,1] ...      Z[1,N]
...    ...    ...      ...
Gene_M Z[M,1] Z[M,N-1] Z[M,N]
====== ====== ======== ======

The complete expression data matrix is specified in the ``sygnal_config.json`` using the ``"all-ratios-file"`` option.

Gene ID thesaurus
~~~~~~~~~~~~~~~~~

The gene IDs used in cMonkey<sub>2</sub> and SYGNAL can be different than the pre-computed databases of gene-gene associations, TF- and miRNA-target gene interactions, etc. This is facilitated through a gene ID thesaurus which maps the gene expression matrix gene IDs to alternative gene IDs. This thesaurus is generated from a comma separated values (CSV) file with the following format:

======= ======================
ExprIDs Mappings
======= ======================
Gene_1  Gene_1;Alt_1;...;Alt_N
...     ...;Alt_1;...;Alt_N
Gene_M  Gene_M;Alt_1;...;Alt_N
==============================

The thesaurus file is specified in the ``sygnal_config.json`` using the ``"synonyms-file"`` option.

TF-target gene regulatory interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To discover putative TF regulatory interactions requires a database of TF-target gene interactions. We provide a pre-computed version of TF-target gene interactions from the `TFBS_DB <http://tfbsdb.systemsbiology.net/>`_. For use in SYGNAL the TF-target gene regulatory interactions should be put into the following JSON format:

.. highlight:: none

::

   {"TF_ID_1": ["Gene_1", ..., "Gene_M"], ..., "TF_ID_N": ["Gene_1", ..., "Gene_M"]}


The TF-target gene interaction files are specified in the ``sygnal_config.json`` using the ``"tfbsdb"`` option.

TF family information
~~~~~~~~~~~~~~~~~~~~~

miRNA-target gene regulatory interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To discover putative miRNA regulatory interactions requires a database of miRNA-target gene interactions. We utilize two different pre-computed miRNA-target gene interaction databases:  `PITA <https://genie.weizmann.ac.il/pubs/mir07/mir07_data.html>`_ and `TargetScan <http://www.targetscan.org/vert_71/>`_. For use in SYGNAL the miNRA-target gene regulatory interactions should be put into the following JSON format:

.. highlight:: none

::

   {"miRNA_ID_1": ["Gene_1", ..., "Gene_M"], ..., "miRNA_ID_N": ["Gene_1", ..., "Gene_M"]}

We have started using miRBase `MIMAT` miRNA IDs due to the constant flux in miRNA names. The miRNA-target gene interaction files are specified in the `sygnal_config.json` using the `"pita"` and `"targetscan"` options.

Promoter and 3' UTR sequences for organism of study
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In addition to enrichment with pre-computed databases SYGNAL also searches promoter and 3' UTR sequences for enriched TF or miRNA binding sites. The sequences are stored as comma separated value (CSV) files in the following format:

.. highlight:: none

::

   Gene_1,Seq_1
   ...,...
   Gene_M,Seq_M

The sequence files are specified in the ``sygnal_config.json`` using the ``"promoterSeq"`` and ``"p3utrSeq"`` options.

Background sequence information for motif callers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Typical run of SYGNAL
---------------------

Typical output from running SYGNAL
----------------------------------
