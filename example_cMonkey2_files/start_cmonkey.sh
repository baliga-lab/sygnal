#!/bin/bash

../../cmonkey2/bin/cmonkey2.sh --organism hsa --pipeline ../../pipeline.json --config ../../config_pita.ini --rsat_organism Homo_sapiens --rsat_dir ../../commonFiles/data/Homo_sapiens --string ../../commonFiles/geneMania.tsv --synonym_file ../../commonFiles/synonymsThesaurus.csv --verbose --case_sensitive --numclusters 534 ../../exprs/mesothelioma_norm_cv_z.tsv
