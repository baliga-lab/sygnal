## Systems Genetic Network AnaLysis (SYGNAL) pipeline
We developed the SYstems Genetic Network AnaLysis (SYGNAL) pipeline to integrate correlative, causal and mechanistic inference approaches into a unified framework that systematically infers the causal flow of information from mutations to TFs and miRNAs to perturbed gene expression patterns across patients.

The SYGNAL pipeline is designed to be cloned into a completed [cMonkey<sub>2</sub>](https://github.com/baliga-lab/cmonkey2) run directory. It will then run the SYGNAL pipeline using the cMonkey<sub>2</sub> result database and put all the results into 'sygnal/output' directory, which it will create if not already present. There is fairly through checkpointing as the full SYGNAL pipeline can take some time to run.

### Documentation

Documentation is available at https://baliga-lab.github.io/sygnal/

### Dependencies
* MEME (http://meme-suite.org/doc/download.html?man_type=web)
* WEEDER (https://github.com/baliga-lab/weeder_patched)
* R (https://cran.r-project.org/)
* rpy2 (http://rpy2.bitbucket.org/)
* AME (http://bioinformatics.org.au/tools/ame/)

Install with:
```
sudo apt-get install r-base r-base-dev
sudo pip install rpy
```
All other programs will have to installed using their installers and following instructions from those software.

### Required data files
As described the SYGNAL pipeline requires user supplied data for patient expression profiles and patient clinical information. It can be run in a modified form without clinical information. Expression data is usually filtered based on differential expression, most variant genes, or expression above a threshold. The reason for this step is that genes with little variance or low expression are not likely to yield much information.

In addition to user supplied expression data and clinical data the pipeline requires additional information in order to run:
* Gene ID thesaurus
* TF-target gene regulatory interactions in JSON format
* TF family information
* miRNA-target gene regulatory interactions in JSON format
* Full expression data including all genes (for TF correlations)
* Expression data and clinical information for replication studies
* Promoter and 3' UTR sequences for organism of study
* Background sequence information for motif callers

### List of parameters
####TF target-gene database construction				
```
tomtom	(4.11.1)	-dist ed -o tmp/tomtom_out -text -thresh 0.001 -min-overlap 6		
fimo	(4.11.1)	--max-stored-scores 1000000 --verbosity 4 --bgfile bgFile.meme -text --thresh 1e-5 		
```				

####SYGNAL pipeline

```
./sygnal.py [--config <cfg.json>] [--out <output-dir>] [--tmp <tmp-dir>] <cmonkey2-rundb>
```

### Configuration

The SYGNAL pipeline is configured by editing the parameters in the file sygnal_config.json.

##### sygnal_config.json

```
{
}
```

#####Motif searching
```
meme	(4.11.1)	 -bfile bgFile.meme -nostatus -text -time 600 -dna -maxsize 9999999 -evt 1e9 -mod zoops -minw 6 -maxw 12 -nmotifs 2		
weederlauncher	(1.4.2 ISB Patch)	HS small T50 S	For upstream seqeunce searching.	
weederlauncher	(1.4.2 ISB Patch)	HS3P small T50	For 3' UTR seqeunce searching.	
tomtom	(4.11.1)	-dist ed -o tmp/tomtom_out -text -thresh 1 -min-overlap 6 -verbosity 1		
ame	(Original Version)	--method mhg --scoring totalhits --verbose 1 --bgformat 2 --bgfile bgFile.meme		
```


### Output and visualizations
