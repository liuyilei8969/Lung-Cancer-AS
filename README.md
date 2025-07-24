# Lung Cancer AS

## Project Overview

This repository contains code for alternative splicing analysis and RBP regulation analysis in lung cancer. It can also be applied to similar analyses in bulk RNA-Seq data.

## Key Features

- **Alternative Splicing Analysis**: A precise pipeline for analyzing splicing events in bulk RNA-Seq samples.
- **RBP Regulation Analysis**: Focus on RNA-binding proteins (RBPs) and their regulatory roles in splicing.

## Dependencies

To use this repository, you need to have the following dependencies installed:

- Quantas
- Olego
- Samtools
- R packages(`Boruta`, `dplyr`, `openxlsx`, `ROCR`, `caret`, `ggplot2`, `survival`, `survminer`)
- Python packages (`pandas`, `Bio`, `collections`, `scipy`, `math`, `statsmodels`)

## Usage

### Create a new Conda environment

```bash
conda create -n RNASeq python=3.9
conda activate RNASeq
conda install -c bioconda quantas olego samtools
cd LungCancerAS/
pip install -r requirements.txt
```

### RNA-Seq analysis & quantification  

In this section, you should prepare FASTQ or FASTQ.GZ files of your samples as input. If you would like to reproduce our results, please download the raw data from the National Genomics Data Center (NGDC), as recommended in the "Data and Code Availability" section of our paper.  

- **Scripts:**   
`QC.sh`: Quality control to ensure your samples is suitable for downstream analysis.  
`Quantification.sh`: Import Quantas pipeline to quantify both gene expression and alternative splicing.  
`Analysis.sh`: Identifies differentially expressed genes and differential splicing events between two groups.     

```bash
# Please modify the paths in the scripts to match your data locations  
bash LungCancerAS/RNA-seqAnalysis/QC.sh
bash LungCancerAS/RNA-seqAnalysis/Quatification.sh
bash LungCancerAS/RNA-seqAnalysis/Analysis.sh
```

### Identify AS biomarkers to separate groups in your data      
  
- **Input files required:**   
`information.xlsx`: An Excel file containing two columns. The first column is "Filename" and the second is "Subtype". You may refer to the example provided in `test/information.xlsx`.   
`AS_matrix.txt`: A text file generated in the previous step, representing the PSI matrix. You can refer to `test/AS_matrix.txt`, which contains PSI values for 500 alternative splicing events across 20 samples used in our study.   
`splicing_diff.txt`: A text file also generated in the previous step, showing differential splicing analysis results between tumor and adjacent normal samples. An example is available in `test/splicing_diff.txt`, which includes results of comparing tumor and normal adjacent for 500 events in our study.   

- **Scripts:**   
`FilterDS.R`: We filter differential spliced events in lung cancer, you can adjust your standard and filter on your own.  
`Boruta.R`: Import boruta to figure out candidate features as biomarkers.  
`Overlap.R`: Intersect two sets to identify biomarkers.  

```r
# Please modify the paths in the scripts to match your data locations 
Rscript LungCancerAS/BiomarkerIdentification/FilterDS.R
Rscript LungCancerAS/BiomarkerIdentification/Boruta.R
Rscript LungCancerAS/BiomarkerIdentification/Overlap.R
```

### Validate your AS biomarkers in your own dataset or public datasets

- **Input files required:**   
`information.xlsx`: An Excel file containing two columns. The first column is "Filename" and the second is "Subtype". You may refer to the example provided in `test/information.xlsx`.  
`AS_matrix.txt`: A text file generated in the previous step, representing the PSI matrix. You can refer to `test/AS_matrix.txt`, which contains PSI values for 500 alternative splicing events across 20 samples used in our study.  
`markers.xlsx`:  An Excel file containing event ID of the markers you have found in the previous step.      
  
- **Scripts:**   
`MDS.R`: Multidimensional scaling for checking biomarkers' behavior.  
`SelfCrossValidation.R`: k-fold cross validation with biomarkers and randomly selected AS events to valuate biomarkers' performance.  
`CrossValidationAcrossDatasets.R`: k-fold cross validation with biomarkers in public datasets.  

```r
# Please modify the paths in the scripts to match your data locations 
Rscript LungCancerAS/BiomarkerVerification/MDS.R
Rscript LungCancerAS/BiomarkerVerification/CrossValidationAcrossDatasets.R
Rscript LungCancerAS/BiomarkerVerification/SelfCrossValidation.R
```

### Evaluate prognostic value of AS biomarkers

- **Input files required:**   
`survival.xlsx`: An Excel spreadsheet containing the following columns: sample names, variable values, OS (overall survival status), and OS.time (overall survival time). We obtained these data from https://xenabrowser.net/datapages/.  
  
- **Scripts:**   
`SingleCovariate.R`: Check whether single AS biomarker has prognostic value.  
`MultiCovariate.R`: Check whether a combination of AS biomarkers has prognostic value.  
`ConstructModel.R`: Construct a Cox model with AS biomarkers.  
`SurvivalCurve.R `: Plot two survival curves to show the difference between two groups separated by AS biomarkers.  

```r
# Please modify the paths in the scripts to match your data locations  
Rscript LungCancerAS/Survival/SingleCovariate.R
Rscript LungCancerAS/Survival/MultiCovariate.R
Rscript LungCancerAS/Survival/ConstructModel.R
Rscript LungCancerAS/Survival/SurvivalCurve.R
```

### Motif enrichment analysis of RBP regulation

- **Input files required:**   
`Hs.seq.all.cass.chrom.can.exon.bed`:  A BED file containing event IDs along with the start and end positions of the alternatively spliced exons. You may refer to the example `test/Hs.seq.all.cass.chrom.can.exon.bed`, which includes 500 events.   
`hg19.fa`: Reference genome. The reference genome used in our study was downloaded from the UCSC Genome Browser, with the accession number GCA_000001405.1.            
`motif-RBP.txt`: A text file containing two columns: the first column is motif, and the second column is the corresponding RBP. You may refer to the example file `test/motif-RBP.txt`.   
`up_events.txt`,`dn_events.txt`: Text files summarizing the results of DS analyses. These files can be generated using the provided scripts in the first section.     

- **Scripts:**   
`GrepDSPosition.sh`: Grep DS exons' position and save as a bed file.  
`GrepDSFasta.py`: Generate fasta file of DS exons and flanking 200 bp introns.  
`FindkmerMotif.py`: Count k-mer motifs in the previous-generated fasta files. You can set the minimum and maximum length of motifs.  
`CountRBPMotif.py`: Count each RBP's motifs in the previous-generated fasta files. You can input your own RBP-motif corresponding relationship.  
`MotifEnrich.py`: Enrichment analysis of k-mer motifs in DS events with all casset exons serving as background.  
`RBPEnrich.py`: Enrichment analysis of each RBP's motifs in DS events with all casset exons serving as background.  

```bash
# Please modify the paths in the scripts to match your data locations  
bash LungCancerAS/MotifEnrichment/GrepDSPosition.sh
```
```python
# Please modify the paths in the scripts to match your data locations  
python LungCancerAS/MotifEnrichment/GrepDSFasta.py
python LungCancerAS/MotifEnrichment/FindkmerMotif.py
python LungCancerAS/MotifEnrichment/CountRBPMotif.py
python LungCancerAS/MotifEnrichment/MotifEnrich.py
python LungCancerAS/MotifEnrichment/RBPEnrich.py
```

### Identify RBP regulatory targets

- **Input files required:**   
`Hs.seq.all.cass.chrom.can.exon.bed`:  A BED file containing event IDs along with the start and end positions of the alternatively spliced exons. You may refer to the example `test/Hs.seq.all.cass.chrom.can.exon.bed`, which includes 500 events.   
`ENCODE_DS/`: A directory containing differential splicing results from ENCODE KD RNA-Seq data. Please download the raw data from https://www.encodeproject.org/ and follow the procedure described in the first section to perform the DS analysis.        
`up_events.txt`,`dn_events.txt`, `up_RBP.txt`, `dn_RBP.txt`: Text files summarizing the results of DS and DE analyses. These files can be generated using the provided scripts in the first section.    

- **Scripts:**   
`GenExS.py`: Generate events and splicing factor relationship matrix based on the DS results of RBP-preturbanced cell line, with 0 representing for no relationship, 1 for activation and -1 for repression.  
`GenDExDS`: Generate a smaller matrix than `GenExS.py`, since only DE splicing factor and DS events are taken into consideration.  
`Count.py`: Summarize how many events are activated or repressed by a splicing factor. If you are only concerned about the consistent regulation-relationship, you can filter the ExS matrix.  
`Enrichment.py`: Enrichment analysis of splicing factors' involvement in regulating DS events.  

```python
# Please modify the paths in the scripts to match your data locations  
python LungCancerAS/RBPTargetEnrich/GenExS.py
python LungCancerAS/RBPTargetEnrich/GenDExDS.py
python LungCancerAS/RBPTargetEnrich/Count.py
python LungCancerAS/RBPTargetEnrich/Enrichment.py

```







