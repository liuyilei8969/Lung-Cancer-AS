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

`QC.sh`: Quality control to ensure your samples is suitable for downstream analysis.  
`Quantification.sh`: Import Quantas pipeline to quantify both gene expression and alternative splicing.  
`Analysis.sh`: Figure out differential expressed genes and differential spliced events between two groups.  

```bash
# please change the path as your sample path in the analysis scripts
bash LungCancerAS/RNA-seqAnalysis/QC.sh
bash LungCancerAS/RNA-seqAnalysis/Quatification.sh
bash LungCancerAS/RNA-seqAnalysis/Analysis.sh
```

### Figure out AS biomarkers to separate two groups in your data

`FilterDS.R`: We filter differential spliced events in lung cancer, you can adjust your standard and filter on your own.  
`Boruta.R`: Import boruta to figure out candidate features as biomarkers.  
`Overlap.R`: Intersect two sets to identify biomarkers.  

```r
# please change the path as your sample path in the analysis scripts
Rscript LungCancerAS/BiomarkerIdentification/FilterDS.R
Rscript LungCancerAS/BiomarkerIdentification/Boruta.R
Rscript LungCancerAS/BiomarkerIdentification/Overlap.R
```

### Verify your AS biomarkers in your own data or public datasets

`MDS.R`: Multidimensional scaling for checking biomarkers' behavior.  
`SelfCrossValidation.R`: k-fold cross validation with biomarkers and randomly selected AS events to valuate biomarkers' performance.  
`CrossValidationAcrossDatasets.R`: k-fold cross validation with biomarkers in public datasets.  

```r
# please change the path as your sample path in the analysis scripts
Rscript LungCancerAS/BiomarkerVerification/MDS.R
Rscript LungCancerAS/BiomarkerVerification/CrossValidationAcrossDatasets.R
Rscript LungCancerAS/BiomarkerVerification/SelfCrossValidation.R
```

### Check whether your AS biomarkers have prognostic value

`SingleCovariate.R`: Check whether single AS biomarker has prognostic value.  
`MultiCovariate.R`: Check whether a combination of AS biomarkers has prognostic value.  
`ConstructModel.R`: Construct a Cox model with AS biomarkers.  
`SurvivalCurve.R `: Plot two survival curves to show the difference between two groups separated by AS biomarkers.  

```r
# please change the path as your sample path in the analysis scripts
Rscript LungCancerAS/Survival/SingleCovariate.R
Rscript LungCancerAS/Survival/MultiCovariate.R
Rscript LungCancerAS/Survival/ConstructModel.R
Rscript LungCancerAS/Survival/SurvivalCurve.R
```

### Part of RBP regulation analysis, here you are intereted in the motifs on the spliced RNA sequences

`GrepDSPosition.sh`: Grep DS exons' position and save as a bed file.  

```bash
# please change the path as your sample path in the analysis scripts
bash LungCancerAS/MotifEnrichment/GrepDSPosition.sh
```

`GrepDSFasta.py`: Generate fasta file of DS exons and flanking 200 bp introns.  
`FindkmerMotif.py`: Count k-mer motifs in the previous-generated fasta files. You can set the minimum and maximum length of motifs.  
`CountRBPMotif.py`: Count each RBP's motifs in the previous-generated fasta files. You can input your own RBP-motif corresponding relationship.  
`MotifEnrich.py`: Enrichment analysis of k-mer motifs in DS events with all casset exons serving as background.  
`RBPEnrich.py`: Enrichment analysis of each RBP's motifs in DS events with all casset exons serving as background.  

```python
# please change the path as your sample path in the analysis scripts
python LungCancerAS/MotifEnrichment/GrepDSFasta.py
python LungCancerAS/MotifEnrichment/FindkmerMotif.py
python LungCancerAS/MotifEnrichment/CountRBPMotif.py
python LungCancerAS/MotifEnrichment/MotifEnrich.py
python LungCancerAS/MotifEnrichment/RBPEnrich.py
```

### Part of RBP regulation analysis, here you are intereted in the RBPs' targets

`GenExS.py`: Generate events and splicing factor relationship matrix based on the DS results of RBP-preturbanced cell line, with 0 representing for no relationship, 1 for activation and -1 for repression.  
`GenDExDS`: Generate a smaller matrix than `GenExS.py`, since only DE splicing factor and DS events are taken into consideration.  
`Count.py`: Summarize how many events are activated or repressed by a splicing factor. If you are only concerned about the consistent regulation-relationship, you can filter the ExS matrix.  
`Enrichment.py`: Enrichment analysis of splicing factors' involvement in regulating DS events.  

```python
# please change the path as your sample path in the analysis scripts
python LungCancerAS/RBPTargetEnrich/GenExS.py
python LungCancerAS/RBPTargetEnrich/GenDExDS.py
python LungCancerAS/RBPTargetEnrich/Count.py
python LungCancerAS/RBPTargetEnrich/Enrichment.py

```

## Test Data

Small sample data for testing purposes can be found in the `test/` directory.







