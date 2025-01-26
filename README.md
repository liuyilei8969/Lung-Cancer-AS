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
`Analysis.sh`: Figure out differential expressed genes and differential spliced events between two sample-groups.

```bash
# please change the path as your sample path in the analysis scripts
bash LungCancerAS/RNA-seqAnalysis/QC.sh
bash LungCancerAS/RNA-seqAnalysis/Quatification.sh
bash LungCancerAS/RNA-seqAnalysis/Analysis.sh
```

### Figure out AS biomarkers to separate two groups in your data

```r
# please change the path as your sample path in the analysis scripts
Rscript LungCancerAS/BiomarkerIdentification/FilterDS.R
Rscript LungCancerAS/BiomarkerIdentification/Boruta.R
Rscript LungCancerAS/BiomarkerIdentification/Overlap.R
```

### Verify your AS biomarkers in your own data or public datasets

```r
# please change the path as your sample path in the analysis scripts
Rscript LungCancerAS/BiomarkerVerification/MDS.R
Rscript LungCancerAS/BiomarkerVerification/CrossValidationAcrossDatasets.R
Rscript LungCancerAS/BiomarkerVerification/SelfCrossValidation.R
```

### Check whether your AS biomarkers have prognostic value

```r
# please change the path as your sample path in the analysis scripts
Rscript LungCancerAS/Survival/SingleCovariate.R
Rscript LungCancerAS/Survival/MultiCovariate.R
Rscript LungCancerAS/Survival/ConstructModel.R
Rscript LungCancerAS/Survival/SurvivalCurve.R
```

### Part of RBP regulation analysis, here you are intereted in the motifs on the spliced RNA sequences

```bash
# please change the path as your sample path in the analysis scripts
bash LungCancerAS/MotifEnrichment/GrepDSPosition.sh
```

```python
# please change the path as your sample path in the analysis scripts
python LungCancerAS/MotifEnrichment/GrepDSFasta.py
python LungCancerAS/MotifEnrichment/FindkmerMotif.py
python LungCancerAS/MotifEnrichment/CountRBPMotif.py
python LungCancerAS/MotifEnrichment/MotifEnrich.py
python LungCancerAS/MotifEnrichment/RBPEnrich.py
```

### Part of RBP regulation analysis, here you are intereted in the RBPs' targets

```python
# please change the path as your sample path in the analysis scripts
python LungCancerAS/RBPTargetEnrich/GenExS.py
python LungCancerAS/RBPTargetEnrich/GenDExDS.py
python LungCancerAS/RBPTargetEnrich/Count.py
python LungCancerAS/RBPTargetEnrich/Enrichment.py

```







