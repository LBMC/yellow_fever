# Divergent Patterns of Clonal CD8+ T Cell Differentiation Give Rise to Distinct Lineages of Memory T cells Following Human Viral Infections

This repository contains the script and package used for the analysis of the paper.

## Installation

To get the code and install the nessessary R package you need to run the following commands:

```sh
git clone git@github.com:LBMC/yellow_fever.git
cd yellow_fever
```

If you are on Ubuntu you will also need the following dependency for R:

```sh
sudo add-apt-repository "deb https://stat.ethz.ch/CRAN/bin/linux/ubuntu xenial/"
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
sudo apt-get update -y
sudo apt-get install -y r-base r-base-dev libcurl4-gnutls-dev libxml2-dev libssl-dev
```

Then go in the repository folder and install the following packages in R:

```R
install.packages("devtools")
devtools::install_github("rhondabacher/SCnorm")
devtools::install_github("lme4/lme4")
devtools::install_github("BatzoglouLabSU/SIMLR")
devtools::install_github("rhondabacher/SCnorm", ref = "v1.1.3")
source('http://bioconductor.org/biocLite.R')
biocLite(c("tximport", "readr", "pCMF", "sva", "ComplexHeatmap", "scran"))
install.packages("R2admb")
install.packages("glmmADMB", 
    repos=c("http://glmmadmb.r-forge.r-project.org/repos",
            getOption("repos")),
    type="source")
devtools::install("pkg/")
```

## Analysis

All the scripts used for the analyses can be found in the `src/` folder.

The scripts are numbered from 0 to 12 in execution order with `1_QC.R` the first script for scRNASeq data analysis.
