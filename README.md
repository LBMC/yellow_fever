# install package with ubuntu:

in a terminal:
```sh
sudo add-apt-repository "deb https://stat.ethz.ch/CRAN/bin/linux/ubuntu xenial/"
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
sudo apt-get update -y
sudo apt-get install -y r-base r-base-dev libcurl4-gnutls-dev libxml2-dev libssl-dev
git clone ...
cd scRNAtools/
R
```

in R:
```R
install.packages("devtools")
devtools::install_github("rhondabacher/SCnorm")
devtools::install_github("lme4/lme4")
devtools::install_github("BatzoglouLabSU/SIMLR")
source('http://bioconductor.org/biocLite.R')
biocLite(c("tximport", "readr", "pCMF", "sva", "ComplexHeatmap"))
install.packages("R2admb")
install.packages("glmmADMB", 
    repos=c("http://glmmadmb.r-forge.r-project.org/repos",
            getOption("repos")),
    type="source")
devtools::install()
```
