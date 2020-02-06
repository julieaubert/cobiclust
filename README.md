# cobiclust

**cobiclust** is an R package dedicated to the biclustering of over-dispersed count data such as those produced by amplicon-based sequencing for example. 

## How to install cobiclust?

**cobiclust** is available on [CRAN](https://cran.r-project.org/web/packages/cobiclust/).

**cobiclust** needs the _cluster_ CRAN R package, so check that it is installed on your computer.

`is.installed <- is.element("cluster", rownames(installed.packages()))`

`if (is.installed == FALSE) install.packages("cluster")`


### For the last stable version, use the CRAN version

`install.packages("cobiclust")`

### Installation from GitHub

To install the **cobiclust** package within R from GitHub, open a R session and:

* install `devtools` with `install.packages("devtools")` (if not installed yet)
* load the `devtools` R package with `library(devtools)`
* run install_github("julieaubert/cobiclust")

## About cobiclust

The cobiclust package has been developped at UMR518 MIA-Paris (INRAE, AgroParisTech, Université Paris-Saclay) by J. Aubert (<julie.aubert@agroparistech.fr>). 
