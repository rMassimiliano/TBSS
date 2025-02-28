# TBSS

This package implements some Tree Based Scan Statistics (TBSS) used in Pharmacoepidemiology. The implementation is R-based and open-source. As such, it can be improved, extended, and re-distributed.  

# Installation
To clone the repository and install the package, you can use the following code:
```
git clone https://gitlab-scm.partners.org/mi475/tbss
https://github.com/rMassimiliano/TBSS
```
Withing the tbss folder  
```
R CMD build TBSS
R CMD INSTALL TBSS
```
Or using `devtools`

```
R> devtools::install_git("https://github.com/rMassimiliano/TBSS", subdir ='TBSS')
```
