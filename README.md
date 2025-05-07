# TBSS

This package implements some Tree-Based Scan Statistics (TBSS) used in Pharmacoepidemiology. It is R-based and open-source, so that it can be improved, extended, and redistributed.  

# Installation
To clone the repository and install the package, you can use the following code:
```
git clone https://github.com/rMassimiliano/TBSS
```
Within the tbss folder  
```
R CMD build TBSS
R CMD INSTALL TBSS
```
Or using `devtools`

```
R> devtools::install_git("https://github.com/rMassimiliano/TBSS", subdir ='TBSS')
```
