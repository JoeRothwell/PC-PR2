# The Principal Component Partial R-squared method (PC-PR2)


#### Background

The PC-PR2 is a statistical method, developed by Fages *et al*. (2014)[^1], used to investigate sources of variability in metabolomics or other omics data. It combines features of principal component analysis and multivariable linear regression analyses. The input is a complete X-matrix of omics data and a corresponding set of descriptive Y-data (subject metadata). The output is the proportion of variation in the omics data attributed to each Y-variable.

The original version of the PCPR-2 code is stored in this repository, as well as the latest script, which is in continuous development.

Test data, consisting of a sample of a transcriptomics dataset, is also included. This consists of five descriptive variables for the 124 subjects (two categorical, three numeric) and 3000 corresponding transcriptomics intensities.

[^1]:Fages et al (2014) Investigating sources of variability in metabolomic data in the EPIC study: the Principal Component Partial R-square (PC-PR2) method. *Metabolomics* 10(6): 1074-1083, DOI: 10.1007/s11306-014-0647-9