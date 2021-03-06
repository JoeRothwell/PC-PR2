# Vectorised script 2019
# For the script to work, name your omics matrix X_DataMatrixScaled and your metadata Z_meta
# The test data are already missing-free and scaled

# The following packages are required:
library(car)
library(MetabolAnalyze)
library(tidyverse)

# Load the data in the X matrix containing NMR spectra and Z matrix containing the list of explanatory variables of interest
#myPath <- "Documents/PCPR2/" Metabolomics_data <- "X_MetaboMatrix.TXT"
#InterestFactors_data <- "Z_FactorMatrix.TXT"
#Metabo_FilePath = paste(myPath,Metabolomics_data, sep="")
#Factors_FilePath = paste(myPath,InterestFactors_data, sep="")
#X_DataMatrix <- read.delim(Metabo_FilePath, row.names = 1, header = TRUE, sep = "\t")
#Z_InterestFactors <- read.delim(Factors_FilePath, sep = "\t", header = TRUE, row.names = 1)

# Center or scale the data, edit the parameter "pareto" or "unit"
#X_DataMatrixCentered = scale(X_DataMatrix, center = TRUE, scale = FALSE)
#X_DataMatrixScaled = scaling(X_DataMatrixCentered, type = "pareto")

# Load test data
X_DataMatrixScaled <- readRDS("X_testdata.rds")
Z_Meta <- readRDS("Y_testdata.rds")

Z_MetaRowN <- nrow(Z_Meta)
Z_MetaColN <- ncol(Z_Meta)
ColNames   <- names(Z_Meta)

# get number of PCs for threshold: set variability desired to be explained
pct_threshold <- 0.8

# Obtain eigenvectors
pct_threshold <- 0.8 # set variability desired to be explained
X_DataMatrixScaled_t <- t(X_DataMatrixScaled)
symMat <- X_DataMatrixScaled %*% X_DataMatrixScaled_t
eigenData    <- eigen(symMat)
eigenVal     <- eigenData$values
eigenVecMat  <- eigenData$vectors
percents_PCs <- eigenVal/sum(eigenVal)

# Get number of PCs required for threshold (force min to 3)
my_counter_2 <- sum(1 - cumsum(rev(percents_PCs)) <= 0.8)
if(my_counter_2 > 3) pc_n <- my_counter_2 else pc_n <- 3

pc_data_matrix <- eigenVecMat[, 1:pc_n ]

# Convert categorical variables to factors (put them in varlist)
#varlist <- c("sex", "smoking.status")
#Z_Meta <- Z_Meta %>% mutate_at(vars(varlist), as.factor)

DataCol <- Z_MetaColN + 1

# Run a linear model with each eigenvector as the response
TotSumSq <- apply(pc_data_matrix, 2, var) * (Z_MetaRowN - 1)
multifit <- lm(pc_data_matrix ~ ., data = Z_Meta)

# Run type 3 ANOVA on each PC
AnovaTab <- Anova(multifit, type=3, singular.ok = F)
SSP      <- AnovaTab$SSP

# Extract sum of squares for each factor, removing intercept column
# Need to take the diagonal of each factor matrix to get sums of squares
Residuals  <- diag(AnovaTab$SSPE)
RR         <- Residuals/TotSumSq

type3mat0 <- sapply(SSP, diag)[, -1]
type3mat  <- cbind(type3mat0, "SumSqResiduals" = Residuals)
ST_ResidualR2 <- cbind("ST_R2" = 1-RR, "ST_Residuals" = RR)

# Make partial R2 matrix and weighted matrix
partialR2mat <- type3mat[, -DataCol] / (type3mat[, -DataCol] + type3mat[, DataCol])
eigenVal     <- eigenVal[1:pc_n]
weight       <- eigenVal/sum(eigenVal) 

partialR2MatWtProp <- cbind(partialR2mat, ST_ResidualR2[, 1]) * weight
colnames(partialR2MatWtProp) <- NULL
pR2Sums <- colSums(partialR2MatWtProp) * 100

# Plot data
#par(mfrow = c(1,2))
bp <- barplot(pR2Sums, ylab = "Weighted Rpartial2", ylim = c(0, max(pR2Sums) * 1.3), 
          xlab = "", col = "red", las=2, 
          main = paste("n =", Z_MetaRowN, ",", "Y-variables =", Z_MetaColN, ",", "components =", pc_n))
axis(1, at = bp, labels = c(ColNames, "R2"), cex.axis = 0.8, las=2) 
rounded <- round(pR2Sums, 3)
text(bp, pR2Sums, labels = rounded, pos = 3, cex = 0.8)
