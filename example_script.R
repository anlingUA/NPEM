### NPEM Example Script ###

## Load required packages
library(np)
library(outliers)

## Load in example data
X <- read.csv("./geneExp.csv",row.names=1)
M <- read.csv("./OTUabnd.csv",row.names=1)
Y <- read.csv("./Diagnosis.csv",row.names=1)

## Put data in correct format
M <- apply(M,c(1,2),function(s) log(s+1))
Y$Diagnosis <- as.factor(Y$Diagnosis)

## Run NPEM with Univariate Mediator
UV.example <- NPEM(X,M,Y,method="UV") 
UV.example$mediation.p


## Run NPEM with Bivariate Mediator with test after iteration
BVS.example <- NPEM(X,M,Y,method="BVS")
BVS.example$mediation.p
