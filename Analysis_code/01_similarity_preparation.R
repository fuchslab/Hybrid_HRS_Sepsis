# load packages
library(tidyverse)
library(rlist)
library(caret)
library(plyr)
library(dplyr)
library(cluster)
library(StatMatch)


# load data
load("dataPreparation.Rdata")

# patient data
patientData$focus1 <- as.factor(patientData$focus1)
patientData$sex <- as.factor(patientData$sex)

### Patient similarity
intVars <- colnames(patientData)
# covariate weights
cov_weights2 <- c(1,3,5,3,2,3,4,3,5,4,4,3,4,4,5,3,4,4) 

names(cov_weights2) <- colnames(patientData)
cov_weights2
# gowers distance
gower_d <- daisy(patientData,
                 metric = "gower",
                 stand = TRUE,
                 weights = cov_weights2[colnames(patientData)])

patientData2 <- patientData
str(patientData)

# integer in numeric 
patientData2$focus1 <- as.character(patientData2$focus1)
patientData2$sex <- as.character(patientData2$sex)
patientData2 <- as.data.frame(patientData2)
str(patientData2)

gd <- gower.dist(data.x = patientData2, KR.corr = FALSE, var.weights = cov_weights2[colnames(patientData2)])
rm(patientData2)

# similarity matrix
gower_mat <- as.matrix(gower_d)
simP <- 1-gower_mat
rownames(simP) <- colnames(simP) <- data_rest$infection_id
summary(c(simP))


#-----------------------
### Therapy similarity
# use Pathogen-Therapy matrix: coverage

# define simT with matches/mismatches
simT <-  matrix(ncol = ncol(coverage), nrow = ncol(coverage))

for(i in 1:nrow(testmatrix)){ #therapies
  for(k in 1:nrow(testmatrix)){
    vec <- c()
    for(j in 1:ncol(testmatrix)){ #pathogens
      if(testmatrix[i,j] == testmatrix[k,j]){
        vec[j] <- 1
      } else if(testmatrix[i,j] == 0 & testmatrix[k,j] == 1 |
                testmatrix[i,j] == 1 & testmatrix[k,j] == 0 |
                testmatrix[i,j] == 0 & testmatrix[k,j] == 2 |
                testmatrix[i,j] == 2 & testmatrix[k,j] == 0 ){
        vec[j] <- 0
      }  else {
        vec[j] <- 0.5 # 1,2 or 2,1
      }
    }  
    simT[i, k] <- mean(vec)
  }
}

colnames(simT) <- rownames(simT) <- colnames(coverage)
summary(c(simT))


save.image("workspace_similarity_preparation.Rdata")


