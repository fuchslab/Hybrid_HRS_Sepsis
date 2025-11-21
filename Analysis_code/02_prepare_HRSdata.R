library(dplyr)
# prepare HRS data

load("workspace_similarity_preparation.Rdata")


### single therapies ###

# choose if to use combination therapies
orig_data <- orig_data_single
# complete dataset, only with therapies
myData_onlyTherapies <- cbind("infection_id" = data_rest$infection_id, orig_data_single)
therapies <-  colnames(orig_data_single)


# add patient data
myData <- cbind(select(data_rest, c(infection_id, sex:gfr)), orig_data) 
myData$sex <- as.factor(myData$sex)
myData$focus1 <- as.factor(myData$focus1)

# prepare train data
myData_scale <- myData
# remove observations with missing focus (too much uncertainty)
myData_scale <- myData_scale[which(!is.na(myData_scale$focus1)),]
ids <- myData_scale$infection_id
myData_scale <- myData_scale[,-1]

# only if patient data is included, do the following:
imp.data <- select(myData_scale, sex:gfr)

# median imputation (we have only numeric variables with missing values)
for(m in 1:ncol(imp.data)){
  if(sum(is.na(imp.data[,m])) > 0){
    imp.data[is.na(imp.data[,m]), m] <- median(imp.data[,m], na.rm = TRUE)
  }
}

# scale patient data
train.means <- colMeans(select(imp.data, -c(sex, focus1)))
train.sds <- matrixStats::colSds(as.matrix(select(imp.data, -c(sex, focus1))))

scale.data <- select(imp.data, -c(sex, focus1))
for(s in 1:ncol(scale.data)){
  scale.data[,s] <- (scale.data[,s] - train.means[s])/train.sds[s]
}
myData_scale[colnames(scale.data)] <- scale.data
#--------------------
# one hot encoding
scaled.data.dummy <- fastDummies::dummy_cols(myData_scale, select_columns = c("sex", "focus1"),
                                             remove_selected_columns = TRUE)
factorVars <- c("sex_female", "sex_male", paste0("focus1_", 1:9))
scaled.data.dummy[,factorVars] <- lapply(scaled.data.dummy[,factorVars[which(factorVars %in% colnames(scaled.data.dummy))]], as.factor)
scaled.data.dummy <- relocate(scaled.data.dummy, c(age:gfr, c("sex_female", "sex_male", paste0("focus1_", 1:9))))
myData_scale <- scaled.data.dummy
myData_scale <- cbind("infection_id"= ids, myData_scale)

# complete dataset, only with therapies
therapies <-  colnames(select(myData_scale, Levofloxacin:PenicillinG))
myData_onlyTherapies <- select(data_rest, c(infection_id, all_of(therapies)))

save.image("HRS_data_single.Rdata")
#--------------------------------------------------------------

### combi therapies ###

rm(list = ls())
# prepare HRS data
load("workspace_similarity_preparation.Rdata")

# choose if to use combination therapies
orig_data <- orig_data_combis
# complete dataset, only with therapies
myData_onlyTherapies <- cbind("infection_id" = data_rest$infection_id, orig_data_combis)
therapies <-  colnames(orig_data_combis)

# add patient data
myData <- cbind(select(data_rest, c(infection_id, sex:gfr)), orig_data) 
myData$sex <- as.factor(myData$sex)
myData$focus1 <- as.factor(myData$focus1)

# prepare train data
myData_scale <- myData
# remove observations with missing focus (too much uncertainty)
myData_scale <- myData_scale[which(!is.na(myData_scale$focus1)),]
ids <- myData_scale$infection_id
myData_scale <- myData_scale[,-1]

# only if patient data is included, do the following:
imp.data <- select(myData_scale, sex:gfr)

# median imputation (we have only numeric variables with missing values)
for(m in 1:ncol(imp.data)){
  if(sum(is.na(imp.data[,m])) > 0){
    imp.data[is.na(imp.data[,m]), m] <- median(imp.data[,m], na.rm = TRUE)
  }
}

# scale patient data
train.means <- colMeans(select(imp.data, -c(sex, focus1)))
train.sds <- matrixStats::colSds(as.matrix(select(imp.data, -c(sex, focus1))))

scale.data <- select(imp.data, -c(sex, focus1))
for(s in 1:ncol(scale.data)){
  scale.data[,s] <- (scale.data[,s] - train.means[s])/train.sds[s]
}
myData_scale[colnames(scale.data)] <- scale.data
#--------------------

# one hot encoding
scaled.data.dummy <- fastDummies::dummy_cols(myData_scale, select_columns = c("sex", "focus1"),
                                             remove_selected_columns = TRUE)
factorVars <- c("sex_female", "sex_male", paste0("focus1_", 1:9))
scaled.data.dummy[,factorVars] <- lapply(scaled.data.dummy[,factorVars[which(factorVars %in% colnames(scaled.data.dummy))]], as.factor)
scaled.data.dummy <- relocate(scaled.data.dummy, c(age:gfr, c("sex_female", "sex_male", paste0("focus1_", 1:9))))
myData_scale <- scaled.data.dummy
myData_scale <- cbind("infection_id"= ids, myData_scale)

save.image("HRS_data_combi.Rdata")

