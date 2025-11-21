# functions for RS
library(dplyr)
library(caret)
library(rlist)
library(mboost)


# HRS function

# prepare data accordingly before:
# input_data: main dataset; should contain a variable infection_id (first column), 
# clinical patient data (variables with core values, vital signs, lab values; already scaled, one hot encoded and without missing values) for one timepoint each, 
# binary therapy variables with AST results
# therapies: variable names (contained in input_data) that represent the therapies
# patientSim: symmetric matrix with patient similarities
# therapySim: symmetric matrix with therapy similarities
# simChoice: "threshold" (for a similarity threshold) or "kN" (for k most similar neighbors)
# k: number of k most similar neighbors
# thresholdP: patient similarity threshold
# thresholdT: therapy similarity threshold


# functions for memory_based
memory1 <- function(input_data, therapies, seed = 123, patientSim, simChoice, k, thresholdP){   
  Eff <- as.matrix(select(input_data, all_of(therapies)))
  for(i in 1:nrow(input_data)){
    # calculate effectiveness for every patient without value for t
    # choose similar patients
    # simChoice should be either threshold or kN
    if(simChoice == "kN"){                # based on k nearest neighbors 
      p_sim <- names(sort(patientSim[i,], decreasing = TRUE)[2:k])
    } else if (simChoice == "threshold"){ # based on threshold
      p_sim <- names(patientSim[i, ][patientSim[i, ] >= thresholdP])
    } 
    # data set with only similar patients
    data_sim <- input_data[input_data$infection_id %in% p_sim,]
    for(t in therapies){
      if(is.na(input_data[i,t])){
        # keep only patients who really got therapy of interest
        data_t <- na.omit(select(data_sim, c(infection_id, t)))
        if(nrow(data_t) > 0){ # only if there is at least one similar patient with given therapy t
          sum1 <- sum2 <- 0
          # effectiveness formula
          # for all similar patients
          for(sp in 1:nrow(data_t)){
            sim_sp <- patientSim[as.character(input_data$infection_id[i]), as.character(data_t[sp,]$infection_id)] 
            sum1 <- sum1 + (as.numeric(data_t[sp,t]) * sim_sp) 
            sum2 <- sum2 + sim_sp
          }
          Eff[i,t] <- sum1/sum2
        }
      }
    }
  }
  return(Eff)
}

memory2 <- function(input_data, therapies, seed = 123, patientSim, therapySim, simChoice, k, thresholdP, thresholdT){
  Eff <- as.matrix(select(input_data, all_of(therapies)))
  for(i in 1:nrow(input_data)){
    # calculate effectiveness for every patient without value for t
    # choose similar patients
    # simChoice should be either threshold or kN
    if(simChoice == "kN"){                # based on k nearest neighbors 
      p_sim <- names(sort(patientSim[i,], decreasing = TRUE)[2:k])
    } else if (simChoice == "threshold"){ # based on threshold
      p_sim <- names(patientSim[i, ][patientSim[i, ] >= thresholdP])
    } 
    # data set with only similar patients
    data_sim <- input_data[input_data$infection_id %in% p_sim,]
    for(t in therapies){
      if(is.na(input_data[i,t])){
        # take similar therapies into account
        t_sim <- names(therapySim[t, ][therapySim[t, ] >= thresholdT])
        # keep only therapies which are available in data
        t_sim <- t_sim[t_sim %in% colnames(input_data)] 
        # keep only patients who really got therapy of interest
        data_t <- select(data_sim, c(infection_id, all_of(t_sim)))
        # remove rows with only NAs
        data_t = data_t[rowSums(is.na(data_t)) != (ncol(data_t)-1), ]
        # remove columns with only NAs
        data_t = data_t[,colSums(is.na(data_t)) != nrow(data_t)]
        
        if(nrow(data_t) > 0){ # only if there is at least one similar patient with given therapy t
          sum1 <- sum2 <- 0
          # effectiveness formula
          # for all similar patients
          for(st in colnames(data_t)[-1]){
            data_st <- na.omit(data_t[,c("infection_id", st)])
            for(sp in 1:nrow(data_st)){
              sim_sp <- patientSim[as.character(input_data$infection_id[i]), as.character(data_st[sp,]$infection_id)] 
              sum1 <- sum1 + (as.numeric(data_st[sp,st]) * sim_sp * therapySim[t, st])
              sum2 <- sum2 + (sim_sp * therapySim[t, st])
            }
          }
          Eff[i,t] <- sum1/sum2
        }
      }
    }
  }
  return(Eff)
}

# HRS: "memory1" or "memory2" or "model" (choice of HRS)
# matrixUpdate: "iterative" (for collective update) or "direct" (for immediate update)
# initialValues: if HRS = "model", which initial values to use: predicted values from "memory1" or "memory2"
# stop_iter: stop after stop_iter iterations
# saveModel: shoudl all final models be saved? TRUE/FALSE

hrs <- function(HRS, input_data, therapies, seed = 123, patientSim, therapySim, simChoice, k, thresholdP, thresholdT,
                matrixUpdate, initialValues, stop_iter, saveModel){
  # chose which HRS to use
  if(HRS == "memory1"){ #Ia
    eff <- memory1(input_data, therapies, seed, patientSim, simChoice, k, thresholdP)
    return(eff)
  } else if(HRS == "memory2"){ #Ib
    eff <- memory2(input_data, therapies, seed, patientSim, therapySim, simChoice, k, thresholdP, thresholdT)
    return(eff)
  } 
  else if(HRS == "model"){ #IIa or IIb (depending on data input)
    if(initialValues == "memory1"){
      input_data_CF <- memory1(input_data, therapies, seed, patientSim, simChoice, k, thresholdP)
    } else if(initialValues == "memory2"){
      input_data_CF <- memory2(input_data, therapies, seed, patientSim, therapySim, simChoice, k, thresholdP, thresholdT)
    }
    
    input_data_orig <- input_data
    # prediction for all therapies, one after the other
    preds <- covs <- betas <- list()
    data_class_new <- list(input_data_CF)
    convergeCheck <- FALSE
    # loop
    i <- 1
    while(convergeCheck == FALSE & i <= stop_iter){ # stop at stop_iter iterations
      data <- data_class_new[[i]]
      data <- cbind(data, input_data_orig[!colnames(input_data_orig) %in% colnames(data)])
      data <- select(data, -infection_id)
      
      preds_i <- covs_i <- betas_i <- models_i <- list()
      if(matrixUpdate == "iterative"){
        data_class_new[[i+1]] <- data_class_new[[i]]
      }
      
      for(t in 1:length(therapies)){
        set.seed(123 + t)
        na_ids <- which(is.na(input_data_orig[, therapies[t]]))
        train <- as.data.frame(data[-na_ids, ])    # independent variables for train
        data.pred <- as.data.frame(data[na_ids, ])      # independent variables for test
        #------
        # remove factor covariates with only one factor level
        rm <- c()
        for(ii in 1:ncol(train)){
          if(is.factor(train[,ii])){
            if(any(table(train[,ii]) < 1))
              rm <- c(rm, colnames(train)[ii])
          }
        }
        train <- select(train, -all_of(rm))
        train[,t] <- as.factor(train[,t])
        
        #------------
        # classification model
        # mboost model
        boost_ctrl <- boost_control(mstop = 500, # initial number of boosting iterations. Default: 100
                                    nu = 0.05, # step length. Default: 0.1
                                    trace = FALSE) # print status information? Default: FALSE
        
        # calculate weights for inbalanced classes
        n0 <- sum(train[,therapies[t]] == 0)
        n1 <- sum(train[,therapies[t]] == 1)
        w0 <- nrow(train) / (2 * n0)
        w1 <- nrow(train) / (2 * n1)
        weights_t <- ifelse(train[,therapies[t]] == 1, w1, w0)
        
        # model
        set.seed(seed)
        boost_model <- mboost::glmboost(as.formula(paste0("`", therapies[t],"`",  " ~ .")),
                                        data = train,
                                        family = Binomial(),
                                        control = boost_ctrl,
                                        center = TRUE, 
                                        weights = weights_t)
        
        # optimization
        # 10-fold cv
        my.cv <- mboost::cv(model.weights(boost_model), type = "kfold",
                            B = 10, prob = 0.8, strata = NULL)
        risk <- cvrisk(boost_model, folds = my.cv)
        
        # final model
        mboost_final <- boost_model[mstop(risk)]
        # save chosen covariate names
        if(is.null(names(mboost_final$coef()))){
          covs_i[[t]] <- NA
        } else {
          covs_i[[t]] <- names(mboost_final$coef())
        }
        names(covs_i)[t] <- therapies[t]
        # save beta coefficients
        betas_i[[t]] <- as.numeric(unlist(mboost_final$coef()))
        names(betas_i)[t] <- therapies[t]
        # save model objects
        models_i[[t]] <- mboost_final
        names(models_i)[t] <- therapies[t]
        
        # final model
        class_final_model <- mboost_final
        # prediction
        pred_na <- predict(object = class_final_model,
                           newdata = data.pred,
                           type = 'response')
        
        #------------------
        # save predictions
        preds_i[[t]] <- pred_na
        names(preds_i)[t] <- therapies[t]
        if(matrixUpdate == "iterative"){
          data_class_new[[i+1]][na_ids, therapies[t]] <- pred_na
        } else {
          # update therapies of interest immediately within current dataset
          data[na_ids, therapies[t]] <- pred_na
        }
      }
      preds[[i]] <- preds_i
      covs[[i]] <- covs_i
      betas[[i]] <- betas_i
      
      # check if convergence is reached for all therapies
      if(i > 1){
        if(identical(preds[[i]], preds[[i-1]])){
          convergeCheck <- TRUE
        }
      }
      
      # in case of direct matrixUpadte, directly update the current data matrix with new predicted results
      if(matrixUpdate == "direct"){
        data_class_new[[i+1]] <- data
      }
      # set iteration to next loop
      i <- i + 1
    }
    if(saveModel == TRUE){
      results <- list(data_loops = data_class_new, preds_loops = preds, covs_loops = covs, betas_loops = betas, final_models = models_i)
    } else {
      results <- list(data_loops = data_class_new, preds_loops = preds, covs_loops = covs, betas_loops = betas)
    }
    return(results)
  }
}

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------

# evaluation function

# fold: number of CV folds


# row-wise validation
validate_classic <- function(HRS, input_data, therapies, seed = 123, patientSim, therapySim, k, thresholdP, simChoice, thresholdT,
                             stop_iter, initialValues, matrixUpdate, fold){ 

  set.seed(seed)
  fold_ids <- createFolds(y = 1:nrow(input_data), k = fold)
  results <- list()
  for(f in 1:fold){
    test_ids <- fold_ids[[f]]
    train.data <- input_data
    train.data[test_ids, therapies] <- NA

    test.data <- input_data[test_ids, therapies] 
    # HRS
    eff <- hrs(HRS = HRS,
               input_data = train.data, 
               therapies = therapies,
               seed = seed,
               patientSim = patientSim,
               therapySim = therapySim,
               k = k,
               thresholdP = thresholdP,
               simChoice = simChoice,
               thresholdT = thresholdT,
               stop_iter = stop_iter,
               initialValues = initialValues,
               matrixUpdate = matrixUpdate, 
               saveModel = FALSE)
    
    # check performance
    test.positions <- c(as.matrix(test.data))
    test.positions <- which(!is.na(test.positions))
    
    if(HRS == "memory1" | HRS == "memory2"){
      prediction <- eff[test_ids, therapies]
      prediction <- c(unlist(prediction))[test.positions]
    } else if(HRS == "model"){
      prediction <- eff$data_loops[[length(eff$data_loops)]][test_ids, therapies]
      prediction <- c(as.matrix(prediction))[test.positions]
    }
    y_true <- c(unlist(test.data))[test.positions]
    compare <- as.data.frame(cbind("y_pred" = prediction, "y_true" = y_true))
    compare$absDiff <- abs(compare$y_pred - compare$y_true)
    
    
    results[[f]] <- list(test_entries = test_ids,
                         eff = eff,
                         compare = compare)
    print(f)
  }
  
  # one model on whole data
  final_eff <- hrs(HRS = HRS,
                   input_data = input_data,
                   therapies = therapies,
                   seed = seed,
                   patientSim = patientSim,
                   therapySim = therapySim,
                   k = k,
                   thresholdP = thresholdP,
                   simChoice = simChoice,
                   thresholdT = thresholdT,
                   stop_iter = stop_iter,
                   initialValues = initialValues,
                   matrixUpdate = matrixUpdate, 
                   saveModel = TRUE)
  
  df_results <- results %>% list.select(compare) %>% bind_rows()
  return(list(results = results, df_results = df_results, final_eff = final_eff))
}

#-----------------------------
# cell-wise validation
validate_RS <- function(HRS, input_data, therapies, seed = 123, patientSim, therapySim, k, thresholdP, simChoice, thresholdT,
                        stop_iter, initialValues, matrixUpdate, fold){ 
  set.seed(seed)
  # fold_ids (ratings which will be set to NA in the training data)
  full_entries <- which(!is.na(as.matrix(input_data[, therapies]))) 
  fold_entries <- createFolds(y = full_entries, k = fold)
  results <- list()
  for(f in 1:fold){
    # set test entries to NA
    test_entries <- full_entries[fold_entries[[f]]]
    train.th <- as.matrix(input_data[, therapies])
    test.data <- train.th[test_entries]
    train.th[test_entries] <- NA
    input_data.f <- input_data
    input_data.f[, therapies] <- train.th
    input_data.f$infection_id <- input_data$infection_id
    # HRS
    eff <- hrs(HRS = HRS,
               input_data = input_data.f,
               therapies = therapies,
               seed = seed,
               patientSim = patientSim,
               therapySim = therapySim,
               k = k,
               thresholdP = thresholdP,
               simChoice = simChoice,
               thresholdT = thresholdT,
               stop_iter = stop_iter,
               initialValues = initialValues,
               matrixUpdate = matrixUpdate, 
               saveModel = FALSE)
    # check performance
    if(HRS == "memory1" | HRS == "memory2"){
      prediction <- as.matrix(eff)[test_entries]
    } else if(HRS == "model"){
      prediction <- as.numeric(as.matrix(eff$data_loops[[length(eff$data_loops)]])[test_entries])
    }
    y_true <- test.data
    compare <- as.data.frame(cbind("y_pred" = prediction, "y_true" = y_true))
    compare$absDiff <- abs(compare$y_pred - compare$y_true)
    
    results[[f]] <- list(test_entries = test_entries,
                         eff = eff,
                         compare = compare)
    print(f)
  }
  
  # one model on whole data
  final_eff <- hrs(HRS = HRS,
                   input_data = input_data,
                   therapies = therapies,
                   seed = seed,
                   patientSim = patientSim,
                   therapySim = therapySim,
                   k = k,
                   thresholdP = thresholdP,
                   simChoice = simChoice,
                   thresholdT = thresholdT,
                   stop_iter = stop_iter,
                   initialValues = initialValues,
                   matrixUpdate = matrixUpdate, 
                   saveModel = TRUE)
  
  
  df_results <- results %>% list.select(compare) %>% bind_rows()
  return(list(results = results, df_results = df_results, final_eff = final_eff))
}



