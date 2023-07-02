## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(foreach)
library(doParallel)
library(DataExplorer)
library(dplyr)
library(pROC)
library(caret)
library(parallel)
library(DMwR)
library(ggplot2)
library(GGally)
library(corrplot)
library(ROSE)


NUM_CORES = 6
options(mc.cores = NUM_CORES)

sink("./Outputs/CurrentOutput.txt")

#Set random seed
set.seed(42)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#For satimage dataset
#
# data <- read.csv("../Datasets/satimage_csv.csv")
# 
# OUTPUT_VAR = "class"
# data[[OUTPUT_VAR]][data[[OUTPUT_VAR]] %in% c(1, 2, 3, 5, 7)] <- "MAJ"
# data[[OUTPUT_VAR]][data[[OUTPUT_VAR]] == 4] <- "MIN"
# ind.cualit <- c(which(names(data) == "class"))

#For adult dataset


data = read.csv("../Datasets/adult.csv", nrows=5000)
data <- data[!apply(data == "?", 1, any), ]

#Removed varuables for simplicity education.num, relationship, fnlwgt, capital.gain and capital.loss
data$occupation<-NULL
data$fnlwgt<-NULL
data$educational.num<-NULL
data$relationship<-NULL

 ind.cualit <- c(which(names(data) == "workclass"),which(names(data)=="education"),which(names(data)=="income"),which(names(data)=="marital.status"):which(names(data)=="gender"),which(names(data)=="native.country"))
 OUTPUT_VAR = "income"
 data = na.omit(data)
 data[[OUTPUT_VAR]] = factor(data[[OUTPUT_VAR]], levels = c("<=50K", ">50K"), labels = c("Low", "High"))
 


#For covid dataset
# 
# data = read.csv("../Datasets/data_tfg.csv")
# ind.cualit <- c(which(names(data) == "SITUACION"),which(names(data)=="SEXO"),which(names(data)=="DM"):which(names(data)=="DC"))
# OUTPUT_VAR = "SITUACION"


# For correlation plot of some characteristics
# minoritary <- data[data$SITUACION == "FALLECIDO", ]
# majoritary <- data[data$SITUACION == "CURADO", ]
# M<-cor(minoritary[, !(1:ncol(minoritary) %in% ind.cualit)])
# corrplot(M, method = "circle")
# 
# M<-cor(majoritary[, !(1:ncol(majoritary) %in% ind.cualit)])
# corrplot(M, method = "circle")



 for(i in ind.cualit){
  
  data[,i] <- as.factor(data[, i])
  

 }


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
barplot(prop.table(table(data[[OUTPUT_VAR]])),
        col = rainbow(3),
        ylim = c(0, 1.01),
        main = "Class Distribution")


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


train.index <- createDataPartition(data[[OUTPUT_VAR]], p = 0.85, list = FALSE)
data.train <- data[ train.index,]
data.test  <- data[-train.index,]



## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

output_lev <- levels(data.train[[OUTPUT_VAR]])

lev_nrow <- c(nrow(data.train[data.train[[OUTPUT_VAR]] == output_lev[1],]), nrow(data.train[data.train[[OUTPUT_VAR]] == output_lev[2],]))

if (lev_nrow[1] < lev_nrow[2]){
  OUTPUT_MIN = output_lev[1]
  OUTPUT_MAJ = output_lev[2]
  nmin = lev_nrow[1]
  nmax = lev_nrow[2]
}else{
  OUTPUT_MIN = output_lev[2]
  OUTPUT_MAJ = output_lev[1]
  nmin = lev_nrow[2]
  nmaj = lev_nrow[1]
}
print(sprintf("%s : %d", OUTPUT_MIN, nmin))
print(sprintf("%s : %d", OUTPUT_MAJ, nmaj))



## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


kappa_mine <- function(TP, TN, FP, FN, lObs) {
    po <- (TP + TN) / lObs
    pe <- ((TP + FP) * (TP + FN) + (FP + TN) * (FN + TN)) / lObs^2
    (po - pe) / (1 - pe)
  }

sensitivity_mine <- function(TP, TN, FP, FN) {
    TP / (TP + FN)
}

specificity_mine <- function(TP, TN, FP, FN) {
    TN / (TN + FP)
  }

# metrics <- function(data, lev = levels(as.factor(data$obs)), model = NULL) {
#   #Calculate TP, TN, FP, FN
#   pred = data[,"pred"]
#   obs = data[,"obs"]
#   
#   TP <- sum(pred == OUTPUT_MIN & obs == OUTPUT_MIN)
#   TN <- sum(pred == OUTPUT_MAJ & obs == OUTPUT_MAJ)
#   FP <- sum(pred == OUTPUT_MIN & obs == OUTPUT_MAJ)
#   FN <- sum(pred == OUTPUT_MAJ & obs == OUTPUT_MIN)
#   
#   # Calculate Cohen's kappa and BalAcc
#   KAPPA <- kappa_mine(TP, TN, FP, FN, length(obs))
#   BAL_ACC <- (sensitivity_mine(TP, TN, FP, FN) + specificity_mine(TP, TN, FP, FN)) / 2
#   
#   # Return metrics as a named vector
#   c(KAPPA = KAPPA, BAL_ACC = BAL_ACC)
# }

metrics <- function(data, lev = levels(as.factor(data$obs)), model = NULL){
    sensi = sensitivity(data[, "pred"],data[, "obs"],positive=OUTPUT_MIN,negative=OUTPUT_MAJ)
    met <- c(
    KAPPA = psych::cohen.kappa(cbind(data[, "obs"],data[, "pred"]))$kappa,
    SENS = sensi,
    BAL_ACC = ( sensi + specificity(data[, "pred"], data[, "obs"],positive=OUTPUT_MIN,negative=OUTPUT_MAJ))/2
    )

  return(met)
}



metrics_all <- function(data, lev = levels(as.factor(data$obs)), model = NULL){
    met <- c(
    ACCURACY = MLmetrics::Accuracy(data[, "pred"], data[, "obs"]),
    SENS = sensitivity(data[, "pred"],data[, "obs"],positive=OUTPUT_MIN,negative=OUTPUT_MAJ),
    SPEC = specificity(data[, "pred"], data[, "obs"],positive=OUTPUT_MIN,negative=OUTPUT_MAJ),
    PPV = posPredValue(data[, "pred"], data[, "obs"],positive=OUTPUT_MIN,negative=OUTPUT_MAJ),
    NPV = negPredValue(data[, "pred"], data[, "obs"],positive=OUTPUT_MIN,negative=OUTPUT_MAJ),
    KAPPA = psych::cohen.kappa(cbind(data[, "obs"],data[, "pred"]))$kappa,
    BAL_ACC = (sensitivity(data[, "pred"],data[, "obs"],positive=OUTPUT_MIN,negative=OUTPUT_MAJ) + specificity(data[, "pred"], data[, "obs"],positive=OUTPUT_MIN,negative=OUTPUT_MAJ))/2,
    F1 = MLmetrics::F1_Score(data[,"pred"],data[,"obs"]),
    AUC = MLmetrics::AUC(data[,"prob"],data[,"obs.prob"]),
    PR = MLmetrics::PRAUC(data[,"prob"],data[,"obs.prob"]),
    MCC = mltools::mcc(data[,"pred"],data[,"obs"]),
    GEOM = sqrt(sensitivity(data[, "pred"],data[, "obs"],positive=OUTPUT_MIN,negative=OUTPUT_MAJ)*specificity(data[, "pred"], data[, "obs"],positive=OUTPUT_MIN,negative=OUTPUT_MAJ))
  )
  return(met)
}

number_metrics = 12
#We add time as a measure we will compute later
metric_names <- c("ACCURACY", "SENS", "SPEC", "PPV", "NPV", "KAPPA", "BAL_ACC", "F1", "AUC", "PR", "MCC", "GEOM", "TIME")


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
calculate_np <- function( nmin, nmaj,  elbow_prob=0.75) 
{
  np <- round(nmin*elbow_prob)
  return(np)
}



## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
calculate_p <- function( np, elbow_prob=0.75, alpha_p= .01) 
{
  p <- ceiling(log(alpha_p)/(log(1-1/np)*np))
  return(p)
}



## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  
calculate_b <- function(np, nmin, nmaj, elbow_prob=0.75, alpha_b= .01) 
{

  b <- ceiling(log(alpha_b)/(log(1-1/nmin)*np))
  return(b)
}





## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



mt <- function(b, n) { ceiling((b-n) / 3) }




## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
get_function_vector <- function(b,function_training){
  function_vector <- c() 
  for(i in 1:b){
    function_vector <- append(function_vector, function_training) 
  } 
  

  return(function_vector)
}



## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


seed_algorithms <- c(

#RFOREST
 function(df.train, metrics, OUTPUT, metric_optimize) {

    
    tC <- trainControl(
      summaryFunction = metrics,
      allowParallel = TRUE,
      classProbs = TRUE
    )
    
    
    method <- "ranger"
    maximize <- T
    
    cl <- makeCluster(NUM_CORES, type="FORK" )
    clusterExport(cl, c("OUTPUT_MIN", "OUTPUT_MAJ"))
    registerDoParallel(cl)
    
    
    rf <- train(
      as.formula(sprintf("%s ~.", OUTPUT)),
      data = df.train,
      method = method,
      metric = metric_optimize,
      maximize = maximize,
      importance = "impurity",
      trControl = tC
    )
    
    stopCluster(cl)
    registerDoSEQ()
    
    return(rf)
  },
#RLOG
 function(df.train, metrics, OUTPUT, metric_optimize) {


    tC <- trainControl(
      summaryFunction = metrics,
      allowParallel = TRUE,
      classProbs = TRUE
    )
    
    method <- "glmnet"
    maximize <- T
    
       
    cl <- makeCluster(NUM_CORES, type="FORK")
    clusterExport(cl, c("OUTPUT_MIN", "OUTPUT_MAJ"))
    registerDoParallel(cl)
    
    
    rlog <- train(as.formula(sprintf("%s ~.", OUTPUT)),
      data = df.train,
      method = "glmnet",
      family = 'binomial',
      metric = metric_optimize,
      maximize = T,
      tuneGrid = expand.grid(
        alpha = seq(0, 1, by = 0.2),
        lambda = seq(0.0001, 1, length = 100)
      ),
      trControl = tC
    )
    
        
    stopCluster(cl)
    registerDoSEQ()

    
    return(rlog)
  },
#SVM
function(df.train, metrics, OUTPUT, metric_optimize) {
  
  
  
  tC <- trainControl(
    summaryFunction = metrics,
    allowParallel = TRUE,
    classProbs = TRUE
  )
  
  # method <- "svmLinear"
  
   cl <- makeCluster(NUM_CORES, type = "FORK")
   clusterExport(cl, c("OUTPUT_MIN", "OUTPUT_MAJ"))
   registerDoParallel(cl)
  # 
  # svm <- train(
  #   as.formula(sprintf("%s ~.", OUTPUT)),
  #   data = df.train,
  #   method = method,
  #   metric = metric_optimize,
  #   maximize = TRUE,
  #   trControl = tC
  # )
  
    svm <- train(as.formula(sprintf("%s ~.", OUTPUT)),
      data = df.train,
      method = "glmnet",
      family = 'binomial',
      metric = metric_optimize,
      maximize = T,
      tuneGrid = expand.grid(
        alpha = seq(0, 1, by = 0.2),
        lambda = seq(0.0001, 1, length = 100)
      ),
      trControl = tC
    )
  
  stopCluster(cl)
  registerDoSEQ()
  
  return(svm)
},

#GBM
 function(df.train, metrics, OUTPUT, metric_optimize) {

    tC <- trainControl(
      summaryFunction = metrics,
      allowParallel = TRUE,
      classProbs = TRUE,
      verboseIter = FALSE
      
    )
    
    method <- "gbm"
    maximize <- T
    
            
    cl <- makeCluster(NUM_CORES, type="FORK")
    clusterExport(cl, c("OUTPUT_MIN", "OUTPUT_MAJ"))
    registerDoParallel(cl)
    
    gbm <- train(as.formula(sprintf("%s ~.", OUTPUT)),
                 data = df.train,
                 method = method,
                 metric = metric_optimize,
                 maximize = maximize,
                 verbose = FALSE,
                 trControl = tC
    )
    
    stopCluster(cl)
    registerDoSEQ()
    
    return(gbm)
  })


alg_names <-c("RANGER", "RLOG", "SVM", "GBM")



## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
imbalancedFold <- function(data, n_folds, target, minority_class) {
  n_samples <- nrow(data)
  n_majority_total <- n_samples - sum(data[[target]] == minority_class)
  n_minority_total <- sum(data[[target]] == minority_class)
  
  n_minority_per_fold <- ceiling(n_minority_total / n_folds)
  n_majority_per_fold <- ceiling(n_samples / n_folds - n_minority_per_fold)

  fold_indices <- list()
  
  for (i in 1:n_folds) {
    
    #We choose the samples of the majoritary class
    majority_indices <- which(data[[target]] != minority_class)
    used_majority_indices <- unlist(fold_indices)
    available_majority_indices <- setdiff(majority_indices, used_majority_indices)
    
    n_available_majority <- length(available_majority_indices)
    n_majority_this_fold <- min(n_majority_per_fold, n_available_majority)
    selected_majority_indices <- sample(available_majority_indices, size = n_majority_this_fold)
    
    #We take the samples of the minoritary class
    
    minority_indices <- which(data[[target]] == minority_class)
    used_minority_indices <- setdiff(used_majority_indices, available_majority_indices)
    available_minority_indices <- setdiff(minority_indices, used_minority_indices)
    n_available_minority <- length(available_minority_indices)
    n_minority_this_fold <- min(n_minority_per_fold, n_available_minority)
    selected_minority_indices <- sample(available_minority_indices, size = n_minority_this_fold)
    
    #We collect and combine both indexes
    selected_indices <- c(selected_majority_indices, selected_minority_indices)

    fold_indices[[i]] <- selected_indices
  }
  
  return(fold_indices)
}



## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

folds <-imbalancedFold(data.train, 5,OUTPUT_VAR, OUTPUT_MIN)



## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

default_q = 0.75


prediction <- function(conj.model, x, q = default_q){ 
  pred <- data.frame(matrix(nrow=nrow(x),ncol=0))
  for(model in conj.model) pred <- cbind(pred, predict(model,x))
  nElems = ncol(pred)
  nElems <- ncol(pred)
  counts <- rowSums(pred == OUTPUT_MAJ)
  pred <- counts
  ifelse(is.na(pred) | pred<q*nElems, OUTPUT_MIN, OUTPUT_MAJ)
}

prediction.fold <- function(ensemble, x, q = default_q){
  pred <- as.data.frame(lapply(ensemble, function(e) prediction(e,x)))
  nElems = ncol(pred)
  nElems <- ncol(pred)
  counts <- rowSums(pred == OUTPUT_MAJ)
  pred <- counts
  ifelse(is.na(pred) | pred<q*nElems, OUTPUT_MIN, OUTPUT_MAJ)
}

prediction.add <- function(old_ensemble, new_element, prediction_current = FALSE, test.set = FALSE){
    if(prediction_current == FALSE){
      prediction_current = as.numeric(unlist(predict(new_element, test.set)) == OUTPUT_MAJ)
    }
    if(length(old_ensemble$ensemble)!=0){
      prev_count =  old_ensemble$model_count
      l <- list("ensemble"= c(old_ensemble$ensemble,list(new_element)), "model_count"= prev_count + prediction_current)
      return(l)
    }
    else{
      l <- list("ensemble"= list(new_element), "model_count"=prediction_current)
      return(l)
    }
    
  
}

prediction.predicted <- function(predicted, q = default_q) {

    nElems <- length(predicted$ensemble) 
    pred <- ifelse(predicted$model_count  < q*nElems, OUTPUT_MIN, OUTPUT_MAJ)
    return(pred)
}




## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
prediction.prob <- function(ensemble, x, q = default_q){
  pred <- data.frame(matrix(nrow=nrow(x),ncol=0))
  for(model in ensemble) pred <- cbind(pred, predict(model,x,type="prob")[[OUTPUT_MAJ]])
  pred<- rowSums(pred)
  pred
}

prediction.fold.prob <- function(ensemble, x, q =default_q){
   #Place all predictions for a sample in each row of a data set
  pred <- as.data.frame(lapply(ensemble, function(e) prediction.prob(e,x)))
  pred <- rowSums(pred)/sum(unlist(lapply(ensemble,length)))
  pred
}



prop.maj = 0.55
metric_optimize = 'BAL_ACC'



## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
train_IPIP <- function( prop.maj, OUTPUT, min_str, maj_str, train.set, test.set, 
                        configuration, prediction, metrics, b, np,  p, mt, TRAIN_METRIC = "BAL_ACC", SELECT_METRIC = "BAL_ACC"){
  
  nmin = sum(train.set[[OUTPUT]] == min_str)
  nmaj = sum(train.set[[OUTPUT]] == maj_str)
  
  

  majSubSize <-  round(np*prop.maj/(1-prop.maj))
  minSubSize <-  np

  #Includes in each position the elements of each partition from 1 to p
  dfs <- list()
  
  minoritary <- subset(train.set, train.set[[OUTPUT]] == min_str)
  majoritary<- subset(train.set, train.set[[OUTPUT]] == maj_str)
  
  for(k in 1:p){
    id.minoritary <- sample(x = 1:nmin, size = minSubSize) #Index for minoritary class for each subset
    id.majoritary <- sample(x= 1:nmaj, size = majSubSize) #Indexes for majoritary class for each subset

    dfs[[k]] <- rbind(minoritary[id.minoritary,],majoritary[id.majoritary,])
  }
  
  
  E <- list() #Final model (ensemble of ensembles)


  for(k in 1:p){
    Ek <- list() # k-esim ensemble
    i <- 0 #Counter for number of tries of enlarging the ensemble
    metric.ensemble  = 0 #Variable for storing the metric of each ensemble accumulated.
    #Balanced partition
    
    df <- dfs[[k]]
    model_i = 1
    
    #We select the elements for the training
    majoritary <- which(df[[OUTPUT]] == maj_str)
    minoritary <- which(df[[OUTPUT]] == min_str)


    while(length(Ek$ensemble)<=b && i<mt(b,length(Ek$ensemble))){
      #We use replace TRUE to make some randomness over the seed learners
      ind.train <- c(
          sample(majoritary, size = majSubSize, replace = TRUE),
          sample(minoritary, size = minSubSize, replace = TRUE) 
      )
      
    
      #We train with the models in configuration with a sequential approach
      
      model <- configuration[[length(Ek$ensemble)+1]](df[ind.train,], metrics, OUTPUT, TRAIN_METRIC)
      
      
      if (length(Ek$ensemble)==0){
          u <- -Inf;
          names(u) <- SELECT_METRIC;
          metrics.ensemble <- u
      }
      
      new_Ek <- prediction.add(Ek, model, test.set = test.set[colnames(test.set)!=OUTPUT])
      pred = prediction.predicted(new_Ek)
      metrics.ensemble.new <- metrics(data.frame(
        obs = test.set[[OUTPUT]],
        pred= as.factor(pred)
      ))
      #We check if the new model changes the result. If it does not we start trying again
      print(metrics.ensemble.new[SELECT_METRIC])
      print(metrics.ensemble[SELECT_METRIC])
      if(metrics.ensemble.new[SELECT_METRIC] <= metrics.ensemble[SELECT_METRIC]){ 
        i <- i+1
      } else{ 
        #If the ensemble tries to enlarge again, we restart the enlarging posibilities.
        Ek <- new_Ek
        metrics.ensemble <- metrics.ensemble.new
        i <- 0
      }
    } # End of the k-esim ensemble building

    E[[length(E)+1]] <- Ek$ensemble
  
  }
return(E);
}




## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
powerset = function(s, test.set){
    len = length(s)
    l = c()
    vector(length=2^len);l[[1]]=numeric()
    counter = 1L
    for(x in 1L:len){
        prediction_current =as.numeric(unlist(predict(s[x], test.set)) == OUTPUT_MAJ)
        prev_count = 0
        for(subset in 1L:counter){
            counter=counter+1
            if(subset == 1){
              l[[counter]] = list("ensemble"= s[x], "model_count"= prediction_current)
            }
            else{
              prev_count = l[[subset]]$model_count
              l[[counter]] = list("ensemble"= c(l[[subset]]$ensemble,s[x]), "model_count"= prev_count + prediction_current)
            }
        }
    }
    return(tail(l, length(l)-1))
}


train_IPIP_exhaustive <- function( prop.maj, OUTPUT, min_str, maj_str, train.set, test.set, configuration, prediction, metrics, np,  p, TRAIN_METRIC = "BAL_ACC", SELECT_METRIC = "BAL_ACC"){


  
  nmin = sum(train.set[[OUTPUT]] == min_str)
  nmaj = sum(train.set[[OUTPUT]] == maj_str)
  
  
  
  majSubSize <-  round(np*prop.maj/(1-prop.maj))
  minSubSize <-  np

  #Incluye en cada posicion los valores de los elementos de dicha particion, de 1 a p
  dfs <- list()
  
  minoritary <- subset(train.set, train.set[[OUTPUT]] == min_str)
  majoritary <- subset(train.set, train.set[[OUTPUT]] == maj_str)

  for(k in 1:p){
    id.minoritary <- sample(x = 1:nmin, size = minSubSize) #Index for minoritary class for each subset
    id.majoritary <- sample(x= 1:nmaj, size = majSubSize) #Indexes for majoritary class for each subset

    dfs[[k]] <- rbind(minoritary[id.minoritary,],majoritary[id.majoritary,])
  }
  

  
  E <- list() #Final model (ensemble of ensembles)


  for(k in 1:p){
      
    cat(sprintf("Ensemble number %d of %d\n", k,p))
    Ek <- list() # k-esim ensemble
    init = Sys.time()
  
    #Balanced partition
    
    df <- dfs[[k]]
    majoritary <- which(df[[OUTPUT]] == maj_str)
    minoritary <- which(df[[OUTPUT]] == min_str)
    ind.train <- c(
          sample(majoritary, size = majSubSize, replace = TRUE),
          sample(minoritary, size = minSubSize, replace = TRUE) 
    )
    train.test <- df[ind.train,]
      Ek <- list()
      for(conf in 1:length(configuration)){

        Ek[[length(Ek)+1]] <- configuration[[conf]](train.test, metrics, OUTPUT, TRAIN_METRIC)

      }
      
      all_sets <- powerset(Ek, as.data.frame(test.set[, -which(names(test.set) == OUTPUT)]))
      
      ##If we want to force bigger sets before smaller ones
      #lengths <- sapply(all_sets, function(x) length(x))
      # Reorder the list based on the lengths
      # all_sets <- all_sets[order(lengths, decreasing = TRUE)]
      

      
      max_metric <- function(model){
        pred = prediction.predicted(model)
        met <- metrics(data.frame(
              obs = test.set[[OUTPUT]],
              pred = as.factor(pred)
          ))[SELECT_METRIC]
          return(met)
      }
      
      
      
      #Run the max metric check in parallel mode
      max_metric_list <- lapply(all_sets, max_metric)
      #print(cbind(lapply(all_sets, function(x) length(x)/2), max_metric_list))
      # find the maximum chosen metric value and corresponding set
      
      max_metric_value <- max(unlist(max_metric_list))
      max_set <- all_sets[[which.max(unlist(max_metric_list))]]$ensemble
      cat(sprintf("Max metric %s of ensemble is %f with length %d\n", SELECT_METRIC, max_metric_value, length(max_set)))
      E[[length(E)+1]] <- max_set
  
  }
return(E);
}




## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
metricPreparedFrame <- function(model, setToTest)
      data.frame(
          obs = setToTest[[OUTPUT_VAR]],
          pred = predict(model,setToTest[,names(setToTest) != OUTPUT_VAR]),
          prob= prediction.fold.prob(list(list(model)), setToTest[,names(setToTest) != OUTPUT_VAR]),
          obs.prob = as.numeric(ifelse(setToTest[[OUTPUT_VAR]] == OUTPUT_MAJ, 1, 0))
        )



## ----warning=FALSE, cache=TRUE-------------------------------------------------------------------------------------------------------------------------------------------------------------
mean_metrics.seed <- list()
mean_metrics.test.seed <- list()
mean_time.seed <- list()


cat("######## Training seed algorithms naive approach ######## \n")
for(alg in 1:length(seed_algorithms)){
    mean_time.seed[[alg]] = 0
    
    metrics.final.seed <- list()
    metrics.test.seed <- list()
    cat(sprintf("Seed algorithm: %d\n", alg))
    for (i in 1:length(folds)) {
        cat(sprintf("Fold %d out of %d\n", i, length(folds)))
        cat("------------------------------------\n")
        
    
        
        
        train.set <- data.train[unlist(folds[i]),]
        test.set <- data.train[-unlist(folds[i]),]
    
        
        start_time <- Sys.time()
        model_seed <- seed_algorithms[[alg]](train.set, metrics, OUTPUT_VAR,metric_optimize)
        end_time <- Sys.time()
        mean_time.seed[[alg]]= mean_time.seed[[alg]] +
        end_time - start_time
        print(end_time - start_time)
        
        
        metrics.final.seed <- append( metrics.final.seed, metrics_all(metricPreparedFrame(model_seed, test.set)))
        metrics.test.seed <- append( metrics.test.seed, metrics_all(metricPreparedFrame(model_seed, data.test)))
        
        cat("------------------------------------\n")
    }
    
  mean_metrics.seed <- append(mean_metrics.seed, 
        apply(matrix(unlist(metrics.final.seed), ncol= number_metrics, byrow=T), 2, mean))
  
  mean_metrics.test.seed <- append(mean_metrics.test.seed, 
        apply(matrix(unlist(metrics.test.seed), ncol= number_metrics, byrow=T), 2, mean))
}





matrix_mean.seed <- matrix(mean_metrics.seed, nrow = length(seed_algorithms), ncol = number_metrics, byrow = TRUE)
matrix_mean.seed<- cbind(matrix_mean.seed, lapply(mean_time.seed, function(x){as.numeric(x, units="secs")}/length(folds)))
colnames(matrix_mean.seed) <-metric_names
rownames(matrix_mean.seed) <- alg_names


matrix_mean.test.seed <- matrix(mean_metrics.test.seed, nrow = length(seed_algorithms), ncol = number_metrics, byrow = TRUE)
matrix_mean.test.seed <- cbind(matrix_mean.test.seed, lapply(mean_time.seed, function(x){as.numeric(x, units="secs")}/length(folds)))
colnames(matrix_mean.test.seed) <-metric_names
rownames(matrix_mean.test.seed) <- alg_names






cat("SEED ALGORITHM NAIVE \n")
cat("On train\n")
print(matrix_mean.seed)
cat("\n")
cat("On test\n")
print(matrix_mean.test.seed)
cat("\n")



## ----warning=FALSE, cache=TRUE-------------------------------------------------------------------------------------------------------------------------------------------------------------
mean_metrics.seed_under <- list()
mean_metrics.test.seed_under <- list()
mean_time.seed_under <- list()


train.seeds_under <- list()


for(i in 1:length(folds)){
  
    train.set <- data.train[unlist(folds[i]),]
    minoritary <- subset(train.set, train.set[[OUTPUT_VAR]] == OUTPUT_MIN)
    majoritary<- subset(train.set, train.set[[OUTPUT_VAR]] == OUTPUT_MAJ)
    
    nmin = nrow(minoritary)
    nmaj = nrow(majoritary)
    
    id.majoritary <- sample(x= 1: nmaj, size = nmin/(1/prop.maj -1 )) #Indexes for majoritary class for each subset
    train.seeds_under[[i]] <- rbind(minoritary,majoritary[id.majoritary,])
  }



cat("######## Training seed algorithms undersampled ######## \n")
for(alg in 1:length(seed_algorithms)){
    mean_time.seed_under[[alg]] = 0
    
    cat(sprintf("Seed algorithm: %d\n", alg))
        
    metrics.final.seed_under <- list()
    metrics.test.seed_under <- list()
    for (i in 1:length(folds)) {

        cat(sprintf("Fold %d out of %d\n", i, length(folds)))
        cat("------------------------------------\n")
        

        
        train.set <- train.seeds_under[[i]]
        test.set <- data.train[-unlist(folds[i]),]
    
        
        start_time <- Sys.time()
        model_seed_under <- seed_algorithms[[alg]](train.set, metrics, OUTPUT_VAR, metric_optimize)
        end_time <- Sys.time()
        mean_time.seed_under[[alg]]= mean_time.seed_under[[alg]] +
        end_time - start_time
        print(end_time - start_time)
        
        metrics.final.seed_under <- append( metrics.final.seed_under, metrics_all(metricPreparedFrame(model_seed_under, test.set)))
        metrics.test.seed_under <- append( metrics.test.seed_under, metrics_all(metricPreparedFrame(model_seed_under, data.test)))
        
        cat("------------------------------------\n")
    }
  
  mean_metrics.seed_under <- append(mean_metrics.seed_under, 
        apply(matrix(unlist(metrics.final.seed_under), ncol= number_metrics, byrow=T), 2, mean))
  
    
  mean_metrics.test.seed_under <- append(mean_metrics.test.seed_under, 
        apply(matrix(unlist(metrics.test.seed_under), ncol= number_metrics, byrow=T), 2, mean))
  
  
}





matrix_mean.seed_under <- matrix(mean_metrics.seed_under, nrow = length(seed_algorithms), ncol = number_metrics, byrow = TRUE)
matrix_mean.seed_under<- cbind(matrix_mean.seed_under, lapply(mean_time.seed_under, function(x){as.numeric(x, units="secs")/length(folds)}))
colnames(matrix_mean.seed_under) <-metric_names
rownames(matrix_mean.seed_under) <- alg_names


matrix_mean.test.seed_under <- matrix(mean_metrics.test.seed_under, nrow = length(seed_algorithms), ncol = number_metrics, byrow = TRUE)
matrix_mean.test.seed_under <- cbind(matrix_mean.test.seed_under, lapply(mean_time.seed_under, function(x){as.numeric(x, units="secs")}/length(folds)))
colnames(matrix_mean.test.seed_under) <-metric_names
rownames(matrix_mean.test.seed_under) <- alg_names


cat("SEED ALGORITHM UNDERSAMPLED \n")
cat("On train\n")
print(matrix_mean.seed_under)
cat("\n")
cat("On test\n")
print(matrix_mean.test.seed_under)
cat("\n")



## ----warning=FALSE, cache=TRUE-------------------------------------------------------------------------------------------------------------------------------------------------------------
mean_metrics.SMOTE <- list()
mean_metrics.test.SMOTE <- list()
mean_time.SMOTE <- list()


train.SMOTE <- list()


for(i in 1:length(folds)){
    
    train.set <- data.train[unlist(folds[i]),]
    
    minoritary <- subset(train.set, train.set[[OUTPUT_VAR]] == OUTPUT_MIN)
    majoritary<- subset(train.set, train.set[[OUTPUT_VAR]] == OUTPUT_MAJ)
    
    nmin = nrow(minoritary)
    
    # #We reduce tha maj class to its half
    #nmaj = nrow(majoritary)/2
    #prop.min = 1-prop.maj
    #oversampling_perc = floor(100*(nmaj/((1/prop.min -1)*nmin) -1))
    
    prop.min = 1
    nmin_new = nmin*prop.min
    maj_class_perc = 100*((prop.maj/(1-prop.maj))*(1 + nmin/nmin_new))
    train.set = SMOTE(as.formula(sprintf("%s ~.", OUTPUT_VAR)) , train.set , k=5, perc.over = prop.min*100, perc.under = maj_class_perc ) 
  
    barplot(prop.table(table(train.set[[OUTPUT_VAR]])),
        col = rainbow(3),
        ylim = c(0, 1.01),
        main = "Class Distribution")
    
    minoritary <- subset(train.set, train.set[[OUTPUT_VAR]] == OUTPUT_MIN)
    majoritary<- subset(train.set, train.set[[OUTPUT_VAR]] == OUTPUT_MAJ)

    train.SMOTE[[i]] <- rbind(minoritary,majoritary)
  }



cat("######## Training seed algorithms SMOTE ######## \n")
for(alg in 1:length(seed_algorithms)){
    mean_time.SMOTE[[alg]] = 0
    
    cat(sprintf("Seed algorithm: %d\n", alg))
    
    metrics.final.SMOTE <- list()
    metrics.test.SMOTE <- list()
    for (i in 1:length(folds)) {

        cat(sprintf("Fold %d out of %d\n", i, length(folds)))
        cat("------------------------------------\n")
        
    

        
        train.set <- train.SMOTE[[i]]
        test.set <- data.train[-unlist(folds[i]),]
    
        
        start_time <- Sys.time()
        model_SMOTE <- seed_algorithms[[alg]](train.set, metrics, OUTPUT_VAR, metric_optimize)
        end_time <- Sys.time()
        mean_time.SMOTE[[alg]]= mean_time.SMOTE[[alg]] +
        end_time - start_time
        print(end_time - start_time)
        
     
        
                
        metrics.final.SMOTE <- append( metrics.final.SMOTE, metrics_all(metricPreparedFrame(model_SMOTE,  test.set)))
        metrics.test.SMOTE <- append( metrics.test.SMOTE, metrics_all(metricPreparedFrame(model_SMOTE, data.test)))
        
        cat("------------------------------------\n")
    }
  
  mean_metrics.SMOTE <- append(mean_metrics.SMOTE, 
      apply(matrix(unlist(metrics.final.SMOTE), ncol= number_metrics, byrow=T), 2, mean))
  
  mean_metrics.test.SMOTE <- append(mean_metrics.test.SMOTE, 
       apply(matrix(unlist(metrics.test.SMOTE), ncol= number_metrics, byrow=T), 2, mean))
}





matrix_mean.SMOTE <- matrix(mean_metrics.SMOTE, nrow = length(seed_algorithms), ncol = number_metrics, byrow = TRUE)
matrix_mean.SMOTE<- cbind(matrix_mean.SMOTE, lapply(mean_time.SMOTE, function(x){as.numeric(x, units="secs")/length(folds)}))


colnames(matrix_mean.SMOTE) <-metric_names
rownames(matrix_mean.SMOTE) <- alg_names


matrix_mean.test.SMOTE <- matrix(mean_metrics.test.SMOTE, nrow = length(seed_algorithms), ncol = number_metrics, byrow = TRUE)
matrix_mean.test.SMOTE <- cbind(matrix_mean.test.SMOTE, lapply(mean_time.SMOTE, function(x){as.numeric(x, units="secs")}/length(folds)))

colnames(matrix_mean.test.SMOTE) <-metric_names
rownames(matrix_mean.test.SMOTE) <- alg_names


cat("SEED ALGORITHM SMOTE\n")
cat("On train\n")
print(matrix_mean.SMOTE)
cat("\n")
cat("On test\n")
print(matrix_mean.test.SMOTE)
cat("\n")



## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
metricPreparedFrameIPIP <- function(ensemble, setToTest)
      data.frame(
          obs = setToTest[[OUTPUT_VAR]],
          pred = as.factor(prediction.fold(ensemble,setToTest[,names(setToTest) != OUTPUT_VAR])),
          prob= prediction.fold.prob(ensemble, setToTest[,names(setToTest) != OUTPUT_VAR]),
          obs.prob = as.numeric(ifelse(setToTest[[OUTPUT_VAR]] == OUTPUT_MAJ, 1, 0))
        )



## ----warning=FALSE, cache=TRUE-------------------------------------------------------------------------------------------------------------------------------------------------------------
mean_metrics.IPIPexhaustMixed <- list()
mean_metrics.test.IPIPexhaustMixed <- list()
mean_time.IPIPexhaustMixed = 0
metrics.final.IPIPexhaustMixed <- list()
metrics.test.IPIPexhaustMixed <- list()


cat("######## Training mixed exhaustive IPIP ######## \n")
 for (i in 1:length(folds)) {
    cat(sprintf("Fold %d out of %d\n",i, length(folds)))
    cat("------------------------------------\n")
      

    
    train.set <- data.train[unlist(folds[i]),]
    test.set <- data.train[-unlist(folds[i]),]
    
    nmin = sum(train.set[[OUTPUT_VAR]] == OUTPUT_MIN)
    nmaj = sum(train.set[[OUTPUT_VAR]] == OUTPUT_MAJ)
    
    np <- calculate_np( nmin, nmaj)
    p <- calculate_p(np)
    b <- calculate_b(np, nmin, nmaj)
    
    
    start_time <- Sys.time()
    ensemble.fold <- train_IPIP_exhaustive(prop.maj, OUTPUT_VAR, OUTPUT_MIN,
          OUTPUT_MAJ, train.set, test.set, seed_algorithms, prediction,  metrics,np, p, metric_optimize, metric_optimize)
    
    end_time <- Sys.time()
    mean_time.IPIPexhaustMixed = mean_time.IPIPexhaustMixed + end_time - start_time
    
    cat(sprintf("Ensembles length in exhaustive MIXED fold: %d\nLength of ensembles:", i))
    cat(unlist(lapply(ensemble.fold,length)))
    cat("\n")
    
    
    metrics.final.IPIPexhaustMixed <- append( metrics.final.IPIPexhaustMixed, 
      metrics_all(metricPreparedFrameIPIP(ensemble.fold, test.set)))
    
    metrics.test.IPIPexhaustMixed <- append( metrics.test.IPIPexhaustMixed, 
      metrics_all(metricPreparedFrameIPIP(ensemble.fold, data.test)))
      cat("------------------------------------\n")

 }


mean_metrics.IPIPexhaustMixed <- append(mean_metrics.IPIPexhaustMixed,
  apply(matrix(unlist(metrics.final.IPIPexhaustMixed), ncol= number_metrics, byrow=T), 2, mean))

matrix_mean.IPIPexhaustMixed <- matrix(mean_metrics.IPIPexhaustMixed, nrow = 1, ncol = number_metrics, byrow = TRUE)
matrix_mean.IPIPexhaustMixed <- cbind(matrix_mean.IPIPexhaustMixed, as.numeric(mean_time.IPIPexhaustMixed, units = "secs")/length(folds))

colnames(matrix_mean.IPIPexhaustMixed) <- metric_names



mean_metrics.test.IPIPexhaustMixed <- append(mean_metrics.test.IPIPexhaustMixed,
  apply(matrix(unlist(metrics.test.IPIPexhaustMixed), ncol= number_metrics, byrow=T), 2, mean))

matrix_mean.test.IPIPexhaustMixed <- matrix(mean_metrics.test.IPIPexhaustMixed, nrow = 1, ncol = number_metrics, byrow = TRUE)
matrix_mean.test.IPIPexhaustMixed <- cbind(matrix_mean.test.IPIPexhaustMixed, as.numeric(mean_time.IPIPexhaustMixed, units = "secs")/length(folds))

colnames(matrix_mean.test.IPIPexhaustMixed) <- metric_names

cat("IPIP EXHAUST MIXED\n")
cat("On train\n")
print(matrix_mean.IPIPexhaustMixed)
cat("\n")
cat("On test\n")
print(matrix_mean.test.IPIPexhaustMixed)
cat("\n")



## ----warning=FALSE, cache=TRUE-------------------------------------------------------------------------------------------------------------------------------------------------------------
mean_metrics.IPIPexhaustRepeat <- list()
mean_metrics.test.IPIPexhaustRepeat <- list()
mean_time.IPIPexhaustRepeat <- list()

cat("######## Training exhaustive repeat IPIP ######## \n")
for(alg in 1:length(seed_algorithms)){
  
   cat(sprintf("Seed algorithm: %d\n", alg))
  
   mean_time.IPIPexhaustRepeat[[alg]] = 0
   metrics.final.IPIPexhaustRepeat <- list()
   metrics.test.IPIPexhaustRepeat <- list()
   
   for (i in 1:length(folds)) {
      cat(sprintf("Fold %d out of %d\n",i, length(folds)))
      cat("------------------------------------\n")
        
      train.set <- data.train[unlist(folds[i]),]
      test.set <- data.train[-unlist(folds[i]),]
      
      nmin = sum(train.set[[OUTPUT_VAR]] == OUTPUT_MIN)
      nmaj = sum(train.set[[OUTPUT_VAR]] == OUTPUT_MAJ)
      
      np <- calculate_np( nmin, nmaj)
      p <- calculate_p(np)
      b <- calculate_b(np, nmin, nmaj)
      
      conf<- get_function_vector(b, seed_algorithms[alg])
      
      start_time <- Sys.time()
      
      ensemble.fold <- train_IPIP_exhaustive(prop.maj, OUTPUT_VAR, OUTPUT_MIN, OUTPUT_MAJ, train.set, test.set, conf, prediction,  metrics,np, p, metric_optimize, metric_optimize)
      
      end_time <- Sys.time()
      mean_time.IPIPexhaustRepeat[[alg]]= mean_time.IPIPexhaustRepeat[[alg]] + end_time - start_time
      
      cat(sprintf("Ensembles length in exhaustive REPEAT fold: %d\nLength of ensembles:", i))
      cat(unlist(lapply(ensemble.fold,length)))
      cat("\n")
    
  
      metrics.final.IPIPexhaustRepeat <- append( metrics.final.IPIPexhaustRepeat, 
      metrics_all(metricPreparedFrameIPIP(ensemble.fold, test.set)))
    
      metrics.test.IPIPexhaustRepeat <- append( metrics.test.IPIPexhaustRepeat, 
      metrics_all(metricPreparedFrameIPIP(ensemble.fold, data.test)))
          cat("------------------------------------\n")

   }
    mean_metrics.IPIPexhaustRepeat <- append(mean_metrics.IPIPexhaustRepeat, 
    apply(matrix(unlist(metrics.final.IPIPexhaustRepeat), ncol= number_metrics, byrow=T), 2, mean))
   
    mean_metrics.test.IPIPexhaustRepeat <- append(mean_metrics.test.IPIPexhaustRepeat, 
    apply(matrix(unlist(metrics.test.IPIPexhaustRepeat), ncol= number_metrics, byrow=T), 2, mean))
}

matrix_mean.IPIPexhaustRepeat <- matrix(mean_metrics.IPIPexhaustRepeat, nrow = length(seed_algorithms), ncol = number_metrics, byrow = TRUE)
matrix_mean.IPIPexhaustRepeat <- cbind(matrix_mean.IPIPexhaustRepeat, lapply(mean_time.IPIPexhaustRepeat, function(x){as.numeric(x, units="secs")/length(folds)}))
colnames(matrix_mean.IPIPexhaustRepeat) <- metric_names
rownames(matrix_mean.IPIPexhaustRepeat) <- alg_names


matrix_mean.test.IPIPexhaustRepeat <- matrix(mean_metrics.test.IPIPexhaustRepeat, nrow = length(seed_algorithms), ncol = number_metrics, byrow = TRUE)
matrix_mean.test.IPIPexhaustRepeat <- cbind(matrix_mean.test.IPIPexhaustRepeat, lapply(mean_time.IPIPexhaustRepeat, function(x){as.numeric(x, units="secs")/length(folds)}))
colnames(matrix_mean.test.IPIPexhaustRepeat) <- metric_names
rownames(matrix_mean.test.IPIPexhaustRepeat) <- alg_names


cat("IPIP EXHAUST REPEATED\n")
cat("On validation\n")
print(matrix_mean.IPIPexhaustRepeat)
cat("\n")
cat("On test\n")
print(matrix_mean.test.IPIPexhaustRepeat)
cat("\n")



## ----warning=FALSE, cache=TRUE-------------------------------------------------------------------------------------------------------------------------------------------------------------

mean_metrics.IPIPrepeated <- list()
mean_metrics.test.IPIPrepeated <- list()
mean_time.IPIPrepeated <- list()

cat("######## Training sequential repeat IPIP ######## \n")
for(alg in 1:length(seed_algorithms)){
  
   cat(sprintf("Seed algorithm: %d\n", alg))
   mean_time.IPIPrepeated[[alg]] = 0
   metrics.final.IPIPrepeated <- list()
   metrics.test.IPIPrepeated <- list()
  
   for (i in 1:length(folds)) {
      cat(sprintf("Fold %d out of %d\n",i, length(folds)))
      cat("------------------------------------\n")
  
      train.set <- data.train[unlist(folds[i]),]
      test.set <- data.train[-unlist(folds[i]),]
      
      nmin = sum(train.set[[OUTPUT_VAR]] == OUTPUT_MIN)
      nmaj = sum(train.set[[OUTPUT_VAR]] == OUTPUT_MAJ)
      
      np <- calculate_np( nmin, nmaj)
      p <- calculate_p(np)
      b <- calculate_b(np, nmin, nmaj)
    
      conf<- get_function_vector(b, seed_algorithms[alg])
      
      start_time <- Sys.time()
      ensemble.fold <- train_IPIP(prop.maj, OUTPUT_VAR, OUTPUT_MIN, OUTPUT_MAJ, train.set, test.set, conf, prediction,  metrics, b, np, p, mt, metric_optimize, metric_optimize)
      end_time <- Sys.time()
      mean_time.IPIPrepeated[[alg]]= mean_time.IPIPrepeated[[alg]] +
              end_time - start_time
      
      
      cat(sprintf("Ensembles length in sequential in fold: %d\nLength of ensembles:", i))
      cat(unlist(lapply(ensemble.fold,length)))
      cat("\n")
      
      metrics.final.IPIPrepeated <- append( metrics.final.IPIPrepeated, 
      metrics_all(metricPreparedFrameIPIP(ensemble.fold, test.set)))
    
      metrics.test.IPIPrepeated <- append( metrics.test.IPIPrepeated, 
      metrics_all(metricPreparedFrameIPIP(ensemble.fold, data.test)))
        cat("------------------------------------\n")

  }

  mean_metrics.IPIPrepeated <- append(mean_metrics.IPIPrepeated, apply(matrix(unlist(metrics.final.IPIPrepeated), ncol= number_metrics, byrow=T), 2, mean))
  mean_metrics.test.IPIPrepeated <- append(mean_metrics.test.IPIPrepeated, apply(matrix(unlist(metrics.test.IPIPrepeated), ncol= number_metrics, byrow=T), 2, mean))

}


matrix_mean.IPIPrepeated <- matrix(mean_metrics.IPIPrepeated, nrow = length(seed_algorithms), ncol = number_metrics, byrow = TRUE)
matrix_mean.IPIPrepeated <- cbind(matrix_mean.IPIPrepeated, lapply(mean_time.IPIPrepeated, function(x){as.numeric(x, units="secs")/length(folds)}))
colnames(matrix_mean.IPIPrepeated) <- metric_names
rownames(matrix_mean.IPIPrepeated) <- alg_names

matrix_mean.test.IPIPrepeated <- matrix(mean_metrics.test.IPIPrepeated, nrow = length(seed_algorithms), ncol = number_metrics, byrow = TRUE)
matrix_mean.test.IPIPrepeated <- cbind(matrix_mean.test.IPIPrepeated, lapply(mean_time.IPIPrepeated, function(x){as.numeric(x, units="secs")/length(folds)}))
colnames(matrix_mean.test.IPIPrepeated) <- metric_names
rownames(matrix_mean.test.IPIPrepeated) <- alg_names


cat("IPIP SEQUENTIAL REPEAT\n")
cat("On train:\n")
print(matrix_mean.IPIPrepeated)
cat("On test:\n")
print(matrix_mean.test.IPIPrepeated)
cat("\n")




## ---- warning=FALSE, cache = TRUE----------------------------------------------------------------------------------------------------------------------------------------------------------
mean_metrics.IPIPseqMixed<- list()
mean_metrics.test.IPIPseqMixed<- list()
mean_time.IPIPseqMixed = 0 
   
conf<- c(seed_algorithms, seed_algorithms)
cat("######## Training sequential mixed IPIP ######## \n")


   metrics.final.IPIPseqMixed <- list()
   metrics.test.IPIPseqMixed <- list()
  
   for (i in 1:length(folds)) {
      cat(sprintf("Fold %d out of %d\n",i, length(folds)))
      cat("------------------------------------\n")
  
      train.set <- data.train[unlist(folds[i]),]
      test.set <- data.train[-unlist(folds[i]),]
      
      nmin = sum(train.set[[OUTPUT_VAR]] == OUTPUT_MIN)
      nmaj = sum(train.set[[OUTPUT_VAR]] == OUTPUT_MAJ)
      
      np <- calculate_np( nmin, nmaj)
      p <- calculate_p(np)
      b <- calculate_b(np, nmin, nmaj)
    

      
      start_time <- Sys.time()
      ensemble.fold <- train_IPIP(prop.maj, OUTPUT_VAR, OUTPUT_MIN, OUTPUT_MAJ, train.set, test.set, conf, prediction,  metrics, b, np, p, mt, metric_optimize, metric_optimize)
      end_time <- Sys.time()
      mean_time.IPIPseqMixed= mean_time.IPIPseqMixed +
              end_time - start_time
      
      
      cat(sprintf("Ensembles length in sequential in fold: %d\nLength of ensembles:", i))
      cat(unlist(lapply(ensemble.fold,length)))
      cat("\n")
    
      
      metrics.final.IPIPseqMixed <- append( metrics.final.IPIPseqMixed, 
      metrics_all(metricPreparedFrameIPIP(ensemble.fold, test.set)))
    
      metrics.test.IPIPseqMixed <- append( metrics.test.IPIPseqMixed, 
      metrics_all(metricPreparedFrameIPIP(ensemble.fold, data.test)))
        cat("------------------------------------\n")




}

mean_metrics.IPIPseqMixed <- append(mean_metrics.IPIPseqMixed, apply(matrix(unlist(metrics.final.IPIPseqMixed), ncol= number_metrics, byrow=T), 2, mean))
matrix_mean.IPIPseqMixed <- matrix(mean_metrics.IPIPseqMixed, nrow = 1, ncol = number_metrics, byrow = TRUE)
matrix_mean.IPIPseqMixed <- cbind(matrix_mean.IPIPseqMixed,lapply(mean_time.IPIPseqMixed, function(x){as.numeric(x, units="secs")/length(folds)}))
colnames(matrix_mean.IPIPseqMixed) <- metric_names

mean_metrics.test.IPIPseqMixed <- append(mean_metrics.test.IPIPseqMixed, apply(matrix(unlist(metrics.test.IPIPseqMixed), ncol= number_metrics, byrow=T), 2, mean))
matrix_mean.test.IPIPseqMixed <- matrix(mean_metrics.test.IPIPseqMixed, nrow = 1, ncol = number_metrics, byrow = TRUE)
matrix_mean.test.IPIPseqMixed <- cbind(matrix_mean.test.IPIPseqMixed,lapply(mean_time.IPIPseqMixed, function(x){as.numeric(x, units="secs")/length(folds)}))
colnames(matrix_mean.test.IPIPseqMixed) <- metric_names


cat("IPIP SEQUENTIAL MIXED\n")
cat("On train:\n")
print(matrix_mean.IPIPseqMixed)
cat("On test:\n")
print(matrix_mean.test.IPIPseqMixed)
cat("\n")


## ---- results= 'hold'----------------------------------------------------------------------------------------------------------------------------------------------------------------------

print_with_asterisks <- function(matrix, best_values, elem) {
    num_spaces = 11
    for(j in 0:nrow(matrix)){
        for (i in 0:ncol(matrix)) {
          if(j==0 && i==0){
            str_to_cat = " "
          }
          else if(i==0){
            str_to_cat = alg_names[j]
          }
          else if(j == 0){
            str_to_cat = metric_names[i]          
          }
          else{
              if(is.nan(matrix[[j,i]]) == FALSE && matrix[[j, i]] == best_values[[i,elem]]) {
                if ((i!=13 && best_values[[i,elem]] == max(unlist(best_values[i,]))) || (i==13&&best_values[[i,elem]]==min(unlist(best_values[i,])))){
                  str_to_cat = sprintf("%f**", matrix[[j,i]])
                } else {
                  str_to_cat = sprintf("%f*", matrix[[j,i]])
                }
              }
              else{
                str_to_cat = sprintf("%f", matrix[[j,i]])
              }
              
              
              
          }
          cat(str_to_cat)
          cat(strrep(" ", num_spaces - nchar(str_to_cat)))
      }
      cat("\n")
    }
    cat("\n")
}

save(matrix_mean.IPIPrepeated, matrix_mean.IPIPseqMixed, matrix_mean.IPIPexhaustRepeat,
     matrix_mean.IPIPexhaustMixed, matrix_mean.seed, matrix_mean.seed_under, matrix_mean.SMOTE,matrix_mean.test.IPIPrepeated, matrix_mean.test.IPIPseqMixed, matrix_mean.test.IPIPexhaustRepeat,matrix_mean.test.IPIPexhaustMixed, matrix_mean.test.seed, matrix_mean.test.seed_under, matrix_mean.test.SMOTE,
     file = "./saved_matrices.RData")

get_better_data <- function(matrix){
  best_data <- list()
  for (col in metric_names){
    if(col == "TIME"){
      best_data <- append(best_data ,min(unlist(matrix[,col]), na.rm=TRUE))
    }
    else{
      best_data <- append(best_data, max(unlist(matrix[,col]), na.rm=TRUE))
    }
  }
  return(best_data)
}


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("----------------------------------------------\n")
cat("Comparing on train ALL\n")

best_values <- sapply( list(
  matrix_mean.IPIPrepeated,
  matrix_mean.IPIPseqMixed,
  matrix_mean.IPIPexhaustRepeat,
  matrix_mean.IPIPexhaustMixed,
  matrix_mean.seed,
  matrix_mean.seed_under,
  matrix_mean.SMOTE
), get_better_data)

cat("IPIP SEQUENTIAL REPEAT\n")
print_with_asterisks(matrix_mean.IPIPrepeated,best_values,1)

cat("IPIP SEQUENTIAL MIXED\n")
print_with_asterisks(matrix_mean.IPIPseqMixed,best_values,2)

cat("IPIP EXHAUST REPEAT\n")
print_with_asterisks(matrix_mean.IPIPexhaustRepeat, best_values, 3)

cat("IPIP EXHAUST MIXED\n")
print_with_asterisks(matrix_mean.IPIPexhaustMixed, best_values, 4)

cat("SEED\n")
print_with_asterisks(matrix_mean.seed, best_values, 5)

cat("SEED UNDERSAMPLED\n")
print_with_asterisks(matrix_mean.seed_under, best_values, 6)

cat("SMOTE\n")
print_with_asterisks(matrix_mean.SMOTE, best_values, 7)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("----------------------------------------------\n")
cat("Comparing on test ALL\n")

best_values <- sapply( list(
  matrix_mean.test.IPIPrepeated,
  matrix_mean.test.IPIPseqMixed,
  matrix_mean.test.IPIPexhaustRepeat,
  matrix_mean.test.IPIPexhaustMixed,
  matrix_mean.test.seed,
  matrix_mean.test.seed_under,
  matrix_mean.test.SMOTE
), get_better_data)

cat("IPIP SEQUENTIAL REPEAT\n")
print_with_asterisks(matrix_mean.test.IPIPrepeated,best_values,1)

cat("IPIP SEQUENTIAL MIXED\n")
print_with_asterisks(matrix_mean.test.IPIPseqMixed,best_values,2)

cat("IPIP EXHAUST REPEAT\n")
print_with_asterisks(matrix_mean.test.IPIPexhaustRepeat, best_values, 3)

cat("IPIP EXHAUST MIXED\n")
print_with_asterisks(matrix_mean.test.IPIPexhaustMixed, best_values, 4)

cat("SEED\n")
print_with_asterisks(matrix_mean.test.seed, best_values, 5)

cat("SEED UNDERSAMPLED\n")
print_with_asterisks(matrix_mean.test.seed_under, best_values, 6)

cat("SMOTE\n")
print_with_asterisks(matrix_mean.test.SMOTE, best_values, 7)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("----------------------------------------------\n")
cat("Comparing on train ONLY IPIP\n")

best_values <- sapply( list(
  matrix_mean.IPIPrepeated,
  matrix_mean.IPIPseqMixed,
  matrix_mean.IPIPexhaustRepeat,
  matrix_mean.IPIPexhaustMixed
), get_better_data)

cat("IPIP SEQUENTIAL REPEAT\n")
print_with_asterisks(matrix_mean.IPIPrepeated,best_values,1)

cat("IPIP SEQUENTIAL MIXED\n")
print_with_asterisks(matrix_mean.IPIPseqMixed,best_values,2)

cat("IPIP EXHAUST REPEAT\n")
print_with_asterisks(matrix_mean.IPIPexhaustRepeat, best_values, 3)

cat("IPIP EXHAUST MIXED\n")
print_with_asterisks(matrix_mean.IPIPexhaustMixed, best_values, 4)



## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("----------------------------------------------\n")
cat("Comparing on test ONLY IPIP\n")

best_values <- sapply( list(
  matrix_mean.test.IPIPrepeated,
  matrix_mean.test.IPIPseqMixed,
  matrix_mean.test.IPIPexhaustRepeat,
  matrix_mean.test.IPIPexhaustMixed
), get_better_data)

cat("IPIP SEQUENTIAL REPEAT\n")
print_with_asterisks(matrix_mean.test.IPIPrepeated,best_values,1)

cat("IPIP SEQUENTIAL MIXED\n")
print_with_asterisks(matrix_mean.test.IPIPseqMixed,best_values,2)

cat("IPIP EXHAUST REPEAT\n")
print_with_asterisks(matrix_mean.test.IPIPexhaustRepeat, best_values, 3)

cat("IPIP EXHAUST MIXED\n")
print_with_asterisks(matrix_mean.test.IPIPexhaustMixed, best_values, 4)

