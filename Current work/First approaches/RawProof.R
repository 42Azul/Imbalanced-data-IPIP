
library(doParallel)
library(DataExplorer)
library(dplyr)
library(pROC)
library(caret)
datos <- read.csv("../data_tfg.csv")


### PREPROCESAMIENTO

ind.cualit <- c(which(names(datos) == "SITUACION"),which(names(datos)=="SEXO"), which(names(datos)=="DM"):which(names(datos)=="DC"))


for(i in ind.cualit) datos[,i] <- as.factor(datos[, i])


# Create a list of indices for k-fold cross-validation
set.seed(123)
folds <- createFolds(datos$SITUACION, k = 5)
### SPLIT ENTRENAMIENTO
  subconjunto <- createDataPartition(y = datos$SITUACION, times = 1, p = 0.7, list = FALSE)
  train.set <- datos[subconjunto, ]
  test.set <- datos[-subconjunto, ]





prediccion <- function(conj.model, x, q = 0.75){ #q=0.75, pero se deberían probar valores como 0.5, 0.25, 0.75...
  pred <- data.frame(matrix(nrow=nrow(x),ncol=0))
  for(modelo in conj.model) pred <- cbind(pred, predict(modelo,x))
  pred <- apply(pred, 1, function(x) prop.table(table(x))["CURADO"])
  ifelse(is.na(pred) | pred<q, "FALLECIDO", "CURADO")
}


#Las métricas, especialmente KAPPA (igual ROC después)

metricas <- function(data, lev = levels(as.factor(data$obs)), model = NULL){
  c(
    ACCURACY = MLmetrics::Accuracy(data[, "pred"], data[, "obs"]),
    SENS = sensitivity(data[, "pred"],data[, "obs"],positive="FALLECIDO",negative="CURADO"),
    SPEC = specificity(data[, "pred"], data[, "obs"],positive="FALLECIDO",negative="CURADO"),
    PPV = posPredValue(data[, "pred"], data[, "obs"],positive="FALLECIDO",negative="CURADO"),
    NPV = negPredValue(data[, "pred"], data[, "obs"],positive="FALLECIDO",negative="CURADO"),
    KAPPA = psych::cohen.kappa(cbind(data[, "obs"],data[, "pred"]))$kappa,
    BAL_ACC = (sensitivity(data[, "pred"],data[, "obs"],positive="FALLECIDO",negative="CURADO") + specificity(data[, "pred"], data[, "obs"],positive="FALLECIDO",negative="CURADO"))/2
  )
}



function_vector_rlog <-c()
for(i in 1:5){
  train_model <- function(df.train, metricas) {
    
    tC <- trainControl(method = 'repeatedcv',summaryFunction = metricas,
                       number = 5,repeats =  5,search = 'random')
    
    method <- "glmnet"
    metric <- "KAPPA"
    maximize <- T
    
    # Entrenamos la Rlog
    rlog <- train(SITUACION ~ .,
                  data = df.train,
                  method = "glmnet",
                  family = 'binomial',
                  metric = "KAPPA",
                  maximize = T,
                  trControl = tC,
    )
    
    return(rlog)
  }
  
  function_vector_rlog <- append(function_vector_rlog, train_model)
}


#Ahora lo mismo pero con SVM

function_vector_svm <-c()
for(i in 1:5){
  train_model <- function(df.train, metricas) {
    
    tC <- trainControl(
      summaryFunction = metricas,
      allowParallel = TRUE,
      classProbs = TRUE
    )
    
    method <- "svmLinear"
    metric <- "KAPPA"
    maximize <- T
    
    # Entrenamos la svm
    svm <- train(
      SITUACION ~ .,
      data = df.train,
      method = method,
      metric = metric,
      maximize = maximize,
      trControl = tC,
    )
    
    return(svm)
  }
  
  function_vector_svm <- append(function_vector_svm, train_model)
}

#Por último con GBM

function_vector_gbm <-c()
for(i in 1:5){
  train_model <- function(df.train, metricas) {
    
    tC <- trainControl(
      summaryFunction = metricas,
      allowParallel = TRUE,
      classProbs = TRUE
    )
    
    method <- "gbm"
    metric <- "KAPPA"
    maximize <- T
    
    # Entrenamos el gradient boosting
    gbm <- train(SITUACION ~ .,
                 data = df.train,
                 method = method,
                 metric = metric,
                 maximize = maximize,
                 trControl = tC,
    )
    
    return(gbm)
  }
  
  function_vector_gbm <- append(function_vector_gbm, train_model)
}



source("train_IPIP.R")

prop.mayoritaria = 0.55





E <- train_IPIP(prop.mayoritaria, "SITUACION", "FALLECIDO", "CURADO", train.set, test.set, 
                function_vector_rlog, prediccion,  metricas)

saveRDS(E,"./TrainCompletoRLOG.rds")



E <- train_IPIP(prop.mayoritaria, "SITUACION", "FALLECIDO", "CURADO", train.set, test.set, 
                function_vector_svm, prediccion,  metricas)

saveRDS(E,"./TrainCompletoSVM.rds")


E <- train_IPIP(prop.mayoritaria, "SITUACION", "FALLECIDO", "CURADO", train.set, test.set, 
                function_vector_gbm, prediccion,  metricas)

saveRDS(E,"./TrainCompletoGBM.rds")


#Ahora veamos los resultados:
  
ensemble.ranger <- readRDS("TrainCompletoRANGER.rds")
ensemble.rlog <- readRDS("TrainCompletoRLOG.rds")
ensemble.svm <- readRDS("TrainCompletoSVM.rds")
ensemble.gbm <- readRDS("TrainCompletoGBM.rds")

#Veamos los modelos que son


unlist(lapply(ensemble.rlog,length))
unlist(lapply(ensemble.svm,length))
unlist(lapply(ensemble.gbm,length))

#Ahora hacemos el ensemble y predecimos:
  
prediccion.final <- function(ensemble, x, q = 0.5){
  # Colocamos en cada fila de un conjunto de datos todas las predicciones para una muestra
  pred <- as.data.frame(lapply(ensemble, function(e) prediccion(e,x)))
  pred <- apply(pred, 1, function(x) prop.table(table(x))["CURADO"])
  ifelse(is.na(pred) | pred<q, "FALLECIDO", "CURADO")
}



metricas.final.rlog<- metricas(data.frame(
  obs = test.set$SITUACION,
  pred= as.factor(prediccion.final(ensemble.rlog, test.set[-1]))
))

metricas.final.svm<- metricas(data.frame(
  obs = test.set$SITUACION,
  pred= as.factor(prediccion.final(ensemble.svm, test.set[-1]))
))

metricas.final.gbm<- metricas(data.frame(
  obs = test.set$SITUACION,
  pred= as.factor(prediccion.final(ensemble.gbm, test.set[-1]))
))



print("Metricas regresion logistica")
metricas.final.rlog

print("Metricas SVM")
metricas.final.svm

print("Metricas GBM")
metricas.final.gbm



