
library(doParallel)
library(DataExplorer)
library(dplyr)
library(pROC)
library(caret)
source("train_IPIP.R")
source("imbalancedFold.R")



kFoldTrain <- function(function_vector, k ){

datos <- read.csv("../data_tfg.csv")


### PREPROCESAMIENTO

ind.cualit <- c(which(names(datos) == "SITUACION"),which(names(datos)=="SEXO"), which(names(datos)=="DM"):which(names(datos)=="DC"))


for(i in ind.cualit) datos[,i] <- as.factor(datos[, i])

set.seed(123)
metricas.final <- list()


##FUNCIONES A USAR

prediccion.final <- function(ensemble, x, q = 0.5){
  # Colocamos en cada fila de un conjunto de datos todas las predicciones para una muestra
  pred <- as.data.frame(lapply(ensemble, function(e) prediccion(e,x)))
  pred <- apply(pred, 1, function(x) prop.table(table(x))["CURADO"])
  ifelse(is.na(pred) | pred<q, "FALLECIDO", "CURADO")
}


prediccion <- function(conj.model, x, q = 0.75){
  pred <- data.frame(matrix(nrow=nrow(x),ncol=0))
  for(modelo in conj.model) pred <- cbind(pred, predict(modelo,x))
  pred <- apply(pred, 1, function(x) prop.table(table(x))["CURADO"])
  ifelse(is.na(pred) | pred<q, "FALLECIDO", "CURADO")
}


#Las métricas, especialmente KAPPA

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

prop.mayoritaria <- 0.55


### SPLIT ENTRENAMIENTO
subconjunto <-imbalancedFold(datos, k,"SITUACION", "FALLECIDO")
for (fold in 1:k) {
  
  train.set <- datos[unlist(subconjunto[fold]),]
  test.set <- datos[-unlist(subconjunto[fold]),]
  
  ensemble.ranger <- train_IPIP(prop.mayoritaria, "SITUACION", "FALLECIDO", "CURADO", train.set, test.set, 
                  function_vector, prediccion,  metricas)
  
  print(unlist(lapply(ensemble.ranger,length)))
  
  metricas.final<- append(metricas(data.frame(
    obs = test.set$SITUACION,
    pred= as.factor(prediccion.final(ensemble.ranger, test.set[-1]))
  )), metricas.final)

  print (metricas.final[fold])

}

mean_metricas <- apply(matrix(unlist(metricas.final), ncol= 7, byrow=T), 2, mean)
print(mean_metricas)
return(mean_metricas)
}