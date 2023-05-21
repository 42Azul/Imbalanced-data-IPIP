E <- list() # Modelo final (ensemble de ensembles)
set.seed(37)

tC <- trainControl(
  summaryFunction = metricas,
  method = "cv",
  number = b,
  allowParallel = TRUE,
  classProbs = TRUE
)
for(k in 1:p){
  Ek <- list() # Ensemble de modelos k-ésimo
  i <- 0 # Contador para el número de intentos de ampliar el ensemble
  # Conjunto de datos perfectamente balanceado:
  df <- dfs[[k]]
  while(length(Ek)<=b && i<mt(length(Ek))){
    # Seleccionamos muestras para entrenar el modelo de random forest
    pob <- which(df$income == "pobre")
    ricos <- which(df$income == "rico")
    
    ind.train <- c(
      sample(pob, size = round(np*prop.mayoritaria/(1-prop.mayoritaria)), replace = TRUE),
      sample(ricos, size = np, replace = TRUE)
    )
    
    cl <- makeCluster(detectCores()-2)
    registerDoParallel(cl)
    
    rf <- train(
      x = df[ind.train,-1],
      num.trees = 200,
      importance = "impurity",
      y = df$income[ind.train],
      method = "ranger",
      metric = "KAPPA",
      maximize = T,
      trControl = tC,
      tuneGrid = hiperparametros
    )
    
    stopCluster(cl)
    
    # Evaluamos el ensemble actual (sin el nuevo modelo)
    
    metricas.ensemble <-
      if (length(Ek)==0){
        u <- -Inf
        names(u) <- "KAPPA"
        u
      } else metricas(data.frame(
        obs = data.test$income,
        pred= prediccion(Ek, data.test[-1])
      ))
    
    Ek[[length(Ek)+1]] <- rf
    # Evaluamos el ensemble formado al añadir el nuevo modelo
    metricas.ensemble.2 <- metricas(data.frame(
      obs = data.test$income,
      pred= prediccion(Ek, data.test[-1])
    ))
    # Comparamos las metricas
    if(metricas.ensemble.2["KAPPA"] <= metricas.ensemble["KAPPA"]){ # Si el ensemble no mejora con el nuevo modelo...
      i <- i+1
      Ek[[length(Ek)]] <- NULL
    } else{ # En caso de ampliar el ensemble, reseteamos las oportunidades de cara a una nueva ampliación
      i <- 0
    }
    
  } # Fin del WHILE (hemos terminado de construir el ensemble k-ésimo)
  
  # Guardamos la información del ensemble k-ésimo
  E[[length(E)+1]] <- Ek
  
} # FIN. Hemos terminado de contruir el ensemble final

saveRDS(E,"./TrainCompleto.rds")