
train_IPIP <- function( p,b, np, prop.mayoritaria, OUTPUT, min_str, may_str, train.set, test.set,dfs, 
                        function_vector, prediccion,  metricas) 
{
  
  may_size <-  round(np*prop.mayoritaria/(1-prop.mayoritaria))
  min_size <-  np
  

E <- list() # Modelo final (ensemble de ensembles)
set.seed(42)

for(k in 1:p){
  Ek <- list() # Ensemble de modelos k-ésimo
  i <- 0 # Contador para el número de intentos de ampliar el ensemble
  
  # Conjunto de datos balanceado:
  df <- dfs[[k]]
  
  while(length(Ek)<=b && i<mt(length(Ek))){
    # Seleccionamos muestras para entrenar el modelo
    
    mayoritaria <- which(df[[OUTPUT]] == may_str)
    minoritaria <- which(df[[OUTPUT]] == min_str)
    ind.train <- c(
      sample(mayoritaria, size = may_size, replace = TRUE),
      sample(minoritaria, size = min_size, replace = TRUE) 
    )

    modelo <- function_vector[[b]](df[ind.train,], metricas)
    

    
    # Evaluamos el ensemble actual (sin el nuevo modelo)
    
    metricas.ensemble <-
      if (length(Ek)==0){
        u <- -Inf;
        names(u) <- "KAPPA";
        u;
      } else metricas(data.frame(
        obs = test.set[[OUTPUT]],
        pred = as.factor(prediccion(Ek, test.set[colnames(test.set)!=OUTPUT]))
    ))

    Ek[[length(Ek)+1]] <- modelo
    
    # Evaluamos el ensemble formado al añadir el nuevo modelo
    metricas.ensemble.new <- metricas(data.frame(
      obs = test.set[[OUTPUT]],
      pred= as.factor(prediccion(Ek, test.set[colnames(test.set)!=OUTPUT]))
    ))
    
    
    # Comparamos las metricas
    if(metricas.ensemble.new["KAPPA"] <= metricas.ensemble["KAPPA"]){ # Si el ensemble no mejora con el nuevo modelo...
      i <- i+1
      Ek[[length(Ek)]] <- NULL
    } else{ # En caso de ampliar el ensemble, reseteamos las oportunidades de cara a una nueva ampliación
      i <- 0
    }
    
  } # Fin del WHILE (hemos terminado de construir el ensemble k-ésimo)
  
  # Guardamos la información del ensemble k-ésimo
  E[[length(E)+1]] <- Ek
  
} # FIN. Hemos terminado de contruir el ensemble final
return(E);
}
