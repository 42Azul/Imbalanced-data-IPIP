
train_IPIP <- function( prop.mayoritaria, OUTPUT, min_str, may_str, train.set, test.set, 
                        function_vector, prediccion,  metricas, prob_codo=0.75, alpha = 0.01, alpha_b= .01) 

  
{
  
  mt <- function(n) { ceiling((b-n) / 3) }
  
  nmin = sum(train.set[[OUTPUT]] == min_str)
  nmay = sum(train.set[[OUTPUT]] == may_str)
  
  

  np <- ceiling(nmin*prob_codo)
  p <- ceiling(log(alpha)/(log(1-1/nmin)*np))
  
  maySubSize <-  round(np*prop.mayoritaria/(1-prop.mayoritaria))
  minSubSize <-  np

  dfs <- list()
  
  minoritario = train.set %>% dplyr::filter(.data[[OUTPUT]] == min_str)
  mayoritario = train.set %>% dplyr::filter(.data[[OUTPUT]] == may_str)

  for(k in 1:p){
    id.minoritaria <- sample(x = 1:nmin, size = minSubSize) #Índices de clase minoritaria para cada subconjunto
    id.mayoritaria <- sample(x= 1:nmay, size = maySubSize) #Índices de la clase mayoritaria para cada subconjunto
    
    dfs[[k]] <- rbind(minoritario[id.minoritaria,],mayoritario[id.mayoritaria,])
  }
  
  
  #La b es el número de modelos por conjunto
  b <- ceiling(log(alpha_b)/(log(1-1/np)*np))

  if(length(function_vector) != b){
    print(sprintf("Different size for function vector:%d and value of b: %d", length(function_vector), b))
    return(NA)
  }
  
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
      sample(mayoritaria, size = maySubSize, replace = TRUE),
      sample(minoritaria, size = minSubSize, replace = TRUE) 
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
