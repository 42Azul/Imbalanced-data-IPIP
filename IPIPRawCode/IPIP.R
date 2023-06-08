
#Function to calculate number of positive samples, np
calculate_np <- function( nmin, nmaj,  elbow_prob=0.75) 
{
  np <- round(nmin*elbow_prob)
  return(np)
}

#Function to calculate number of partitions, p
calculate_p <- function( np, prob_codo=0.75, alpha_p= .01) 
{
  p <- ceiling(log(alpha_p)/(log(1-1/np)*np))
  return(p)
}


#Function to calculate maximum ensemble size, b
calculate_b <- function(np, nmin, nmaj, elbow_prob=0.75, alpha_b= .01) 
{
  
  b <- ceiling(log(alpha_b)/(log(1-1/nmin)*np))
  return(b)
}


#Function to calculate maximum number of consecutive attempts to enlarge ensemble
mt <- function(b, n) { ceiling((b-n) / 3) }



# Function to update the ensemble with a new model adn enlarge its prediction values 

prediction.add <- function(old_ensemble, new_element, prediction_current = FALSE, val.test = FALSE){
  # if no current prediction is passed, generate prediction using new model
  if(prediction_current == FALSE){
    prediction_current = as.numeric(unlist(predict(new_element, val.test)) == OUTPUT_MAJ)
  }
  # If the old ensemble is not empty, update the ensemble and the model count
  if(length(old_ensemble$ensemble)!=0){
    prev_count =  old_ensemble$model_count
    l <- list("ensemble"= c(old_ensemble$ensemble,list(new_element)), "model_count"= prev_count + prediction_current)
    return(l)
  }
  else{  # else initialize ensemble with the new model and the prediction count
    l <- list("ensemble"= list(new_element), "model_count"=prediction_current)
    return(l)
  }
}


# Function to get the final predicted class by comparing the model count with a threshold (q times the ensemble length)
prediction.predicted <- function(predicted, q = default_q) {
  
  nElems <- length(predicted$ensemble) 
  pred <- ifelse(predicted$model_count  < q*nElems, OUTPUT_MIN, OUTPUT_MAJ)
  return(pred)
}





train_IPIP_seq <- function( prop.maj, OUTPUT, min_str, maj_str, train.set, val.test, 
                        configuration, prediction, metrics, b, np,  p, mt, TRAIN_METRIC = "BAL_ACC", SELECT_METRIC = "BAL_ACC"){
  
  nmin = sum(train.set[[OUTPUT]] == min_str)
  nmaj = sum(train.set[[OUTPUT]] == maj_str)
  
  
  
  majSubSize <-  round(np*prop.maj/(1-prop.maj))
  minSubSize <-  np
  
  #Icludes from 1 to p all ensembles of the set
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
      
      new_Ek <- prediction.add(Ek, model, val.test = val.test[colnames(val.test)!=OUTPUT])
      pred = prediction.predicted(new_Ek)
      metrics.ensemble.new <- metrics(data.frame(
        obs = val.test[[OUTPUT]],
        pred= as.factor(pred)
      ))
      #We check if the new model changes the result. If it does not we start trying again
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

#Creates the powerset of the trained models as well as increases their count in a dynamic approach.
powerset = function(s, val.test){
  len = length(s)
  l = c()
  vector(length=2^len);l[[1]]=numeric()
  counter = 1L
  for(x in 1L:len){
    prediction_current =as.numeric(unlist(predict(s[x], val.test)) == OUTPUT_MAJ)
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


train_IPIP_exhaustive <- function( prop.maj, OUTPUT, min_str, maj_str, train.set, val.test, configuration, prediction, metrics, np,  p, TRAIN_METRIC = "KAPPA", SELECT_METRIC = "KAPPA"){
  
  
  
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
    
    all_sets <- powerset(Ek, as.data.frame(val.test[, -which(names(val.test) == OUTPUT)]))
    
    ##If we want to force bigger sets before smaller ones
    #lengths <- sapply(all_sets, function(x) length(x))
    # Reorder the list based on the lengths
    # all_sets <- all_sets[order(lengths, decreasing = TRUE)]
    
    
    
    max_metric <- function(model){
      pred = prediction.predicted(model)
      met <- metrics(data.frame(
        obs = val.test[[OUTPUT]],
        pred = as.factor(pred)
      ))[SELECT_METRIC]
      return(met)
    }
    
    
    
    #Run the max metric check in parallel mode
    max_metric_list <- lapply(all_sets, max_metric)
    
    max_metric_value <- max(unlist(max_metric_list))
    max_set <- all_sets[[which.max(unlist(max_metric_list))]]$ensemble
    cat(sprintf("Max metric %s of ensemble is %f with length %d\n", SELECT_METRIC, max_metric_value, length(max_set)))
    E[[length(E)+1]] <- max_set
    
  }
  return(E);
}



trainIPIP <- function(train.set, val.set, OUTPUT_VAR, OUTPUT_MIN, OUTPUT_MAJ, conf, metrics, metric_optimize, prop.maj = 0.55, exhaustive = FALSE,  alpha = .01, metric_select = metric_optimize){
  
  nmin = sum(train.set[[OUTPUT_VAR]] == OUTPUT_MIN)
  nmaj = sum(train.set[[OUTPUT_VAR]] == OUTPUT_MAJ)
  
  np <- calculate_np( nmin, nmaj)
  p <- calculate_p(np, alpha)
  b <- calculate_b(np, nmin, nmaj,alpha)
  
  ensemble.fold = 0
  
  if(exhaustive == TRUE) return(train_IPIP_exhaustive(prop.maj, OUTPUT_VAR, OUTPUT_MIN,
                                         OUTPUT_MAJ, train.set, val.set, seed_algorithms, metrics,np, p, metric_optimize, metric_select))
  
  return(train_IPIP_seq(prop.maj, OUTPUT_VAR, OUTPUT_MIN, OUTPUT_MAJ, train.set, val.set, conf,  metrics, b, np, p, mt, metric_optimize, metric_select))
}

