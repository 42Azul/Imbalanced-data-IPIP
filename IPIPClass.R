library(R6)

#train set:       Set from which train the inside ensembles
#val.set:         Validation set for ensemble selection
#OUTPUT_VAR:      Variable where is our objective located
#OUTPUT_MIN:      Value of the minoritary class in the output characteristic
#OUTPUT_MAJ:      Value of the majoritary class in the output characteristic
#conf:            Vector of algorithms to be executed calling the TRAIN_METRIC metric and the train set
#metrics:         Metric function which calculates, given a prediction in the form of data[,obs] data[,pred]
#TRAIN_METRIC: Metric will be used to train the model
#prop_maj:        Proportion of the majoritary class
#exhaustive:      Boolean variable which checks if the training is made in sequential or exhaustive mode
#SELECT_METRIC:   Second metric stablished in order to check if one ensemble is better than other. They are equal by default.

#Returns -> List of trained models. Prediction over those models is to be made with the prediction.fold functions.

IPIP <- R6Class("IPIP",
                public = list(
                  OUTPUT_VAR = NULL,
                  OUTPUT_MIN = NULL,
                  OUTPUT_MAJ = NULL,
                  conf = NULL,
                  metrics = NULL,
                  TRAIN_METRIC = NULL,
                  prop_maj = 0.55,
                  exhaustive = FALSE,
                  alpha_p = 0.01,
                  alpha_b = 0.01,
                  SELECT_METRIC = NULL,
                  np = NULL,
                  p = NULL,
                  b = NULL,
                  ensemble_model = NULL,
                  percentage_per_sample = 0.7,
                  q = 0.7,
                  
                  initialize = function(OUTPUT_VAR, OUTPUT_MIN, OUTPUT_MAJ, conf, metrics, 
                                        TRAIN_METRIC, prop_maj = 0.55, exhaustive = FALSE, 
                                        alpha_p = 0.01, alpha_b = 0.01, SELECT_METRIC = TRAIN_METRIC) {
                    self$OUTPUT_VAR <- OUTPUT_VAR
                    self$OUTPUT_MIN <- OUTPUT_MIN
                    self$OUTPUT_MAJ <- OUTPUT_MAJ
                    self$conf <- conf
                    self$metrics <- metrics
                    self$TRAIN_METRIC <- TRAIN_METRIC
                    self$prop_maj <- prop_maj
                    self$exhaustive <- exhaustive
                    self$alpha_p <- alpha_p
                    self$alpha_b <- alpha_b
                    self$SELECT_METRIC <- SELECT_METRIC
                  },
                  train = function(train.set, val.set) {
                    nmin = sum(train.set[[self$OUTPUT_VAR]] == self$OUTPUT_MIN)
                    nmaj = sum(train.set[[self$OUTPUT_VAR]] == self$OUTPUT_MAJ)
                    self$np <- self$calculate_np(nmin, nmaj)
                    self$p <- self$calculate_p()
                    self$b <- self$calculate_b(nmin)
                    
                    if(self$exhaustive == TRUE) {
                      self$ensemble_model <- self$train_IPIP_exhaustive(train.set, val.set)
                    } else {
                      self$ensemble_model <- self$train_IPIP_seq(train.set, val.set)
                    }
                    invisible(self)
                  },
                  
                  predict = function(x) {
                    prediction.fold(self$ensemble_model, x, self$default_q)
                  },
                  
                  #Function to calculate number of positive samples, np
                  calculate_np= function(nmin, nmaj) 
                  {
                    np <- round(nmin*self$percentage_per_sample)
                    return(np)
                  },
                  #Function to calculate number of partitions, p
                  calculate_p = function() 
                  {
                    p <- ceiling(log(self$alpha_p)/(log(1-1/self$np)*self$np))
                    return(p)
                  },
                  
                  #Function to calculate maximum ensemble size, b
                  calculate_b = function( nmin) 
                  {
                    
                    b <- ceiling(log(self$alpha_b)/(log(1-1/nmin)*self$np))
                    return(b)
                  },
                  
                  
                  #Function to calculate maximum number of consecutive attempts to enlarge ensemble
                  mt = function(b, n) { ceiling((b-n) / 3) },
                  
                  
                  
                  # Function to update the ensemble with a new model adn enlarge its prediction values 
                  
                  prediction.add = function(old_ensemble, new_element, prediction_current = FALSE, val.test = FALSE){
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
                  },
                  
                  
                  # Function to get the final predicted class by comparing the model count with a threshold (q times the ensemble length)
                  prediction.predicted = function(predicted){
                    
                    nElems <- length(predicted$ensemble) 
                    pred <- ifelse(predicted$model_count  < self$q*nElems, self$OUTPUT_MIN, self$OUTPUT_MAJ)
                    return(pred)
                  },
                  
                  #Functions used to predict a final trained model
                  
                  #Function to predict over a single ensemble
                  prediction = function(conj.model, x){ 
                    pred <- data.frame(matrix(nrow=nrow(x),ncol=0))
                    for(model in conj.model) pred <- cbind(pred, predict(model,x))
                    nElems = ncol(pred)
                    nElems <- ncol(pred)
                    counts <- rowSums(pred == OUTPUT_MAJ)
                    pred <- counts
                    ifelse(is.na(pred) | pred<self$q*nElems, self$OUTPUT_MIN, self$OUTPUT_MAJ)
                  },
                  
                  #Function to predict into a fold of many ensembles
                  prediction.fold = function(ensemble, x){
                    pred <- as.data.frame(lapply(ensemble, function(e) prediction(e,x)))
                    nElems = ncol(pred)
                    nElems <- ncol(pred)
                    counts <- rowSums(pred == self$OUTPUT_MAJ)
                    pred <- counts
                    ifelse(is.na(pred) | pred<self$q*nElems, self$OUTPUT_MIN, self$OUTPUT_MAJ)
                  },
                  
                  
                  #Functions to calculate the probability of the predictions
                  prediction.prob = function(ensemble, x){
                    pred <- data.frame(matrix(nrow=nrow(x),ncol=0))
                    for(model in ensemble) pred <- cbind(pred, predict(model,x,type="prob")[[self$OUTPUT_MAJ]])
                    pred<- rowSums(pred)
                    pred
                  },
                  
                  prediction.fold.prob = function(ensemble, x){
                    #Place all predictions for a sample in each row of a data set
                    pred <- as.data.frame(lapply(ensemble, function(e) prediction.prob(e,x)))
                    pred <- rowSums(pred)/sum(unlist(lapply(ensemble,length)))
                    pred
                  },
                  
                  
                  
                  
                  
                  train_IPIP_seq = function(  train.set, val.test){
                    
                    nmin = sum(train.set[[self$OUTPUT_VAR]] == self$OUTPUT_MIN)
                    nmaj = sum(train.set[[self$OUTPUT_VAR]] == self$OUTPUT_MAJ)
                    
                    
                    
                    majSubSize <-  round(self$np * self$prop_maj /(1-self$prop_maj ))
                    minSubSize <-  self$np
                    #Icludes from 1 to p all ensembles of the set
                    dfs <- list()
                    
                    minoritary <- subset(train.set, train.set[[self$OUTPUT_VAR]] == self$OUTPUT_MIN)
                    majoritary<- subset(train.set, train.set[[self$OUTPUT_VAR]] == self$OUTPUT_MAJ)
                    
                    for(k in 1:self$p){
                      id.minoritary <- sample(x = 1:nmin, size = minSubSize) #Index for minoritary class for each subset
                      id.majoritary <- sample(x= 1:nmaj, size = majSubSize) #Indexes for majoritary class for each subset
                      
                      dfs[[k]] <- rbind(minoritary[id.minoritary,],majoritary[id.majoritary,])
                    }
                    
                    
                    E <- list() #Final model (ensemble of ensembles)
                    
                    
                    for(k in 1:self$p){
                      Ek <- list() # k-esim ensemble
                      i <- 0 #Counter for number of tries of enlarging the ensemble
                      metric.ensemble  = 0 #Variable for storing the metric of each ensemble accumulated.
                      #Balanced partition
                      
                      df <- dfs[[k]]
                      model_i = 1
                      
                      #We select the elements for the training
                      majoritary <- which(df[[self$OUTPUT_VAR]] == self$OUTPUT_MAJ)
                      minoritary <- which(df[[self$OUTPUT_VAR]] == self$OUTPUT_MIN)
                      
                      
                      while(length(Ek$ensemble)<=self$b && i<self$mt(self$b,length(Ek$ensemble))){
                        #We use replace TRUE to make some randomness over the seed learners
                        ind.train <- c(
                          sample(majoritary, size = majSubSize, replace = TRUE),
                          sample(minoritary, size = minSubSize, replace = TRUE) 
                        )
                        
                        
                        #We train with the models in configuration with a sequential approach
                        
                        model <- self$conf[[length(Ek$ensemble)+1]](df[ind.train,], self$metrics, self$OUTPUT_VAR, self$TRAIN_METRIC)
                        
                        
                        if (length(Ek$ensemble)==0){
                          u <- -Inf;
                          names(u) <- self$SELECT_METRIC;
                          metrics.ensemble <- u
                        }
                        
                        new_Ek <- self$prediction.add(Ek, model, val.test = val.test[colnames(val.test)!=self$OUTPUT_VAR])
                        pred = self$prediction.predicted(new_Ek)
                        metrics.ensemble.new <- self$metrics(data.frame(
                          obs = val.test[[self$OUTPUT_VAR]],
                          pred= as.factor(pred)
                        ))
                        #We check if the new model changes the result. If it does not we start trying again
                        if(metrics.ensemble.new[self$SELECT_METRIC] <= metrics.ensemble[self$SELECT_METRIC]){ 
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
                  },
                  
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
                  },
                  
                  
                  train_IPIP_exhaustive = function(train.set, val.test){
                    
                    
                    
                    nmin = sum(train.set[[self$OUTPUT_VAR]] == self$OUTPUT_MIN)
                    nmaj = sum(train.set[[self$OUTPUT_VAR]] == self$OUTPUT_MAJ)
                    
                    
                    
                    majSubSize <-  round(self$np*self$prop_maj/(1-self$prop_maj))
                    minSubSize <-  self$np
                    
                    #Incluye en cada posicion los valores de los elementos de dicha particion, de 1 a p
                    dfs <- list()
                    
                    minoritary <- subset(train.set, train.set[[self$OUTPUT_VAR]] == self$OUTPUT_MIN)
                    majoritary <- subset(train.set, train.set[[self$OUTPUT_VAR]] == self$OUTPUT_MAJ)
                    
                    for(k in 1:self$p){
                      id.minoritary <- sample(x = 1:nmin, size = minSubSize) #Index for minoritary class for each subset
                      id.majoritary <- sample(x= 1:nmaj, size = majSubSize) #Indexes for majoritary class for each subset
                      
                      dfs[[k]] <- rbind(minoritary[id.minoritary,],majoritary[id.majoritary,])
                    }
                    
                    
                    
                    E <- list() #Final model (ensemble of ensembles)
                    
                    
                    for(k in 1:self$p){
                      
                      cat(sprintf("Ensemble number %d of %d\n", k,p))
                      Ek <- list() # k-esim ensemble
                      init = Sys.time()
                      
                      #Balanced partition
                      
                      df <- dfs[[k]]
                      majoritary <- which(df[[self$OUTPUT_VAR]] == self$OUTPUT_MAJ)
                      minoritary <- which(df[[self$OUTPUT_VAR]] == self$OUTPUT_MIN)
                      ind.train <- c(
                        sample(majoritary, size = majSubSize, replace = TRUE),
                        sample(minoritary, size = minSubSize, replace = TRUE) 
                      )
                      train.test <- df[ind.train,]
                      Ek <- list()
                      for(conf in 1:length(self$conf)){
                        
                        Ek[[length(Ek)+1]] <- self$conf[[conf]](train.test, self$metrics, self$OUTPUT_VAR, self$TRAIN_METRIC)
                        
                      }
                      
                      all_sets <- self$powerset(Ek, as.data.frame(val.test[, -which(names(val.test) == self$OUTPUT_VAR)]))
                      
                      ##If we want to force bigger sets before smaller ones
                      #lengths <- sapply(all_sets, function(x) length(x))
                      # Reorder the list based on the lengths
                      # all_sets <- all_sets[order(lengths, decreasing = TRUE)]
                      
                      
                      
                      max_metric = function(model){
                        pred = self$prediction.predicted(model)
                        met <- metrics(data.frame(
                          obs = val.test[[self$OUTPUT_VAR]],
                          pred = as.factor(pred)
                        ))[self$SELECT_METRIC]
                        return(met)
                      }
                      
                      
                      
                      #Run the max metric check in parallel mode
                      max_metric_list <- lapply(all_sets, max_metric)
                      
                      max_metric_value <- max(unlist(max_metric_list))
                      max_set <- all_sets[[which.max(unlist(max_metric_list))]]$ensemble
                      cat(sprintf("Max metric %s of ensemble is %f with length %d\n", self$SELECT_METRIC, max_metric_value, length(max_set)))
                      E[[length(E)+1]] <- max_set
                      
                    }
                    return(E);
                  }))
                  