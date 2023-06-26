source("IPIPClass.R")
myIPIPModel <- IPIP$new(OUTPUT_VAR, OUTPUT_MIN, OUTPUT_MAJ, seed_algorithms, metrics, metric_optimize)
myIPIPModel$exhaustive = TRUE
myIPIPModel$train(data.train, data.test)
