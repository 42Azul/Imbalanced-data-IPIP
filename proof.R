source("IPIPClass.R")
myIPIPModel <- IPIP$new(OUTPUT_VAR, OUTPUT_MIN, OUTPUT_MAJ, c(seed_algorithms, seed_algorithms), metrics, metric_optimize, exhaustive = FALSE)
myIPIPModel$train(data.train, data.test)
