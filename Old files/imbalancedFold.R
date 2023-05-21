imbalancedFold <- function(data, n_folds, target, minority_class) {
  n_samples <- nrow(data)
  n_majority_total <- n_samples - sum(data[[target]] == minority_class)
  n_minority_total <- sum(data[[target]] == minority_class)
  
  n_minority_per_fold <- ceiling(n_minority_total / n_folds)
  n_majority_per_fold <- ceiling(n_samples / n_folds - n_minority_per_fold)

  fold_indices <- list()
  
  for (i in 1:n_folds) {
    # Select majority class samples
    majority_indices <- which(data[[target]] != minority_class)
    used_majority_indices <- unlist(fold_indices)
    available_majority_indices <- setdiff(majority_indices, used_majority_indices)
    
    n_available_majority <- length(available_majority_indices)
    n_majority_this_fold <- min(n_majority_per_fold, n_available_majority)
    selected_majority_indices <- sample(available_majority_indices, size = n_majority_this_fold)
    
    # Select minority class samples
    minority_indices <- which(data[[target]] == minority_class)
    used_minority_indices <- setdiff(used_majority_indices, available_majority_indices)
    available_minority_indices <- setdiff(minority_indices, used_minority_indices)
    n_available_minority <- length(available_minority_indices)
    n_minority_this_fold <- min(n_minority_per_fold, n_available_minority)
    selected_minority_indices <- sample(available_minority_indices, size = n_minority_this_fold)
    
    # Combine selected indices and add to fold_indices
    selected_indices <- c(selected_majority_indices, selected_minority_indices)

    fold_indices[[i]] <- selected_indices
  }
  
  return(fold_indices)
}
