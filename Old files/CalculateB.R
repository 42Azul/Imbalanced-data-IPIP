
calculate_b <- function( train.set,min_str, may_str, prob_codo=0.75, alpha_b= .01) 
{
  nmin = sum(train.set[[OUTPUT]] == min_str)
  nmay = sum(train.set[[OUTPUT]] == may_str)
  
  
  
  np <- ceiling(nmin*prob_codo)
  b <- ceiling(log(alpha_b)/(log(1-1/nmin)*np))
  return(b)
}
  