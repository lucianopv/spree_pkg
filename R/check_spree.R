# This script contains the checks of arguments that have be done for the
# spree function.

spree_check1 <- function(population_data, sample_data){
  if(!all(colnames(population_data) %in% colnames(sample_data))){
    stop('Population and sample data must contain same variable names.')
  }
}

spree_check2 <- function(population_domains, sample_domains, population_data, sample_data){
if (!(population_domains %in% colnames(population_data))) {
  stop(paste0("The domain variable ", population_domains, " is not contained in population_data.
                Please provide valid variable name for population_domains."))
  if (!(sample_domains %in% colnames(sample_data))) {
  }
    stop(paste0("The domain variable ", sample_domains, " is not contained in sample_data.
                Please provide valid variable name for sample_domains."))
}
}




