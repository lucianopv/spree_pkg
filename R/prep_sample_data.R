#'Prepare sample data
#'
#'This function helps to prepare the survey data set as required to use the
#' function \code{spree}. For example, the sum of rows and columns must be equal.
#' @param population_domains population domains.
#' @param sample_domains sample domains.
#' @param population_data a data frame with population data.
#' @param sample_data a data frame with sample data.
#' @return survey data set ready to use in \code{spree}
#' @importFrom dplyr n_distinct
#' @export
NULL

prep_sample_data<-function(population_domains,
                           sample_domains,
                           population_data,
                           sample_data){

  population_data <- as.data.frame(population_data)
  sample_data <- as.data.frame(sample_data)

  #population_data<-population_data[order(population_domains),]
  #sample_data<-sample_data[order(sample_domains),]

  sample_domains_vec <- sample_data[, sample_domains]
  population_domains_vec <- population_data[, population_domains]

  obs_dom <- population_domains_vec %in% sample_domains_vec
  out_smp <- population_domains_vec[!(population_domains_vec %in% sample_domains_vec)]

  sample_data[, sample_domains] <- NULL
  population_data[, population_domains] <- NULL

  sample_data[, population_domains] <- population_domains_vec[obs_dom == TRUE]
  population_data[, population_domains] <- population_domains_vec
  population_data<-as.data.frame(population_data[,population_domains])
  colnames(population_data)<-population_domains
  data <- merge(sample_data, population_data, by = population_domains, all = TRUE)
  data[is.na(data)]<-1
  #data[data==0]<-1

  summary_info<-list(cat("Out of sample domains: \n"),
                     cat("Percentage: ", round((n_distinct(out_smp)/n_distinct(population_domains_vec))*100,2),"%","\n"),
                     cat("Number: ", n_distinct(out_smp), "\n"),
                     cat("Domains: ", out_smp, "\n")
                     )


  return(data=data)
}
