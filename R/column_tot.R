#'Adjust column margins
#'
#'This function helps to prepare column margins of the survey data as required
#' to use the function \code{spree}.
#' @param data_true input survey data with to be adjusted.
#' @param domains a character string containing the name of a variable that
#' indicates domains in the \code{data_true}.
#' @param total_true a data frame with sample data.
#' @return survey data set ready to use in \code{spree}
#' @export
#'

column_tot<-function(data_true, domains, total_true){
  data_true<-as.data.frame(select(data_true, -one_of(domains)))

  props<-prop.table(colSums(data_true, na.rm = T))
  column_total <- NULL
  for (i in 1:ncol(data_true)) {
    column_total <- rbind(column_total,as.vector(c(round((total_true)*props[i]))))

}

  return(column_total=c(column_total))
}

