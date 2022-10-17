#' census data
#'
#' This data set contains aggregated information for 81 Cantons of Costa Rica from
#' the Population and Housing Census 2011 of Costa Rica.
#' The data set was created by the authors with information provided by
#' The National Institute of Statistics and Censuses (INEC - Instituto
#' Nacional de Estadística y Censos) of Costa Rica
#'
#' @format A data frame with 81 observations and 4 variables specifying the number
#' of deprivations that a household has in the dimension ``health" of the
#' multidimensional poverty index (MPI) if Costa Rica:
#' \describe{
#' \item{Canton}{factor; id of each domain (cantons).}
#' \item{Zero}{numeric; total of households with zero deprivations.}
#' \item{One}{numeric; total of households with one deprivations.}
#' \item{Two}{numeric; total of households with two deprivations.}
#' \item{Three_Four}{numeric; total of households with three or four deprivations.}

#' }
#' @docType data
"census"
