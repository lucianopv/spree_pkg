\name{census}
\alias{census}
\docType{data}
\title{Aggregated census data}
\usage{data("census")}
\format{
A data frame with 81 observations and 4 variables specifying the number of deprivations that a household has in the dimension ``health" of the multidimensional poverty index (MPI) of Costa Rica.
  \describe{
  \item{Canton}{factor; id of each domain (cantons).}
  \item{Zero}{numeric; total of households with zero deprivations.}
  \item{One}{numeric; total of households with one deprivations.}
  \item{Two}{numeric; total of households with two deprivations.}
  \item{Three_Four}{numeric; total of households with three or four deprivations.}
  }
}
\details{
The data set was created by the authors with information provided by The National Institute of Statistics and Censuses (INEC - Instituto Nacional de Estadística y Censos) of Costa Rica.
}
\usage{
census
}

\examples{
data(census)
str(census)
}
\keyword{datasets}
