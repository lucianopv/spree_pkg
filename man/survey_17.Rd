\name{survey_17}
\alias{survey_17}
\docType{data}
\title{Aggregated survey data}
\usage{data("survey_17")}
\format{
A data frame with 80 observations and 4 variables specifying the number of deprivations that a household has in the dimension ``health" of the multidimensional poverty index (MPI) of Costa Rica.
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
survey_17
}

\examples{
data(survey_17)
str(survey_17)
}
\keyword{datasets}
