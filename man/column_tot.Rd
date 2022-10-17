\name{column_tot}
\alias{column_tot}
\title{Adjust the column margins of the sample data}
\usage{
column_tot(data_true, domains, total_true)
}
\arguments{
\item{data_true}{input survey data with to be adjusted.}

\item{domains}{a character string containing the name of a variable that indicates domains in the \code{data_true}.}

\item{total_true}{a data frame with sample data.}
}
\value{
survey data set ready to use in \code{spree}
}
\description{
This function helps to prepare column margins of the survey data as required to use the function \code{spree}.
}
