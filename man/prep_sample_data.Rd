\name{prep_sample_data}
\alias{prep_sample_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prepare sample data
%%  ~~function to do ... ~~
}
\description{
This function helps to prepare the survey data set as required to use the
function \code{spree}. For example, the sum of rows and columns must be equal.
}
\usage{
prep_sample_data(population_domains,
                 sample_domains,
                 population_data,
                 sample_data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{population_domains}{population domains.}

\item{sample_domains}{sample domains.}

\item{population_data}{a data frame with population data.}

\item{population_data}{a data frame with sample data.}
}
\value{
survey data set ready to use in \code{spree}.
}
\examples{
# Example 1: Check and adjust sample data as census data.
data("census")
data("survey_17")

survey_17f <- prep_sample_data(population_domains = "Canton",
sample_domains = "Canton", population_data =
census, sample_data = survey_17)
}

