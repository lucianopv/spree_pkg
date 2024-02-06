\name{spree}
\alias{spree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Structure Preserving Estimation
%%  ~~function to do ... ~~
}
\description{
Function \code{spree} updates counts and proportions arranged by domains (areas) and a categorical variable (at least 2 categories). This function allows the implementation of the most populare SPREE-type estimations, based on \cite{Purcell, N. J. and L. Kish (1980)}, \cite{Zhang, L.-C. and R. L. Chambers (2004)}, and \cite{Luna-Hernández, A. (2016)}
}
\usage{
spree(population_domains,
                 sample_domains,
                 population_data,
                 sample_data,
                 row_margins = NULL,
                 col_margins = NULL,
                 type = c("SPREE", "GSPREE","MSPREE"),
                 B = 100,
                 method = c("ML", "IWLS"),
                 design_effect = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{population_data}{a data frame that with exactly same size as \code{sample_data} (A rows and J columns). It contains also \code{population_domains}, \code{row_margins}, and \code{col_margins}. No NAs are allowed. If this file contains zeros, they are treated as structural zeros, so they will remain as zeros in the output.}

\item{sample_data}{a data frame that with exactly same size as \code{population_data}. (A rows and J columns). It contains also \code{population_domains},\code{row_margins}, and \code{col_margins}.}

\item{population_domains}{a character string containing the name of a variable that indicates domains in the sample data. The variable can be numeric or a factor but needs to be of the same class as the variable named in \code{sample_domains}.}

\item{sample_domains}{a character string containing the name of a variable that indicates domains in the population data. The variable can be numeric or a factor but needs to be of the same class as the variable named in \code{population_domains}.}

\item{row_margins}{a numeric vector with A elements containing the ‘true’ domain sizes in t1. NULL is set as default and the sum by rows in sample_data is used.}

\item{col_margins}{numeric vector with J elements containing the ‘true’ totals for each category. The sum by columns in sample_data is used as default.}

\item{type}{a character string. SPREE version to implement: (i) ‘SPREE’; (ii) ‘GSPREE’; (iii) ‘MSPREE’. If this argument remains empty, the first option is selected. The option ‘MSPREE’ requires J >= 3. If J = 2, the ‘GSPREE’ will be applied.}

\item{B}{number of bootstraps used in the MSE estimation.}

\item{method}{a character string. Specify the method to obtain estimates of
the beta coefficient. This  argument is only required when \code{type} =
‘MSPREE’. The alternatives are (i) ‘ML’ and (ii) ‘IWLS’.}

\item{design_effect}{a character string. If type =‘MSPREE’ and method= ‘IWLS’
are selected, a design effect for the column totals can be provided. A
 vector containing 1’s is used as default.}
}
}
\details{
For GSPREE Poisson likelihood is used to obtain Maximum Likelihood estimates of the beta coefficient. This option can be also set for MSPREE when specifying ‘ML’ in \code{method}. As an alternative for the MSPSREE, Iterative Weighted Least Squares can be implemented. In this case, a variance estimate is required which is obtain using simple random samplig without replacement and multiplying it by a design effect.
For MSE estimation, a parametric bootstrap is implemented. Population compositions are repeatedly generated from a multinomial distribution.
}
\value{
An object of class "spree", with the updated population counts and respective uncertainty measures.
}
\references{
Luna-Hernandez (2016) Multivariate Structure Preserving Estimation
for Population Compositions. Ph.D. thesis, University of Southampton.\cr \cr
Purcell, N.J. & Kish, L. (1980) Postcensal estimates for local areas
(or domains). International Statistical Review/Revue Internationale de
Statistique, 48, 3– 18.\cr \cr
Zhang, L.C. & Chambers, R.L. (2004) Small area estimates for
cross-classifications. Journal of the Royal Statistical Society.
Series B: Statistical Methodology, 66, 479– 496. \cr \cr
}
\examples{
\donttest{
# Loading census and sample data
data("census")
data("survey_17")
data("P2017c")

# Prepare sample data
survey_17f <- prep_sample_data(population_domains = "Canton",
                                sample_domains = "Canton", population_data =
                                census, sample_data = survey_17)


# Adjust column margins, so sums by rows and columns are equal
ctotals_17 <- column_tot(data_true= survey_17f, domains = "Canton", total_true= sum(P2017c))

# SPREE
SPREE_17 <-spree(population_domains = "Canton",
                  sample_domains = "Canton", population_data = census,
                  sample_data = survey_17f, row_margins = P2017c,
                  col_margins = ctotals_17, type = "SPREE")

# GSPREE
GSPREE_17 <-spree(population_domains = "Canton",
                  sample_domains = "Canton", population_data = census,
                  sample_data = survey_17f, row_margins = P2017c,
                  col_margins = ctotals_17, type = "GSPREE")

# MSPREE with ML
MSPREE_ML_17 <-spree(population_domains = "Canton",
                  sample_domains = "Canton", population_data = census,
                  sample_data = survey_17f, row_margins = P2017c,
                  col_margins = ctotals_17, type = "MSPREE", method = "ML")

# MSPREE with IWLS
MSPREE_IWLS_17 <-spree(population_domains = "Canton",
                        sample_domains = "Canton", population_data = census,
                        sample_data = survey_17f, row_margins = P2017c,
                        col_margins = ctotals_17, type = "MSPREE",
                        method = "IWLS", design_effect = c(1,1,1,1))

# Summary of the updated point estimates
summary(MSPREE_IWLS_17$updated_point)

# Plots to compare methods
plots_spree<- compare_spree(SPREE_17,MSPREE_IWLS_17, GSPREE_17, MSPREE_ML_17)
}
