\name{f.predict.SPREE}
\alias{f.predict.SPREE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Predict SPREE
%%  ~~function to do ... ~~
}
\description{
Contains the function f.predict.SPREE, which produces the point estimates of the tables of counts and within-area proportions.
}
\usage{
f.predict.SPREE(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{fitted} {This is an output from f.fit.SPREE. Brings information of the type of estimator implemented ("spree","gspree","mspree") and estimated coefficients, among others.}
\item{new.proxy} {Optional argument. If provided, estimates will be obtained updating the association structure of this table with the coefficients estimated in the previous step. This is useful when the set of areas used for fitting the model is different from the set of areas for which predictions are required. For instance, when there are out-of-sample areas (those are eliminated from the estimation step but it would be possible to obtain predictions if there is proxy information for them). If NULL, estimates will be obtained updating the association structure of the proxy table provided in f.fit.SPREE.}
\item{row.margin} {Required. This is a vector containing the totals for the same set of areas (and in the same order) as those in the proxy data (either the proxy provided in f.fit.SPREE or new.proxy). Estimates will satisfy this row margin.}
\item{col.margin} {Optional. This is a vector of J columns containing the totals aggregating over all areas. If provided, estimates will satisfy this column margin.}

}
\details{

- fitted must be  a valid output from f.fit.SPREE.

- if new.proxy is provided, it will only contain the relevant counts (no domain identifier), with the columns representing the same categories and in the same order as in the proxy and survey data used in argument fitted. The set of areas may differ from that used in the fitted argument. Cannot contain zeroes, NAs or negative values. Zeroes in the proxy table can be substituted by a small count (e.g., 3) prior to using this function.

- row.margin and col.margin are vectors of a dimension compatible with the proxy that will be used to obtain estimates. Neither can contain zeroes, NAs or negative values.
- the sum of row.margin and col.margin (if provided) must be exactly the same. }
\value{
The output is a list containing two tables of point estimates: counts (in the total scale) and proportions (in the within-area proportion scale).
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\donttest{
# Loading proxz and sample data
data("census")
data("survey_17")

# Substituting the 0's in the proxy by 3.

censusa <- census
censusa[censusa==0] <- 3
census <- censusa
rm(censusa)

# Increasing the size of the counts in the last two columns to avoid issues with MSE estimation later on.

survey_12$Two <- survey_12$Two*10
survey_12$Three_Four <- survey_12$Three_Four*10


# Filtering a common set of areas in proxy and survey and eliminating the area identifier:


old.proxy <- census[1:70,-1]
survey.data <- survey_12[1:70,-1]

# Fixing an issue with the margins, as they don't sum exactly the same amount.

data("P2012c")

column_tot_new<-function(data_true, domains, total_true){
  data_true<-as.data.frame(select(data_true, -one_of(domains)))

  props<-prop.table(colSums(data_true, na.rm = T))
  column_total <- NULL
  for (i in 1:ncol(data_true)) {
    column_total <- rbind(column_total,as.vector(c((total_true)*props[i])))
  }
  return(column_total=c(column_total))
}

row.mar <-P2012c[1:70]
col.mar <- column_tot_new(data_true= survey_12, domains = "Canton", total_true= sum(P2012c[1:70]))
ctotals_12 <- column_tot_new(data_true= survey_12, domains = "Canton", total_true= sum(P2012c))



# SPREE estimates

# SPREE doesn't require the estimation of a coefficient. However, to make the treatment analogous to that of GSPREE and MSPREE, we also require a f.fit.SPREE step. method has no effect in this case.

res1 <- f.fit.SPREE(proxy=old.proxy, survey=survey.data,
                    estimator="spree", method="ml")

res1

# SPREE estimates for the same proxy as in res. Only row margins are used.

res1.1 <- f.predict.SPREE(fitted=res1, new.proxy=NULL, row.margin=row.mar, col.margin=NULL)

res1.1

# SPREE estimates for the same proxy as in res. Row and column margins are used.


res1.2 <- f.predict.SPREE(fitted=res1, new.proxy=NULL, row.margin=row.mar, col.margin=col.mar)

res1.2

# SPREE estimates for the whole set of areas using the new.proxy argument. Row and column margins are used.


res1.3 <- f.predict.SPREE(fitted=res1, new.proxy=census[,-1], row.margin=P2012c, col.margin=ctotals_12)

res1.3

}
}
