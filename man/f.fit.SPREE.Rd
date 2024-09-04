\name{f.fit.SPREE}
\alias{f.fit.SPREE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Structure PREserving Estimation
%%  ~~function to do ... ~~
}
\description{
Function \code{f.fit.SPREE} *produces estimates of the coefficients of the chosen SPREE estimator. Implemented SPREE, GSPREE and MSPREE. Estimates can be obtained using Maximum Likelihood or IWLS.* Produce updates counts and proportions arranged by domains (areas) and a categorical variable (at least 2 categories). This function allows the implementation of the most popular SPREE-type estimations, based on \cite{Purcell, N. J. and L. Kish (1980)}, \cite{Zhang, L.-C. and R. L. Chambers (2004)}, and \cite{Luna-Hernández, A. (2016)}
}
\usage{
f.fit.SPREE(proxy,
            survey=NULL,
            estimator=c("spree","gspree", "mspree"),
            method=c("ml", "iwls"),
            iwls.opt = list(nsa.=NULL, ini.prop =NULL,                         deff=NULL))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{proxy}{a data frame that with exactly same size as \code{sample_data} (A rows and J columns). It contains also \code{population_domains}, \code{row_margins}, and \code{col_margins}. No NAs are allowed. If this file contains zeros, they are treated as structural zeros, so they will remain as zeros in the output.}
\item{survey} {This is the AxJ table obtained from the survey data. It is **not** required if estimator is "spree". If method is "ml", survey contains the observed survey counts in each category. If method is "iwls", survey contains the estimated population counts using the expansion factors. Please notice that when iwls is used, the survey table is internally transformed to the within-area proportions, so the expansion factors will only affect the estimate if the sampling design is *informative*, in the sense that the estimated within-area proportions are different from the observed within-area proportions in the sample.}
\item{estimator} {Can take the values "spree" (no update of structure),"gspree" (one coefficient) or "mspree" (J-1xJ-1 independent coefficients).}
\item{method} {Can take the values "ml" for Maximum Likelihood (Poisson or Multinomial are equivalent) or  "iwls" for Iterative Weighted Least Squares. ML estimates are obtained using the observed sample counts and assuming a Poisson/Multinomial distribution (both are equivalent). IWLS estimates are obtained using estimated counts obtained with sampling weights, and allow for the specification of category-specific design effects.}
\item{iwls.opt} {This is a list containing the additional elements required for method "iwls". Has components nsa. (vector of realized area sample sizes), ini.prop (AxJ matrix of within-area proportions which will be used to initialize the algorithm), and deff (vector of size J containing column-specific design effects. If missing, it is set to 1).}
%%     ~~Describe \code{x} here~~
}
\details{
- The tables proxy and survey only contain the relevant counts. No variable identifying the domains is included in either table.

- Both tables contain the exact same set of areas and the exact same categories, and these are in the same order in both tables. Notice that in the case where the set of areas differ between both sources (for instance because the survey doesn't produce estimates for all areas of interest), it is assumed that the user has filtered the common set of areas prior to using this function. I will ensure this by using a merge of the two sources as a check.

- No zeroes are allowed in the proxy table. No NA or negative values are allowed on either the proxy or the survey table. The survey table *can* contain zeroes but a large number of those in a given category may compromise the quality of the estimates. Zeroes in the proxy table can be substituted by a small count (e.g., 3) prior to using this function.

- If IWLS is used, arguments ini.prop and nsa. must be provided on iwls.opt. Notice that the area sample sizes nsa. are required in this case as survey contains expanded (population-level) counts.
Please use functions f.fit.SPREE and f.predict.SPREE to obtain a "spree" estimate of the within-area proportions and provide it as ini.prop.

- If IWLS is used, the user can provide a vector of column-specific design effects. If not provided, a deff of 1 is assumed for each column. If all columns have the same deff, IWLS is analogous to the fitting of an overdispersed Poisson/Multinomial model with dispersion=deff. The model coefficients, as well as the predicted counts and proportions will be the same as those assuming deff=1 but their estimated variances may be larger.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class "f.fit.SPREE", with:

\item{coeff} {Estimated coefficients. When estimator is "spree", it takes the value 1. When estimator is "gspree" is one constant. When estimator is "mspree" it is a list with two matrices of coefficients: beta (JxJ matrix of coefficients summing 1 by row and column) and Bp (reparametrization as in eq. 2.19, where the diagonal elements are free and the off-diagonal sum zero by row and column).}

\item{var.coeff} {Estimated variance of the coefficients.}
\item{proxy} {Copy of the proxy argument.}
\item{estimator} {Copy of the estimator argument.}
\item{proxy.alpha.aj} {Log-linear representation of the proxy table.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
\author{
Luna-Hernández A.
Arias-Salazar A.
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

#SPREE doesn't require the estimation of a coefficient. However, to make #the treatment analogous to that of GSPREE and MSPREE, we also require a #f.fit.SPREE step. method has no effect in this case.

res1 <- f.fit.SPREE(proxy=old.proxy, survey=survey.data,
                    estimator="spree", method="ml")

res1

#SPREE estimates for the same proxy as in res. Only row margins are used. }
}
