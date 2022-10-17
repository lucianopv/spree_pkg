\name{compare_spree}
\alias{compare_spree}
\title{Visual comparison between SPREE-type estimators
%%  ~~function to do ... ~~
}
\description{
This functions helps to compare visually the different estimators in terms of the point and uncertainty estimations. The instruction requires a minimum of two lists which are the outputs from using the function spree.
}
\usage{
compare_spree(spree_type1, spree_type2, spree_type3=NULL,spree_type4=NULL)
}

\examples{
plots_spree<- compare_spree(SPREE_17,MSPREE_IWLS_17)
}
