% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robomit_functions.R
\name{o_beta_rsq_viz}
\alias{o_beta_rsq_viz}
\title{Visualization of beta*s over a range of maximum R-squares}
\usage{
o_beta_rsq_viz(y, x, con, m = "none", w = NULL, id = "none", time = "none", delta = 1,
type, data)
}
\arguments{
\item{y}{Name of the dependent variable (as string).}

\item{x}{Name of the independent treatment variable (i.e., variable of interest; as string).}

\item{con}{Name of related control variables. Provided as string in the format: "w + z +...".}

\item{m}{Name of unrelated control variables (m; see Oster 2019; as string; default is m = "none").}

\item{w}{weights (only for weighted estimations). Warning: For weighted panel models R can report different R-square than Stata, leading deviation between R and Stata results.}

\item{id}{Name of the individual id variable (e.g. firm or farm; as string). Only applicable for fixed effect panel models.}

\item{time}{Name of the time id variable (e.g. year or month; as string). Only applicable for fixed effect panel models.}

\item{delta}{delta for which beta*s should be estimated (default is delta = 1).}

\item{type}{Model type (either \emph{lm} or \emph{plm}; as string).}

\item{data}{Dataset.}
}
\value{
Returns ggplot2 object, which depicts beta*s over a range of maximum R-squares.
}
\description{
Estimates and visualizes beta*s, i.e., the bias-adjusted treatment effects (or correlations) (following Oster 2019) over a range of maximum R-squares.
}
\details{
Estimates and visualizes beta*s, i.e., the bias-adjusted treatment effects (or correlations) (following Oster 2019) over a range of maximum R-squares. The range of maximum R-squares starts from the R-square of the controlled model rounded up to the next 1/100 to 1. The function supports linear cross-sectional (see \emph{lm} objects in R) and fixed effect panel (see \emph{plm} objects in R) models.
}
\examples{
# load data, e.g. the in-build mtcars dataset
data("mtcars")
data_oster <- mtcars

# preview of data
head(data_oster)

# load robomit
require(robomit)

# estimate and visualize beta*s over a range of maximum R-squares
o_beta_rsq_viz(y = "mpg",            # dependent variable
               x = "wt",             # independent treatment variable
               con = "hp + qsec",    # related control variables
               delta = 1,            # delta
               type = "lm",          # model type
               data = data_oster)    # dataset
}
\references{
Oster, E. (2019). Unobservable Selection and Coefficient Stability: Theory and Evidence. Journal of Business & Economic Statistics, 37, 187-204.
}
