\name{ANOVEX_test}
\alias{ANOVEX_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ANalysis Of Variability for heavy-tailed EXtremes (ANOVEX) test
}
\description{
Test for tail equality between samples.
}
\usage{
ANOVEX_test(X,kn=trunc(min(sapply(X,length))/10),L=2,alpha=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{List of J samples. J must be greater than 2.}
  \item{kn}{Number of upper order statistics. Default: 10\% of the smallest sample size.}
  \item{L}{Number of extreme quantiles. Must be greater than 2 (default). This argument may be omitted if alpha is given instead.}
  \item{alpha}{Levels of the L extreme quantiles (between 0 and 1). Default: 1-(1:L)/n, where n is the average sample size.}
}
\details{
The value kn is denoted \deqn{n(1-\beta_n)} in the paper given as reference.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{statistic}{Value of the test statistic.}
  \item{p.value}{p-value of the test.}
  \item{method}{Name of the test.}
  \item{data.name}{Name of the data.}
%% ...
}
\references{
Girard, S., Opitz, T., and Usseglio-Carleve, A. (2024). ANOVEX: ANalysis Of Variability for heavy-tailed EXtremes, \emph{Electronic Journal of Statistics}, 18(2): 5258-5303.
}
\author{
Antoine Usseglio-Carleve
}
