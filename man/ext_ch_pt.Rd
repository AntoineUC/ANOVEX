\name{ext_ch_pt}
\alias{ext_ch_pt}
\title{Extreme change-point detection}
\usage{
ext_ch_pt(Y,X=1:length(Y),kn=trunc(length(Y)/20),method="ANOVEX",grid,plot=FALSE)
}
\description{
Returns the location of a change-point in the extreme value properties of the underlying distribution for a heavy-tailed data series.
}
\arguments{
  \item{Y}{
Data vector. May be a time series.
}
  \item{X}{
Covariate vector. If Y is a time series, let X=1:length(Y).
}
  \item{kn}{
Number of upper order statistics. Default: 5\% of the data.
}
  \item{method}{
Change-point algorithm used. May be ANOVEX (default), KL, LR-GPD, LR-Pareto or Quantile.
}
  \item{grid}{
Grid on which the test statistics are computed. Default: seq(quantile(X,0.15,type=1),quantile(X,0.85,type=1),length=sum((as.numeric(names(table(X)))<=quantile(X,0.85))*(as.numeric(names(table(X)))>=quantile(X,0.15)))).
}
  \item{plot}{
Plot of the test statistics if TRUE. Default: FALSE.
}
}
\references{
Dupuis, D.J., Sun, Y., and Wang, H. J. (2015). Detecting change-points in extremes. \emph{Statistics and Its Interface}, 8: 19-31.

Girard, S., Opitz, T., and Usseglio-Carleve, A. (2024). ANOVEX: ANalysis Of Variability for heavy-tailed EXtremes, \emph{Electronic Journal of Statistics}, 18(2): 5258-5303.

Girard, S., Opitz, T., Usseglio-Carleve, A., and Yan, C. (2025). Changepoint identification in heavy-tailed distributions, preprint.
}
\author{
Antoine Usseglio-Carleve
}
