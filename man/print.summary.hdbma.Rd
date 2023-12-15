\name{print.summary.hdbma}
\alias{print.summary.hdbma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print
}
\description{
print the results from summary.hdbma.
}
\usage{
\method{print}{summary.hdbma}(x, ..., digit = x$digit, method = x$method, RE = x$RE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
the results from summary.hdbma.
}
  \item{\dots}{
further arguments passed to or from other methods.
}
  \item{digit}{
the number of digit to be shown.
}
  \item{method}{
the method to be shown. See the reference.
}
  \item{RE}{
if true, present results for relative effects. Default is True.
}
}
\value{
No return value, called for to print results and plots.
}
\references{
Yu, Q., Hagan, J., Wu, X., Richmond-Bryant, J., and Li, B., 2023, High-Dimensional Bayesian Mediation Analysis with Adaptive Laplace Priors. Submitted.
}
\author{
Qingzhao Yu and Bin Li
}
\examples{
#Examples at summary.hdbma.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
