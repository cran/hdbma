\name{hdbma}
\alias{hdbma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
High-Dimensional Bayesian Mediation Analysis
}
\description{
We use the adaptive lasso priors for the Bayesian mediation analysis. Significant exposure variables, mediators are identified and their effects infered.
}
\usage{
hdbma(pred, m, y, refy = rep(NA, ncol(data.frame(y))),
predref = rep(NA, ncol(data.frame(pred))), fpy = NULL,
deltap = rep(0.001, ncol(data.frame(pred))), fmy = NULL,
deltam = rep(0.001, ncol(data.frame(m))), fpm = NULL,
mref = rep(NA, ncol(data.frame(m))), cova = NULL, mcov = NULL, mclist = NULL,
inits = NULL, n.chains = 1, n.iter = 1100, n.burnin = 100, n.thin = 1,
mucv = NULL, Omegacv = NULL, mu0.1 = NULL, Omega0.1 = NULL, mu1.1 = NULL,
Omega1.1 = NULL, mu0.a = NULL, Omega0.a = NULL, mu1.a = NULL, Omega1.a = NULL,
mu0.b = NULL, Omega0.b = NULL, mu1.b = NULL, Omega1.b = NULL, mu0.c = NULL,
Omega0.c = NULL, mu1.c = NULL, Omega1.c = NULL, preci = 1e-06, tmax = Inf,
multi = NULL, filename = NULL, deltax = 1, r1 = 1, partial = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pred}{
the vector/matrix of exposure(s)/predictor(s).
}
  \item{m}{
a data frame contains all potential mediators and covariates.
}
  \item{y}{
the vector/matrix of outcome(s).
}
  \item{refy}{
the reference group of y if the outcome is categorical.
}
  \item{predref}{
the reference group of pred if the exposure/predictor is categorical.
}
  \item{fpy}{
the transformation function of predictor(s) (pred) to explain y. [[1]] list all continuous predictors to be transformed, then following items list the transformation functions for each predictor in list [[1]] in that order.
}
  \item{deltap}{
the vector of changing amount in predictors.
}
  \item{fmy}{
the transformation functions of mediators (m) to explain y, [[1]] list all continuous mediators in m to be transformed, then following items list the transformation functions for each mediator in list [[1]] in that order.
}
  \item{deltam}{
the vector of changing amount in mediators.
}
  \item{fpm}{
the transformation functions of predictors (pred) to explain mediators (m), [[1]] is a matrix, the first column indicator the mediators to be explained, the second column are the continuous predictors to be transformed; then transformation functions are listed in the following items by the row order of [[1]].
}
  \item{mref}{
the reference group of m if any of them is categorical. By default, the reference group is the first one in alphebetic order.
}
  \item{cova}{
the covariates for the outcome.
}
  \item{mcov}{
the data frame with all covariates for mediators
}
  \item{mclist}{
the list of all covariates for mediators. If mclist is NULL but mcov is not, use all covariates in mcov for all mediators. Otherwise the first item of mclist lists all mediators that are using different covariates, the following items give the columns of covariates in mcov for the mediators in order of mclist[[1]]. Use NA is no covariates are to be used. If a mediator is not listed in mclist[[1]], use all covariates in mcov.
}
  \item{inits}{
to specify the starting values of parameters. Default is NULL. See R2jags:jags.
}
  \item{n.chains}{
number of Markov chains (default: 1). See R2jags:jags.
}
  \item{n.iter}{
number of total iterations per chain (including burn in; default: 1100). See R2jags:jags.
}
  \item{n.burnin}{
length of burn in, i.e. number of iterations to discard at the beginning. Default is 100. If n.burnin is 0, jags() will run 100 iterations for adaption. See R2jags:jags.
}
  \item{n.thin}{
thinning rate. Must be a positive integer. Default is 1. See R2jags:jags.
}
  \item{mucv}{
the prior mean for the variables in cova. Default is 0.
}
  \item{Omegacv}{
the prior precision for the variables in cova. Default is 0 preci.
}
  \item{mu0.1}{
the prior mean for the intercept in the prediction model for all mediators. Default is 0.
}
  \item{Omega0.1}{
the prior precision for the intercept in the prediction model for all mediators. Default is preci.
}
  \item{mu1.1}{
a vector of the size of meditators include the prior mean for the slope of the exposures in the prediction model for all mediators. Default is rep(0,P).
}
  \item{Omega1.1}{
the prior precision matrix (P*P) for the slope of the exposures in the prediction model for all mediators. Default is a diagoal matrix with preci.
}
  \item{mu0.a}{
the prior mean for the intercept in the prediction model for all continuous mediators. Default is rep(0,p1), p1 is the number of continuous mediators.
}
  \item{Omega0.a}{
the prior precision for the intercept in the prediction model for all continuous mediators. Default is diag(preci).
}
  \item{mu1.a}{
the prior mean for the slope of the exposures in the prediction model for all continuous mediators. Default is rep(0,p1), p1 is the number of continuous mediators.
}
  \item{Omega1.a}{
the prior precision for the slope of the exposures in the prediction model for all continuous mediators. Default is diag(preci).
}
  \item{mu0.b}{
the prior mean for the intercept in the prediction logit model for all binary mediators. Default is rep(0,p2), p2 is the number of bianry mediators.
}
  \item{Omega0.b}{
the prior precision matrix for the intercept in the prediction logit model for all binary mediators. Default is diag(preci).
}
  \item{mu1.b}{
the prior mean for the slope of exposure(s) in the prediction logit model for all binary mediators. Default is rep(0,p2).
}
  \item{Omega1.b}{
the prior precision matrix for the slope of exposure(s) in the prediction logit model for all binary mediators. Default is diag(preci).
}
  \item{mu0.c}{
the prior mean for the intercept in the prediction logit model for all categorical mediators. Default is 0 for all (array(0,p3,cat1,nmc), where p3 is the number of categorical mediators, cat1 is the maximum number of categories of all categorical mediators, and nmc is the number of mcov.
}
  \item{Omega0.c}{
the prior precision matrix for intercept in the prediction logit model for all categorical mediators. Default is preci at the diagnal matrix of nmc*nmc dimension.
}
  \item{mu1.c}{
the prior mean for the slope of exposures in the prediction logit model for all categorical mediators. Default is 0 for all (array(0,p3,cat1,c1), where p3 is the number of categorical mediators, cat1 is the maximum number of categories of all categorical mediators, and c1 is the number of exposures.
}
  \item{Omega1.c}{
the prior precision matrix for the slope of exposures in the prediction logit model for all categorical mediators. Default is preci at the diagnal matrix of c1*c1 dimension.
}
  \item{preci}{
the prior precision level. Default is 0.000001.
}
  \item{tmax}{
the maximum time to event for survival analysis. Default is Inf.
}
  \item{multi}{
in the survival analysis only. If true, calculate the multiplicative effect of survival time.
}
  \item{filename}{
the directory and filename for the bugs model. If is NULL, the function will generate the bugs model using default functions.
}
  \item{deltax}{
the change unit in the exposures to calculate the mediation effects. Default is 1.
}
  \item{r1}{
the penalty parameter is \deqn{\lambda^*=\frac{\lambda}{(\alpha\hat{\beta})^r}.}. Default is 1.
}
  \item{partial}{
if true, do the partial lasso, \deqn{\lambda^*=\frac{\lambda}{\alpha^r}.} Default is FALSE.
}
}
\details{
The function will automatically catch the types of the outcome and compile the bugs model. Results will be summarized use the summary function. Please see examples under summary.
}
\value{
A hdbma object is returned with the following items. Results are summarized using the summary.hdbma function.
  \item{aie1,ade1,ate1,\ldots,ate4}{
the average indirect effect (ie), direct effect (de) and total effect (te) from the four different methods. See the references.
}
  \item{sims.list}{
the simulation results from the bugs model.
}
  \item{data0}{
the organized data set that was analyzed.
}
  \item{omu3,omu4}{
the total effects from method 3 or 4 for survival outcome only.
}
}
\references{
Yu, Q., Hagan, J., Wu, X., Richmond-Bryant, J., and Li, B., 2023, High-Dimensional Bayesian Mediation Analysis with Adaptive Laplace Priors. Submitted.
}
\author{
Qingzhao Yu and Bin Li
}
\examples{
#Check summary.hdbma.
}
