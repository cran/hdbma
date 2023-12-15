\name{summary.hdbma}
\alias{summary.hdbma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary for hdbma results
}
\description{
This function summarize the results from hdbma objects with estimations, standard errors, and confidence intervals.
}
\usage{
\method{summary}{hdbma}(object, ..., plot = TRUE, RE = TRUE,
quant = c(0.025, 0.25, 0.5, 0.75, 0.975), digit = 4, method = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
the hdbma object from the hdbma function.
}
  \item{\dots}{
further arguments passed to or from other methods.
}
  \item{plot}{
if true, plot the estimation summaries. Default is True.
}
  \item{RE}{
if true, present results for relative effects. Default is True.
}
  \item{quant}{
the  set of quantiles to be shown in the sumamry results. Default is c(0.025, 0.25, 0.5, 0.75, 0.975).
}
  \item{digit}{
the number of digit to be shown.
}
  \item{method}{
the method to be shown. See the reference.
}
}
\value{
 \item{result1,\ldots,result4 }{the inference results for estimated mediation effects from methods 1 to 4.}
 \item{result1.re,\ldots,result4.re }{the inference results for estimated relative effects from methods 1 to 4.}
%% ...
}
\references{
Yu, Q., Hagan, J., Wu, X., Richmond-Bryant, J., and Li, B., 2023, High-Dimensional Bayesian Mediation Analysis with Adaptive Laplace Priors. Submitted.
}
\author{
Qingzhao Yu and Bin Li
}
\examples{
data("weight_behavior")
#for all the examples, needs to remove or increase n.iter and n.burnin
\donttest{
#binary predictor
test.b.c<- hdbma(pred=weight_behavior[,3], m=weight_behavior[,c(4,14,12,13)],
                 y=weight_behavior[,1],n.iter=10,n.burnin = 1)
summary(test.b.c)

##use covariate for y
test.b.c.2<- hdbma(pred=weight_behavior[,3], m=weight_behavior[,12:14],
                     y=weight_behavior[,1],cova=weight_behavior[,2],n.iter=10,n.burnin = 1)
summary(test.b.c.2)

#categorical predictor
test.ca.c<- hdbma(pred=weight_behavior[,4], m=weight_behavior[,12:14],
                     y=weight_behavior[,1],n.iter=10,n.burnin = 1)
summary(test.ca.c)

#use covariate for mediators
test.b.c.3<- hdbma(pred=weight_behavior[,3], m=weight_behavior[,c(9,12:14)],
                       y=weight_behavior[,1],mcov=weight_behavior[,c(2,5)],
                       mclist = list(1,2),n.iter=10,n.burnin = 1)
summary(test.b.c.3)

#use continuous predictor
test.c.c<- hdbma(pred=weight_behavior[,2], m=weight_behavior[,12:14],
                       y=weight_behavior[,1],n.iter=10,n.burnin = 1)
summary(test.c.c,method=3)

#use transfered continuous predictor
test.c.c.2<- hdbma(pred=weight_behavior[,2], m=weight_behavior[,12:14],
                     y=weight_behavior[,1],fpy=list(1,c("x","x^2")),
                     n.iter=10,n.burnin = 1)
summary(test.c.c.2,method=1)

#multiple predictors
test.m.c<- hdbma(pred=weight_behavior[,2:4], m=weight_behavior[,12:14],
                 y=weight_behavior[,1],n.iter=10,n.burnin = 1)
summary(test.m.c,RE=FALSE)

##binary outcome
test.m.b<- hdbma(pred=weight_behavior[,2:4], m=weight_behavior[,12:14],
                     y=weight_behavior[,15],cova=weight_behavior[,5],
                     n.iter=10,n.burnin = 1)
summary(test.m.b,method=2)

##categorical outcome
weight_behavior[,14]=as.factor(weight_behavior[,14])
test.m.c<- hdbma(pred=weight_behavior[,2:4], m=weight_behavior[,12:13],
                     y=weight_behavior[,14],cova=weight_behavior[,5],
                     n.iter=10,n.burnin = 1)
summary(test.m.c,method=2)
summary(test.m.c,method=1)

##time-to-event outcome
##Surv class outcome (survival analysis)
#use a simulation
set.seed(1)
N=100

alpha=0.5
x=rnorm(N,0,1)
x=ifelse(x>0,1,0)
e1=rnorm(N,0,1)
M=alpha*x+e1
lambda=0.01
rho=1
beta=1.2
c=-1
rateC=0.001
v=runif(n=N)
Tlat =(- log(v) / (lambda * exp(c*x+M*beta)))^(1 / rho)
C=rexp(n=N, rate=rateC)
time=pmin(Tlat, C)
status <- as.numeric(Tlat <= C)

test.m.t.1<- hdbma(pred=x, m=M,y=Surv(time,status),inits=function(){
  list(r=1,lambda=0.01)
},n.iter=10,n.burnin = 1)
summary(test.m.t.1,RE=FALSE)
}
}
