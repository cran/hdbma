hdbma<-function(pred, m, y, refy = rep(NA, ncol(data.frame(y))),
                    predref = rep(NA, ncol(data.frame(pred))), fpy = NULL,
                    deltap = rep(0.001,ncol(data.frame(pred))), fmy = NULL,
                    deltam = rep(0.001,ncol(data.frame(m))), fpm = NULL,
                    mref = rep(NA, ncol(data.frame(m))),cova = NULL,
                    mcov = NULL, mclist=NULL, inits=NULL,
                    n.chains = 1, n.iter = 1100,n.burnin=100,n.thin = 1,
                    mucv=NULL, Omegacv=NULL,
                    mu0.1=NULL,  Omega0.1=NULL, mu1.1=NULL, Omega1.1=NULL,
                    mu0.a=NULL, Omega0.a=NULL, mu1.a=NULL, Omega1.a=NULL,
                    mu0.b=NULL, Omega0.b=NULL, mu1.b=NULL, Omega1.b=NULL,
                    mu0.c=NULL, Omega0.c=NULL, mu1.c=NULL, Omega1.c=NULL,
                    preci=0.000001,tmax=Inf,multi=NULL,filename=NULL,deltax=1,r1=1,
                    partial=FALSE)
{
data_org<-function(pred, m, y, refy = rep(NA, ncol(data.frame(y))),
                    predref = rep(NA, ncol(data.frame(pred))), fpy = NULL,
                    deltap = rep(0.001,ncol(data.frame(pred))), fmy = NULL,
                    deltam = rep(0.001,ncol(data.frame(m))), fpm = NULL,
                    mref = rep(NA, ncol(data.frame(m))),cova = NULL,
                    mcov = NULL, mclist=NULL)  #mcov is the data frame with all covariates for TVs, mind is the indicator for covariates
  #if mclist is null but not mcov, mcov is applied to all tvs.
{ns.dev<-function (x, df = NULL, knots = NULL, qnots=NULL,intercept = FALSE, Boundary.knots = range(x),derivs1=0)
{
  nx <- names(x)
  x <- as.vector(x)
  nax <- is.na(x)
  if (nas <- any(nax))
    x <- x[!nax]
  if (!missing(Boundary.knots)) {
    Boundary.knots <- sort(Boundary.knots)
    outside <- (ol <- x < Boundary.knots[1L]) | (or <- x >
                                                   Boundary.knots[2L])
  }
  else outside <- FALSE
  if (!is.null(df) && is.null(knots)) {
    nIknots <- df - 1L - intercept
    if (nIknots < 0L) {
      nIknots <- 0L
      warning(gettextf("'df' was too small; have used %d",
                       1L + intercept), domain = NA)
    }
    knots <- if (nIknots > 0L) {
      knots <- seq.int(0, 1, length.out = nIknots + 2L)[-c(1L,
                                                           nIknots + 2L)]
      stats::quantile(x[!outside], knots)
    }
  }
  else {if(is.null(df) && is.null(knots) && !is.null(qnots))
    knots<-quantile(x[!outside], qnots)
  nIknots <- length(knots)}
  Aknots <- sort(c(rep(Boundary.knots, 4L), knots))
  if (any(outside)) {
    basis <- array(0, c(length(x), nIknots + 4L))
    if (any(ol)) {
      k.pivot <- Boundary.knots[1L]
      xl <- cbind(1, x[ol] - k.pivot)
      tt <- splineDesign(Aknots, rep(k.pivot, 2L), 4, c(0,
                                                        1),derivs=rep(derivs1,2L))
      basis[ol, ] <- xl %*% tt
    }
    if (any(or)) {
      k.pivot <- Boundary.knots[2L]
      xr <- cbind(1, x[or] - k.pivot)
      tt <- splineDesign(Aknots, rep(k.pivot, 2L), 4, c(0,1),derivs=rep(derivs1,2L))
      basis[or, ] <- xr %*% tt
    }
    if (any(inside <- !outside))
      basis[inside, ] <- splineDesign(Aknots, x[inside],
                                      4,derivs=rep(derivs1,length(x[inside])))
  }
  else basis <- splineDesign(Aknots, x, 4,derivs=rep(derivs1,length(x)))
  const <- splineDesign(Aknots, Boundary.knots, 4, c(2, 2),derivs=rep(derivs1,length(Boundary.knots)))
  if (!intercept) {
    const <- const[, -1, drop = FALSE]
    basis <- basis[, -1, drop = FALSE]
  }
  qr.const <- qr(t(const))
  basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1L:2L),
                                                     drop = FALSE])
  n.col <- ncol(basis)
  if (nas) {
    nmat <- matrix(NA, length(nax), n.col)
    nmat[!nax, ] <- basis
    basis <- nmat
  }
  dimnames(basis) <- list(nx, 1L:n.col)
  a <- list(degree = 3L, knots = if (is.null(knots)) numeric() else knots,
            Boundary.knots = Boundary.knots, intercept = intercept)
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("ns", "basis", "matrix")
  basis
}



bs.dev<-function (x, df = NULL, knots = NULL, degree = 3, intercept = FALSE,
                  Boundary.knots = range(x),derivs1=0)
{
  nx <- names(x)
  x <- as.vector(x)
  nax <- is.na(x)
  if (nas <- any(nax))
    x <- x[!nax]
  if (!missing(Boundary.knots)) {
    Boundary.knots <- sort(Boundary.knots)
    outside <- (ol <- x < Boundary.knots[1L]) | (or <- x >
                                                   Boundary.knots[2L])
  }
  else outside <- FALSE
  ord <- 1L + (degree <- as.integer(degree))
  if (ord <= 1)
    stop("'degree' must be integer >= 1")
  if (!is.null(df) && is.null(knots)) {
    nIknots <- df - ord + (1L - intercept)
    if (nIknots < 0L) {
      nIknots <- 0L
      warning(gettextf("'df' was too small; have used %d",
                       ord - (1L - intercept)), domain = NA)
    }
    knots <- if (nIknots > 0L) {
      knots <- seq.int(from = 0, to = 1, length.out = nIknots +
                         2L)[-c(1L, nIknots + 2L)]
      stats::quantile(x[!outside], knots)
    }
  }
  Aknots <- sort(c(rep(Boundary.knots, ord), knots))
  if (any(outside)) {
    warning("some 'x' values beyond boundary knots may cause ill-conditioned bases")
    derivs <- 0:degree
    scalef <- gamma(1L:ord)
    basis <- array(0, c(length(x), length(Aknots) - degree -
                          1L))
    if (any(ol)) {
      k.pivot <- Boundary.knots[1L]
      xl <- cbind(1, outer(x[ol] - k.pivot, 1L:degree,
                           "^"))
      tt <- splineDesign(Aknots, rep(k.pivot, ord), ord,
                         derivs+derivs1)
      basis[ol, ] <- xl %*% (tt/scalef)
    }
    if (any(or)) {
      k.pivot <- Boundary.knots[2L]
      xr <- cbind(1, outer(x[or] - k.pivot, 1L:degree,
                           "^"))
      tt <- splineDesign(Aknots, rep(k.pivot, ord), ord,
                         derivs+derivs1)
      basis[or, ] <- xr %*% (tt/scalef)
    }
    if (any(inside <- !outside))
      basis[inside, ] <- splineDesign(Aknots, x[inside],
                                      ord,derivs=rep(derivs1,length(x[inside])))
  }
  else basis <- splineDesign(Aknots, x, ord, derivs=rep(derivs1,length(x)))
  if (!intercept)
    basis <- basis[, -1L, drop = FALSE]
  n.col <- ncol(basis)
  if (nas) {
    nmat <- matrix(NA, length(nax), n.col)
    nmat[!nax, ] <- basis
    basis <- nmat
  }
  dimnames(basis) <- list(nx, 1L:n.col)
  a <- list(degree = degree, knots = if (is.null(knots)) numeric(0L) else knots,
            Boundary.knots = Boundary.knots, intercept = intercept)
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("bs", "basis", "matrix")
  basis
}


x2fx<-function(x,func) #x is the list of original numerical vector, func is a vector of character functions.
{ # eg.  func <- c("x","x+1","x+2","x+3","log(x)")
  func.list <- list()
  #test.data <- matrix(data=rep(x,length(func)),length(x),length(func))
  #test.data <- data.frame(test.data)
  result<-NULL
  for(i in 1:length(func)){
    func.list[[i]] <- function(x){}
    body(func.list[[i]]) <- parse(text=func[i])
  }
  #result <- mapply(do.call,func.list,lapply(test.data,list))
  col_fun<-NULL
  z<-1
  for (i in 1:length(func.list))
  {res<-as.matrix(func.list[[i]](x))
  result<-cbind(result,res)
  col_fun<-cbind(col_fun,c(z,z+ncol(res)-1))
  z<-z+ncol(res)}
  list(values=as.matrix(result),col_fun=as.matrix(col_fun))
}

x2fdx<-function(x,func)  #x is the list of original numerical vector, func is a vector of character functions.
{ fdx<-NULL              # eg.  func <- c("x","x+1","x+2","x+3","log(x)")
for(i in 1:length(func)){
  if (length(grep("ifelse",func[i]))>0)
  {str<-unlist(strsplit(func[i],","))
  fun1<-D(parse(text=str[2]), "x")
  fun2<-D(parse(text=unlist(strsplit(str[3],")"))),"x")
  x1<-eval(fun1)
  x2<-eval(fun2)
  if(length(x1)==1)
    x1<-rep(x1,length(x))
  if(length(x2)==1)
    x2<-rep(x2,length(x))
  fun3<-paste(str[1],"x1,x2)",sep=",")
  fdx<-cbind(fdx,eval(parse(text=fun3)))
  }
  else if(length(grep("ns",func[i]))>0)
  {temp<-paste("ns.dev",substring(func[i],3,nchar(func[i])-1),",derivs1=1)",sep="")
  fdx<-cbind(fdx,eval(parse(text=temp)))}
  else if(length(grep("bs",func[i]))>0)
  {temp<-paste("bs.dev",substring(func[i],3,nchar(func[i])-1),",derivs1=1)",sep="")
  fdx<-cbind(fdx,eval(parse(text=temp)))}
  else{
    dx2x <- D(parse(text=func[i]), "x")
    temp<-eval(dx2x)
    if(length(temp)==1)
      fdx<-cbind(fdx,rep(temp,length(x)))
    else fdx<-cbind(fdx,temp)}
}
as.matrix(fdx)
}

bin_cat <- function(M2, M1, cat1, catref = 1)
  #turn a categorical variable to binary dummy variables
  #M2 is the original data frame
  #cat1 is the column number of the categorical variable in M2
  #catref is the reference group
{a <- factor(M2[, cat1])
b <- sort(unique(a[a != catref]))
d<-NULL
e<-rep(1,nrow(M2))
for(i in 1:length(b))
{d=cbind(d,ifelse(a==b[i],1,0))
e=ifelse(a==b[i],i+1,e)
}
d[is.na(M2[,cat1]),]=NA
e[is.na(M2[,cat1])]=NA
M2[,cat1]=e
xnames=colnames(M1)
M1=cbind(M1,d)
colnames(M1)=c(xnames,paste(colnames(M2)[cat1],b,sep="."))
list(M1=M1,M2=M2,cat=c(ncol(M1)-ncol(d)+1,ncol(M1)))
}

order_char<-function(char1,char2)  #find the position of char2 in char1
{a<-1:length(char1)
b<-NULL
for (i in 1:length(char2))
  b<-c(b,a[char1==char2[i]])
b
}

###start the main code
#clean up the outcomes:y_type=type of outcomes;
#y_type=2:binary, 3:category, 1:continuous, 4:time-to-event
#consider only 1 outcome for now

#consider only complete case
data.temp=cbind(pred,m,y)
if (!is.null(mcov) & !is.null(mclist))
  data.temp=cbind(pred,m,y,mcov)
if(!is.null(cova))
  data.temp=cbind(data.temp,cova)
choose.temp=complete.cases(data.temp)
if(ncol(data.frame(pred))==1)
  pred=pred[choose.temp]
else
  pred=pred[choose.temp,]
if(ncol(data.frame(y))==1)
  y=y[choose.temp]
else
  y=y[choose.temp,]
if(ncol(data.frame(m))==1)
  m=m[choose.temp]
else
  m=m[choose.temp,]
if(!is.null(cova)){
  if(ncol(data.frame(cova))==1)
    cova=cova[choose.temp]
  else
    cova=cova[choose.temp,]}
if(!is.null(mcov)){
  if(ncol(data.frame(mcov))==1)
    mcov=mcov[choose.temp]
  else
    mcov=mcov[choose.temp,]}

if (!is(y,"Surv"))
{if (nlevels(droplevels(as.factor(y))) == 2) {#binary
  y_type <- 2
  if (!is.na(refy))
    y <- ifelse(y == refy, 0, 1)
  else {
    refy <- levels(droplevels(as.factor(y)))[1]
    y <- ifelse(as.factor(y) == refy,0, 1)
  }
}
  else if (is.character(y) | is.factor(y)) {#categorical
    y_type <- 3
    y <- droplevels(y)
    if (is.na(refy))
      refy <- levels(as.factor(y))[1]
    a <- factor(y)
    b <- sort(unique(a[a != refy]))
    e <- rep(1,nrow(as.matrix(y)))
    for(j in 1:length(b))
      e=ifelse(a==b[j],j+1,e)
    e[is.na(y)]=NA
    y=factor(e)
  }
  else  #continuous
    y_type = 1}
else
{y_type<-4
y=cbind(y[,1],y[,2])}


#clean up m and pred
mnames <- colnames(m)
if (!is.null(cova)) {
  if (is.null(colnames(cova)))
    cova_names = paste("cova",1:ncol(data.frame(cova)))
  else
    cova_names = colnames(cova)
  cova=data.frame(cova)
  colnames(cova)=cova_names
}

if (!is.null(mcov)) {
  if (length(grep("for.m", names(mcov))) == 0)
    mcov_names = colnames(mcov)
  else mcov_names = colnames(mcov[[1]])
  mcov=data.frame(mcov)
}
pred_names = names(pred)

##prepare for the predictor(s)
pred1 <- data.frame(pred) #original format
pred1.0 <- NULL           #original format with binarized categorical predictor
pred1.0_names <-NULL
pred2 <- NULL             #all transformed
pred2_names = NULL
pred3 <- NULL             #transformed continuous predictors with pred+delta(pred)
pred3_names = NULL
pred.cont.der <- NULL     #derivative of the transformation function for cont predictors
if (is.null(pred_names))
  pred_names = paste("pred",1:ncol(pred1),sep='')
colnames(pred1) = pred_names
binpred = NULL   #binary predictor in pred2
catpred = NULL   #categorical predictor in pred2, each row is for one categorical predictor
contpred = NULL  #continuos predictor in pred2, each row is for a continuous predictor
binpred1 = NULL  #binary predictor in pred1
catpred1 = NULL  #categorical predictor in pred1
contpred1 = NULL #continuous predictor in pred1
binpred1.0 = NULL  #binary predictor in pred1.0
catpred1.0 = NULL  #categorical predictor in pred1.0
contpred1.0 = NULL #continuous predictor in pred1.0
contpred3 = NULL         #index for pred3 and pred.cont.dev
npred = ncol(pred1)
n1=nrow(pred1)

for (i in 1:npred)
  if (nlevels(droplevels(as.factor(pred1[,i]))) == 2) { #binary predictor
    if (!is.na(predref[i]))
    {pred2 <- cbind(pred2,ifelse(pred1[, i] == predref[i],0, 1))
    pred1.0 <- cbind(pred1.0, ifelse(pred1[, i] == predref[i],0, 1))
    #pred3 <- cbind(pred3,as.factor(ifelse(pred1[, i] == predref[i],0, 1)))
    pred1[,i] <- ifelse(pred1[, i] == predref[i],0, 1)}
    else {
      temp.pred <- as.factor(pred1[, i])
      pred2 <- cbind(pred2,ifelse(temp.pred == levels(droplevels(temp.pred))[1], 0, 1))
      pred1.0 <- cbind(pred1.0,ifelse(temp.pred == levels(droplevels(temp.pred))[1], 0, 1))
      #pred3 <- cbind(pred3,as.factor(ifelse(temp.pred == levels(droplevels(temp.pred))[1], 0, 1)))
      pred1[,i] <- ifelse(temp.pred == levels(droplevels(temp.pred))[1], 0, 1)
    }
    binpred1 = c(binpred1, i)
    pred2_names=c(pred2_names,pred_names[i])
    pred1.0_names=c(pred1.0_names,pred_names[i])
    binpred = c(binpred,ncol(pred2))
    binpred1.0 = c(binpred1.0,ncol(pred1.0))
  }
else if (is.character(pred1[, i]) | is.factor(pred1[, i])) { #category predictor
  pred1[, i] = droplevels(pred1[, i])
  if(!is.null(pred2))
    colnames(pred2)=pred2_names
  if(!is.null(pred1.0))
    colnames(pred1.0)=pred1.0_names
  #    catn = catn + 1
  if (!is.na(predref[i]))
  {pred.temp1 <- bin_cat(pred1, pred2, i, predref[i])
  pred.temp2 <- bin_cat(pred1, pred1.0, i, predref[i])}
  else
  {pred.temp1 <- bin_cat(pred1, pred2, i, levels(as.factor(pred1[,i]))[1])
  pred.temp2 <- bin_cat(pred1, pred1.0, i, levels(as.factor(pred1[,i]))[1])}
  pred2 = pred.temp1$M1
  pred1.0 = pred.temp2$M1
  pred1 = pred.temp1$M2
  #pred3 = cbind(pred3,pred.temp1$M1[,pred.temp1$cat[1]:pred.temp1$cat[2]])
  catpred1 = c(catpred1,i)
  catpred = rbind(catpred,pred.temp1$cat)
  catpred1.0 = rbind(catpred1.0,pred.temp2$cat)
  pred2_names = colnames(pred.temp1$M1)
  pred1.0_names = colnames(pred.temp2$M1)
}
else #consider the transformation of continuous x
{contpred1 = c(contpred1, i)
pred1.0=cbind(pred1.0,pred1[,i])
pred1.0_names<-c(pred1.0_names,pred_names[i])
contpred1.0= c(contpred1.0, ncol(pred1.0))
if(!(i %in% fpy[[1]]))  #fpy has the transformation functions for pred to y, [[1]] list all cont pred to be transformed
{pred2<-cbind(pred2,pred1[,i])
pred.cont.der<-cbind(pred.cont.der,rep(1,n1))
contpred3 = rbind(contpred3,c(ncol(pred.cont.der),ncol(pred.cont.der)))
pred3 = cbind(pred3,pred1[,i]+deltap[i]) #deltap is the changing amount for predictors
contpred = rbind(contpred,c(ncol(pred2),ncol(pred2)))
pred2_names=c(pred2_names,pred_names[i])
pred3_names=c(pred3_names,pred_names[i])
}
else if (i %in% fpy[[1]])
{p=match(i,fpy[[1]])
temp.pred=x2fx(pred1[,i],fpy[[1+p]])$values
pred2<-cbind(pred2,temp.pred)
pred.cont.der<-cbind(pred.cont.der,x2fdx(pred1[,i],fpy[[1+p]]))
if(is.null(pred3))
  contpred3 = rbind(contpred3,c(1,ncol(pred.cont.der)))
else
  contpred3 = rbind(contpred3,c(ncol(pred3)+1,ncol(pred.cont.der)))
pred3<-cbind(pred3,x2fx(pred1[,i]+deltap[i],fpy[[1+p]])$values)
l=length(fpy[[1+p]])
contpred=rbind(contpred,c(ncol(pred2)-ncol(temp.pred)+1,ncol(pred2)))
pred2_names<-c(pred2_names,paste(pred_names[i],1:l,sep="."))
pred3_names=c(pred3_names,paste(pred_names[i],1:l,sep="."))
}
}

colnames(pred1.0)=pred1.0_names
colnames(pred2)=pred2_names
if(!is.null(pred3))
{colnames(pred3)=pred3_names
colnames(pred.cont.der)=pred3_names}

## prepare mediators for y
m1 <- data.frame(m) #original format
m2 <- NULL             #all transformed
m2_names = NULL
m3 <- NULL             #transformed continuous mediators with mediator+delta(med)
m3_names = NULL
m.cont.der <- NULL     #derivative of the transformation function for cont mediators
if (is.null(mnames))
  mnames = paste("m",1:ncol(m1),sep='')
colnames(m1) = mnames
binm = NULL
catm = NULL
contm = NULL
binm1 = NULL
catm1 = NULL
contm1 = NULL
contm3 = NULL         #index for m3 and m.cont.dev
nm = ncol(m1)
n2=nrow(m1)
for (i in 1:nm)
  if (nlevels(droplevels(as.factor(m1[,i]))) == 2) { #binary mediator
    if (!is.na(mref[i]))
    {m2 <- cbind(m2,ifelse(m1[, i] == mref[i],0, 1))
    #m3 <- cbind(m3,as.factor(ifelse(m1[, i] == mref[i],0, 1)))
    m1[,i] <- ifelse(m1[, i] == mref[i],0, 1)}
    else {
      temp.m <- as.factor(m1[, i])
      m2 <- cbind(m2,ifelse(temp.m == levels(droplevels(temp.m))[1], 0, 1))
      #m3 <- cbind(m3,as.factor(ifelse(temp.m == levels(droplevels(temp.m))[1], 0, 1)))
      m1[,i] <- ifelse(temp.m == levels(droplevels(temp.m))[1], 0, 1)
    }
    binm1 = c(binm1, i)
    m2_names=c(m2_names,mnames[i])
    colnames(m2)=m2_names
    binm = c(binm,ncol(m2))
  }
else if (is.character(m1[, i]) | is.factor(m1[, i])) { #category mediator
  m1[, i] = droplevels(as.factor(m1[, i]))
  if (!is.na(mref[i]))
    m.temp1 <- bin_cat(m1, m2, i, mref[i])
  else
    m.temp1 <- bin_cat(m1, m2, i, levels(as.factor(m1[,i]))[1])
  m2 = m.temp1$M1
  m1 = m.temp1$M2
  #m3 = cbind(m3,m.temp1$M1[,m.temp1$cat[1]:m.temp1$cat[2]])
  catm1 = c(catm1,i)
  catm = rbind(catm,m.temp1$cat)
  m2_names = c(m2_names, colnames(m.temp1$M1)[m.temp1$cat[1]:m.temp1$cat[2]])
  colnames(m2)=m2_names
}
else #consider the transformation of continuous m
{contm1 = c(contm1, i)
if(!(i %in% fmy[[1]]))  #fmy has the transformation functions for m to y, [[1]] list all cont mediators in m to be transformed
{m2<-cbind(m2,m1[,i])
m.cont.der<-cbind(m.cont.der,rep(1,n1))
contm3 = rbind(contm3,c(ncol(m.cont.der),ncol(m.cont.der)))
m3 = cbind(m3,m1[,i]+deltam[i]) #deltam is the changing amount for mediators
m2_names=c(m2_names,mnames[i])
colnames(m2)=m2_names
m3_names=c(m3_names,mnames[i])
contm=rbind(contm,c(ncol(m2),ncol(m2)))
}
else if (i %in% fmy[[1]])
{p=match(i,fmy[[1]])
temp.m=x2fx(m1[,i],fmy[[1+p]])$values
m2<-cbind(m2,temp.m)
m.cont.der<-cbind(m.cont.der,x2fdx(as.matrix(m)[,i],fmy[[1+p]]))
if(is.null(m3))
  contm3 = rbind(contm3,c(1,ncol(m.cont.der)))
else
  contm3 = rbind(contm3,c(ncol(m3)+1,ncol(m.cont.der)))
m3<-cbind(m3,x2fx(m1[,i]+deltap[i],fmy[[1+p]])$values)
contm = rbind(contm,c(ncol(m2)-ncol(temp.m)+1,ncol(m2)))
l=length(fmy[[1+p]])
m2_names<-c(m2_names,paste(mnames[i],1:l,sep="."))
colnames(m2)=m2_names
m3_names=c(m3_names,paste(mnames[i],1:l,sep="."))
}
}

colnames(m2)=m2_names
if(!is.null(m3))
{colnames(m3)=m3_names
colnames(m.cont.der)=m3_names}

## prepare predictors for mediators
pm=rep(0,n1) #the first row are all 0s
fpm.2=fpm  #the transformed predictors in pm
binp=NULL
catp=NULL
contp=NULL
j=1
names.pm=("zero")
if(!is.null(binpred))
  for (i in binpred)
  {pm=cbind(pm,pred2[,i])
  j=j+1
  binp=c(binp,j)
  names.pm=c(names.pm,pred_names[i])}
if(!is.null(catpred))
  for (i in 1:nrow(catpred))
  {pm=cbind(pm,pred2[,catpred[i,1]:catpred[i,2]])
  catp=rbind(catp,c(j+1,j+1+catpred[i,2]-catpred[i,1]))
  j=j+1+catpred[i,2]-catpred[i,1]
  names.pm=c(names.pm,colnames(pred2)[catpred[i,1]:catpred[i,2]])
  }
if(!is.null(contpred1))
  for (i in contpred1)
  {pm=cbind(pm,pred1[,i])
  j=j+1
  contp=rbind(contp,rep(j,3))
  names.pm=c(names.pm,pred_names[i])}
pm.der=matrix(1,n1,ncol(pm))
pm.idx=as.list(rep(1,ncol(m1))) #the predictors in pm to predict the ith mediator
for (i in 1:ncol(m1))
  pm.idx[[i]]=2:ncol(pm)

if(!is.null(fpm))
{k=unique(fpm[[1]][,2]) #the first column is for the mediators,
#the second column the continuous predictors to be transformed
for (l in k){
  temp<-(2:length(fpm))[fpm[[1]][,2]==l]
  allfun=fpm[[temp[1]]]
  if (length(temp)>1)
    for(i in 2:length(temp))
      allfun<-c(allfun,fpm[[temp[i]]])
  unifun<-unique(allfun)
  unifun1<-unifun[unifun!="x"]
  unifun2<-c("x",unifun1)
  d_d<-x2fx(pred1[,l],unifun1)
  d.der<-x2fdx(pred1[,l],unifun1)
  d<-as.matrix(d_d$values)
  names.pm<-c(names.pm,paste(pred_names[l],1:ncol(d),sep="."))
  pm.der<-cbind(pm.der,d.der)
  place=match(l,contpred1)
  contp[place,2]=j+1
  contp[place,3]=j+ncol(d)
  pm<-cbind(pm,d)
  for(i in temp)
  {ttemp<-order_char(unifun1,fpm[[i]])
  fpm.2[[i]]=j+ttemp #what does this do?
  pm.indx[[fpm[[1]][i-1,1]]]=c(pm.indx[[fpm[[1]][i-1,1]]],j+ttemp)
  if(length(order_char('x',fpm[[i]]))==0)
    pm.indx[[fpm[[1]][i-1,1]]]= (pm.indx[[fpm[[1]][i-1,1]]])[pm.indx[[fpm[[1]][i-1,1]]]!=l]
  }
  j=j+ncol(d)
}}

colnames(pm)=names.pm
colnames(pm.der)=names.pm
pm.ind=matrix(1,length(pm.idx),max(sapply(pm.idx,length)))
for (i in 1:nrow(pm.ind))
  pm.ind[i,1:length(pm.idx[[i]])]=pm.idx[[i]]

p2=ifelse(is.null(binm1),0,length(binm1))
p3=ifelse(is.null(catm1),0,length(catm1))
p1=ifelse(is.null(contm1),0,length(contm1))

#prepare for covariates of mediators
if(is.null(mcov))
{mcov=data.frame(intercept=rep(1,nrow(m1)))
mind=matrix(T,p1+p2+p3,1)}
else
{mcov=data.frame(intercept=rep(1,nrow(m)),mcov)
mind=matrix(T,p1+p2+p3,ncol(mcov))
mcov_names=colnames(mcov)
if (!is.null(mclist))
{if (is.character(mclist[[1]]))
  mclist[[1]]=match(mclist[[1]],mnames)
mcov=data.frame(mcov,no=rep(0,nrow(mcov)))  #add a column of 0 in mcov
mind=matrix(rep(1:(ncol(mcov)-1),each=p1+p2+p3),p1+p2+p3,ncol(mcov)-1)
for (i in 1:length(mclist[[1]]))
{if(sum(is.na(mclist[[i+1]]))>=1)
  temp=1
else if(is.character(mclist[[i+1]]))
  temp=c(1,match(mclist[[i+1]],mcov_names))
else
  temp=c(1,mclist[[i+1]]+1)
mind[mclist[[1]][i],(1:(ncol(mcov)-1))[-temp]]=ncol(mcov)
}}
}  #use all covariates for all tv if mclist is NULL. Otherwise, the first item of mclist lists all tvs
#that are using different mcov, the following items gives the mcov for the tv in order. use NA is no mcov to be used. If not specified in mclist, use all mcov.

results = list(N=nrow(data.frame(y)), y_type=y_type, y=y, pred1=pred1, pred1.0=pred1.0,
               pred2=pred2, pred3=pred3, cova=cova,
               pred.cont.der=pred.cont.der, binpred2=binpred, catpred2=catpred,
               contpred2=contpred, binpred1=binpred1, catpred1=catpred1,contpred1=contpred1,
               contpred1.0=contpred1.0, binpred1.0=binpred1.0, catpred1.0=catpred1.0,
               contpred3=contpred3, npred=npred,
               m1=m1, m2=m2, m3=m3, m.cont.der=m.cont.der, binm2=binm, catm2=catm,
               contm2=contm,binm1=binm1, catm1=catm1, contm1=contm1, contm3=contm3,
               nm=nm, pm=pm, pm.der=pm.der, pm.idx=pm.idx, pm.ind=pm.ind, fpm.2=fpm.2,
               binp=binp, catp=catp, contp=contp,p1=p1,p2=p2,p3=p3,mcov=mcov,mind=mind)
return(results)
}

### build the bugs model if it is not defined
  jags_model <- function (contm=c('ncontm', 'ycontm'),binm=c('nbinm','ybinm'),
                          catm=c('ncatm','ycatm'),
                          cova=c('ncova','ycova'),
                          ytype=c('contc','binc','catc','survc','contnnc',
                                  'binnc','catnc','survnc')) {
    contm <- match.arg(contm)
    binm <- match.arg(binm)
    catm <- match.arg(catm)
    cova <- match.arg(cova)

    # raw script
    script <-
      "model {

  for(i in 1:N){
    $yfunc

    $contm

    $binm

    $catm
}

$ypriors

$cntmprior

$bmprior

$cmprior

  var4 ~ dgamma(1,0.1)
  prec4 <-1/var4
  }"

  # define macros
  macros <- list(list("$yfunc",
                      switch(ytype,
                             contc='mu_y[i] <- beta0 + inprod(c, x[i,]) + inprod(beta,M1[i,]) + inprod(eta,cova[i,])
    y[i] ~ dnorm(mu_y[i],prec4)',
                             contnc='mu_y[i] <- beta0 + inprod(c, x[i,]) + inprod(beta,M1[i,])
                             y[i] ~ dnorm(mu_y[i],prec4)',
                             binc='logit(mu_y[i]) <- beta0 + inprod(c, x[i,]) + inprod(beta,M1[i,]) + inprod(eta,cova[i,])
    y[i] ~ dbern(mu_y[i])',
                             binnc='logit(mu_y[i]) <- beta0 + inprod(c, x[i,]) + inprod(beta,M1[i,])
    y[i] ~ dbern(mu_y[i])',
                             catc='mu_y1[i,1] <- 1
for (k in 2:caty)
{mu_y1[i,k] <- exp(beta0[k-1] + inprod(c[k-1,], x[i,]) + inprod(beta[k-1,],M1[i,]) + inprod(eta[k-1,],cova[i,]))}
sum_y[i] <- sum(mu_y1[i,1:caty])
for (l in 1:caty)
{mu_y[i,l] <- mu_y1[i,l]/sum_y[i]}
y[i] ~ dcat(mu_y[i,])',
                             catnc='mu_y1[i,1] <- 1
for (k in 2:caty)
{mu_y1[i,k] <- exp(beta0[k-1] + inprod(c[k-1,], x[i,]) + inprod(beta[k-1,],M1[i,]))}
sum_y[i] <- sum(mu_y1[i,1:caty])
for (l in 1:caty)
{mu_y[i,l] <- mu_y1[i,l]/sum_y[i]}
y[i]~dcat(mu_y[i,])',
                             survc=   'elinpred[i] <- exp(inprod(c, x[i,]) + inprod(beta,M1[i,]) + inprod(eta,cova[i,]))
                             base[i] <- lambda*r*pow(y[i,1], r-1)
                             loghaz[i] <- log(base[i]*elinpred[i])
                             phi[i] <- 100000-y[i,2]*loghaz[i]-log(exp(-lambda*pow(y[i,1],r)*elinpred[i])-exp(-lambda*pow(tmax,r)*elinpred[i])) +log(1-exp(-lambda*pow(tmax,r)*elinpred[i]))
                             zero[i] ~ dpois(phi[i])'                             ,
                             survnc= 'elinpred[i] <- exp(inprod(c, x[i,]) + inprod(beta,M1[i,]))
                             base[i] <- lambda*r*pow(y[i,1], r-1)
                             loghaz[i] <- log(base[i]*elinpred[i])
                             phi[i] <- 100000-y[i,2]*loghaz[i]-log(exp(-lambda*pow(y[i,1],r)*elinpred[i])-exp(-lambda*pow(tmax,r)*elinpred[i])) +log(1-exp(-lambda*pow(tmax,r)*elinpred[i]))
                             zero[i] ~ dpois(phi[i])')),
                 list("$contm",
                      switch(contm,
                             ycontm='for (j in 1:p1){
      mu_M1[i,contm[j]] <- inprod(alpha0.a[j,mind[contm[j],]],mcov[i,mind[contm[j],]])+inprod(alpha1.a[j,1:c1],x1[i,])
      M2[i,contm[j]] ~ dnorm(mu_M1[i,contm[j]],prec1[j])
      for (k in contm1[j,1]:contm1[j,2]){
        mu_M1_c[i,k] <- inprod(alpha0[k,mind[contm[j],]],mcov[i,mind[contm[j],]])+inprod(alpha1[k,1:c1],x1[i,])
        M1[i,k] ~ dnorm(mu_M1_c[i,k],prec2[k])
        mu_M1_c1[i,k] <- inprod(alpha0[k,mind[contm[j],]],mcov[i,mind[contm[j],]])+inprod(alpha1[k,1:c1],x1[i,]+deltax)
      }
    }
',
                             ncontm='')),
                 list("$cntmprior",
                      switch(contm,
                             ycontm='  for(j in 1:P)
    {alpha1[j,1:c1] ~ dmnorm(mu1.1[j,1:c1], Omega1.1[1:c1, 1:c1])
     alpha0[j,1:nmc] ~ dmnorm(mu0.1[j,1:nmc], Omega0.1[1:nmc, 1:nmc])
    }
     for (i in 1:P){
       var2[i] ~ dgamma(0.1,0.1)
       prec2[i] <- 1/var2[i]
     }

  for(j in 1:p1)
    {alpha1.a[j,1:c1] ~ dmnorm(mu1.a[j,1:c1], Omega1.a[1:c1, 1:c1])
     alpha0.a[j,1:nmc] ~ dmnorm(mu0.a[j,1:nmc], Omega0.a[1:nmc, 1:nmc])
     for(k in contm1[j,1]:contm1[j,2])
      {temp1[k]=mean(mu_M1_c1[1:N,k]-mu_M1_c[1:N,k])/deltax}
     var1[j] ~ dgamma(1,0.1)
     prec1[j] <- 1/var1[j]
  }',
                             ncontm='')),
                 list("$binm",
                      switch(binm,
                             ybinm="  for (k in 1:p2){
    logit(mu_M1[i,binm[k]]) <- inprod(alpha0.b[k,mind[binm[k],]],mcov[i,mind[binm[k],]])+inprod(alpha1.b[k,1:c1],x1[i,])
    M2[i,binm[k]] ~ dbern(mu_M1[i,binm[k]])
    logit(mu_M1_b1[i,binm[k]]) <- inprod(alpha0.b[k,mind[binm[k],]],mcov[i,mind[binm[k],]])+inprod(alpha1.b[k,1:c1],x1[i,]+deltax)
  }",
                             nbinm='')),
                 list("$bmprior",
                      switch(binm,
                             nbinm='',
                             ybinm='  for(j in 1:p2)
                             {alpha1.b[j,1:c1] ~ dmnorm(mu1.b[j,1:c1], Omega1.b[1:c1, 1:c1])
                              alpha0.b[j,1:nmc] ~ dmnorm(mu0.b[j,1:nmc], Omega0.b[1:nmc, 1:nmc])
                              temp1[binm1[j]]<-mean(mu_M1_b1[1:N,binm[j]]-mu_M1[1:N,binm[j]])/deltax
                             }')),
                 list("$catm",
                      switch(catm,
                             ycatm=" for (j in 1:p3){
      mu_Mc[i,j,1] <- 1 #baseline is the 1st category
      mu_Mc1[i,j,1] <- 1 #baseline is the 1st category
      for (k in 2:cat2[j]){
        mu_Mc[i,j,k] <- exp(inprod(alpha0.c[j,k-1,mind[catm[j],]],mcov[i,mind[catm[j],]])+inprod(alpha1.c[j,k-1,1:c1],x1[i,]))
        mu_Mc1[i,j,k] <- exp(inprod(alpha0.c[j,k-1,mind[catm[j],]],mcov[i,mind[catm[j],]])+inprod(alpha1.c[j,k-1,1:c1],x1[i,]+deltax))
      }
      sum_Mc[i,j]  <- sum(mu_Mc[i,j,1:cat2[j]])
      sum_Mc1[i,j] <- sum(mu_Mc1[i,j,1:cat2[j]])
      for (l in 1:cat2[j])
      {mu_Mc0[i,j,l] <- mu_Mc[i,j,l]/sum_Mc[i,j]
       mu_Mc01[i,j,l] <- mu_Mc1[i,j,l]/sum_Mc1[i,j]}
       M2[i,catm[j]] ~ dcat(mu_Mc0[i,j,1:cat2[j]])
      }",
                             ncatm='')),
                 list("$cmprior",
                      switch(catm,
                             ncatm='',
                             ycatm='    for (i in 1:p3){
    for(j in 1:cat1)
      {alpha1.c[i,j,1:c1] ~ dmnorm(mu1.c[j,1:c1], Omega1.c[1:c1, 1:c1])
       alpha0.c[i,j,1:nmc] ~ dmnorm(mu0.c[j,1:nmc], Omega0.c[1:nmc, 1:nmc])}
    for (l in 1:(cat2[i]-1))
      {temp1[catm1[i,1]+l-1]<-mean(mu_Mc01[1:N,i,l]-mu_Mc0[1:N,i,l])/deltax}
  }')),
                 list("$ypriors",
                      switch(ytype,
                             contc='  for (i in 1:P)
  {beta[i] ~ dunif(-u[1+c2+i]*sqrt(var4),u[1+c2+i]*sqrt(var4))
   u[1+c2+i] ~ dgamma(2,(lambda1/abs(para1[i]+0.001)^r1))}
  beta0 ~ dunif(-u[1]*sqrt(var4),u[1]*sqrt(var4))
  u[1] ~ dgamma(2,lambda1/para0^r1)
  for(j in 1:c2)
   {c[j] ~ dunif(-u[1+j]*sqrt(var4),u[1+j]*sqrt(var4))
    u[1+j] ~ dgamma(2,lambda1/(para2[j])^r1)}
  eta[1:cv1] ~ dmnorm(mucv[1:cv1], Omegacv[1:cv1,1:cv1])
  lambda1 ~ dgamma(0.01,0.01)',
                             contnc='  for (i in 1:P)
  {beta[i] ~ dunif(-u[1+c2+i]*sqrt(var4),u[1+c2+i]*sqrt(var4))
   u[1+c2+i] ~ dgamma(2,lambda1/abs(para1[i]+0.001)^r1)}
  beta0 ~ dunif(-u[1]*sqrt(var4),u[1]*sqrt(var4))
  u[1] ~ dgamma(2,lambda1/para0^r1)
  for(j in 1:c2)
   {c[j] ~ dunif(-u[1+j]*sqrt(var4),u[1+j]*sqrt(var4))
    u[1+j] ~ dgamma(2,lambda1/(para2[j])^r1)}
  lambda1 ~ dgamma(0.01,0.01)',
                             binc='  for (i in 1:P)
  {beta[i] ~ dunif(-u[1+c2+i]*sqrt(var4),u[1+c2+i]*sqrt(var4))
   u[1+c2+i] ~ dgamma(2,lambda1/abs(para1[i]+0.001)^r1)}
  beta0 ~ dunif(-u[1]*sqrt(var4),u[1]*sqrt(var4))
  u[1] ~ dgamma(2,lambda1/para0^r1)
  for(j in 1:c2)
   {c[j] ~ dunif(-u[1+j]*sqrt(var4),u[1+j]*sqrt(var4))
    u[1+j] ~ dgamma(2,lambda1/(para2[j])^r1)}
  eta[1:cv1] ~ dmnorm(mucv[1:cv1], Omegacv[1:cv1,1:cv1])
  lambda1 ~ dgamma(0.01,0.01)',
                             binnc='  for (i in 1:P)
  {beta[i] ~ dunif(-u[1+c2+i]*sqrt(var4),u[1+c2+i]*sqrt(var4))
   u[1+c2+i] ~ dgamma(2,lambda1/abs(para1[i]+0.001)^r1)}
  beta0 ~ dunif(-u[1]*sqrt(var4),u[1]*sqrt(var4))
  u[1] ~ dgamma(2,lambda1/para0^r1)
  for(j in 1:c2)
   {c[j] ~ dunif(-u[1+j]*sqrt(var4),u[1+j]*sqrt(var4))
    u[1+j] ~ dgamma(2,lambda1/(para2[j])^r1)}
  lambda1 ~ dgamma(0.01,0.01)',
                             catc='  for(j in 1:(caty-1))
  {for (i in 1:P)
   {beta[j,i] ~ dunif(-u[j,1+c2+i]*sqrt(var4),u[j,1+c2+i]*sqrt(var4))
    u[j,1+c2+i] ~ dgamma(2,lambda1/abs(para1[j,i]+0.001)^r1)}
  beta0[j] ~ dunif(-u[j,1]*sqrt(var4),u[j,1]*sqrt(var4))
  u[j,1] ~ dgamma(2,lambda1/para0[j]^r1)
  for(l in 1:c2)
   {c[j,l] ~ dunif(-u[j,1+l]*sqrt(var4),u[j,1+l]*sqrt(var4))
    u[j,1+l] ~ dgamma(2,lambda1/(para2[j,l])^r1)}
  eta[j,1:cv1] ~ dmnorm(mucv[1:cv1], Omegacv[1:cv1,1:cv1])}
  lambda1 ~ dgamma(0.01,0.01)',
                             catnc='  for(j in 1:(caty-1))
  {for (i in 1:P)
   {beta[j,i] ~ dunif(-u[j,1+c2+i]*sqrt(var4),u[j,1+c2+i]*sqrt(var4))
    u[j,1+c2+i] ~ dgamma(2,lambda1/abs(para1[j,i]+0.001)^r1)}
  beta0[j] ~ dunif(-u[j,1]*sqrt(var4),u[j,1]*sqrt(var4))
  u[j,1] ~ dgamma(2,lambda1/para0[j]^r1)
  for(l in 1:c2)
   {c[j,l] ~ dunif(-u[j,1+l]*sqrt(var4),u[j,1+l]*sqrt(var4))
    u[j,1+l] ~ dgamma(2,lambda1/(para2[j,l])^r1)}
  lambda1 ~ dgamma(0.01,0.01)',
                             survc='  for (i in 1:P)
  {beta[i] ~ dunif(-u[1+c2+i]*sqrt(var4),u[1+c2+i]*sqrt(var4))
   u[1+c2+i] ~ dgamma(2,lambda1/abs(para1[i]+0.001)^r1)}
  beta0 ~ dunif(-u[1]*sqrt(var4),u[1]*sqrt(var4))
  u[1] ~ dgamma(2,lambda1/para0^r1)
  for(j in 1:c2)
   {c[j] ~ dunif(-u[1+j]*sqrt(var4),u[1+j]*sqrt(var4))
    u[1+j] ~ dgamma(2,lambda1/(para2[j])^r1)}
  r~dunif(0,10)  # dunif(0.5,1.5)
  lambda~dgamma(1,0.01)
  eta[1:cv1] ~ dmnorm(mucv[1:cv1], Omegacv[1:cv1,1:cv1])
  lambda1 ~ dgamma(0.01,0.01)',
                             survnc= '  for (i in 1:P)
  {beta[i] ~ dunif(-u[1+c2+i]*sqrt(var4),u[1+c2+i]*sqrt(var4))
   u[1+c2+i] ~ dgamma(2,lambda1/abs(para1[i]+0.001)^r1)}
  beta0 ~ dunif(-u[1]*sqrt(var4),u[1]*sqrt(var4))
  u[1] ~ dgamma(2,lambda1/para0^r1)
  for(j in 1:c2)
   {c[j] ~ dunif(-u[1+j]*sqrt(var4),u[1+j]*sqrt(var4))
    u[1+j] ~ dgamma(2,lambda1/(para2[j])^r1)}
  r~dunif(0,10)  # dunif(0.5,1.5)
  lambda~dgamma(1,0.01)
  lambda1 ~ dgamma(0.01,0.01)')) # dunif(1.0E-8,3.0E-7)
  )
  # apply macros
  for (m in seq(macros)) {
    script <- gsub(macros[[m]][1], macros[[m]][2], script, fixed=TRUE)
  }
  script
  }

  ### Incomplete Gamma function
  Igamma<-function(z,u)pgamma(u,z)*gamma(z)
  ### Truncated Weibull distribution moments
  weib.trunc.mean<-function(vec, right)vec[-1]*Igamma(1/vec[1]+1,(right/vec[-1])^vec[1])/(1-exp(-(right/vec[-1])^vec[1]))

  matrix.prod<-function(m1)
  {return(m1[1,]%*%t(m1[-1,]))}

  expp<-function(v1,v2)
  {v1^v2}


  data0<- data_org(pred=pred, m=m, y=y, refy=refy, predref=predref, fpy=fpy,
                   deltap=deltap, fmy=fmy, deltam=deltam, fpm=fpm, mref=mref,
                   cova=cova,mcov = mcov, mclist=mclist)
  y.type=data0$y_type #1 for continuous outcome, 4 for time-to-event, 2 for binary, 3 for categorical
  N=data0$N
  x=data0$pred2
  c2=ncol(x)
  x1=data0$pred1.0
  c1=ncol(x1)    #c1 is the number of predictors, k-class categorical variable are considered as k-1 predictors
  y=data0$y
  M1=data0$m2
  M2=data0$m1
  M3=data0$m3
  #alpha=rep(1,ncol(M1))
  contm=data0$contm1
  contm1=data0$contm2
  contm3=data0$contm3
  p1=data0$p1
  binm=data0$binm1
  binm1=data0$binm2
  p2=data0$p2
  p3=data0$p3
  if(p3>0){
  cat1=max(data0$catm2[,2]-data0$catm2[,1]+1)
  cat2=data0$catm2[,2]-data0$catm2[,1]+2
  catm=data0$catm1
  catm1=data0$catm2}
  else{
    cat1=NULL
    cat2=NULL
    catm=NULL
    catm1=NULL
  }
  P=ncol(data0$m2)
  cova=data0$cova
  mcov=data0$mcov
  mind=data0$mind
  nmc=ncol(mcov)
 # pm=data0$pm
 # pm.ind=data1$pm.ind

#  if(is.null(mu))
#    mu=rep(0,P)
#  if(is.null(Omega))
#    Omega=diag(preci,P)

#with linear regression
  if(y.type==3)
  {caty=nlevels(y)}

  if(!partial){
  if(is.null(cova))
    temp=cbind(M1,x)
  else
    temp=cbind(M1,x,cova)
if(y.type==1)
  para=abs((summary(glm(y~.,data=data.frame(temp))))$coefficient[,1])
else if (y.type==2)
  para=abs((summary(glm(y~.,data=data.frame(temp),family=binomial(link = "logit"))))$coefficient[,1])
else if(y.type==4)
  para=abs((summary(coxph(Surv(y)~.,data=data.frame(temp))))$coefficient[,1])
else if(y.type==3)
  {para=NULL
  for (i in 2:caty)
  {temp.t=data.frame(temp[y==levels(y)[c(1,i)],])
   y.t=y[y==levels(y)[c(1,i)]]
   para=rbind(para,abs((summary(glm(y.t~.,data=temp.t,family=binomial(link = "logit"))))$coefficient[,1]))}
  }
}

#with ridge regression
#if(y.type==1)
#  fit1=cv.glmnet(x=cbind(M1,x,cova),y=y,alpha=0)
#  else if (y.type==2)
#    fit1=cv.glmnet(x=cbind(M1,x,cova),y=y,alpha=0,family="binomial")
#  else if (y.type==4)
#    fit1=cv.glmnet(x=cbind(M1,x,cova),y=y,alpha=0,family="cox")
#para=abs(coef(fit1,s="lambda.min")[,1])
  if(partial){
    if(y.type!=3){
    para0=1
    para1=rep(1,ncol(M1))
    para2=rep(1,ncol(as.matrix(x)))}
    else{
      para0=rep(1,caty-1)
      para1=matrix(1,caty-1,ncol(M1))
      para2=matrix(1,caty-1,ncol(as.matrix(x)))
    }
  }
  else{
    if(y.type<3){
    para0=para[1]
    para1=para[2:(ncol(M1)+1)]
    para2=para[(ncol(M1)+2):(ncol(M1)+1+ncol(as.matrix(x)))]}
    else if(y.type==3){
      para0=para[,1]
      para1=para[,2:(ncol(M1)+1)]
      para2=para[,(ncol(M1)+2):(ncol(M1)+1+ncol(as.matrix(x)))]}
    else if(y.type==4){
      para0=1
      para1=para[1:(ncol(M1))]
      para2=para[(ncol(M1)+1):(ncol(M1)+ncol(as.matrix(x)))]}
  }
para.a=rep(0,ncol(M1))

if(p1>0){
  if(is.null(mu0.1))
    mu0.1=matrix(0,P,nmc)
  if(is.null(Omega0.1))
    Omega0.1=diag(preci,nmc)
  if(is.null(mu1.1))
    mu1.1=matrix(0,P,c1)
  if(is.null(Omega1.1))
    Omega1.1=diag(preci,c1)
  if(is.null(mu0.a))
    mu0.a=matrix(0,p1,nmc)
  if(is.null(Omega0.a))
    Omega0.a=diag(preci,nmc)
  if(is.null(mu1.a))
    mu1.a=matrix(0,p1,c1)
  if(is.null(Omega1.a))
    Omega1.a=diag(preci,c1)
  for (j in 1:p1)
    for (k in contm1[j,1]:contm1[j,2])
     {temp.alpha=lm(M1[,k]~.-1,data=data.frame(x1,mcov[,mind[contm[j],]]))
      temp.1=predict(temp.alpha)
      temp.2=predict(temp.alpha,newdata=data.frame(x1+deltax,mcov[,mind[contm[j],]]))
      para.a[k]=mean(temp.2-temp.1)/deltax
     }
}

if(p2>0){
  if(is.null(mu0.b))
    mu0.b=matrix(0,p2,nmc)
  if(is.null(Omega0.b))
    Omega0.b=diag(preci,nmc)
  if(is.null(mu1.b))
    mu1.b=matrix(0,p2,c1)
  if(is.null(Omega1.b))
    Omega1.b=diag(preci,c1)
  for (k in 1:p2)
  {temp.alpha=glm(M1[,binm1[k]]~.-1,data=data.frame(x1,mcov[,mind[binm[k],]]),family=binomial(link = "logit"))
   temp.1=predict(temp.alpha,type="response")
   temp.2=predict(temp.alpha,newdata=data.frame(x1+deltax,mcov[,mind[binm[k],]]),type="response")
   para.a[binm1[k]]=mean(temp.2-temp.1)/deltax
  }
}

if(p3>0){
   if(is.null(mu0.c))
    mu0.c=matrix(0,cat1,nmc)
  if(is.null(Omega0.c))
    Omega0.c=diag(preci,nmc)
  if(is.null(mu1.c))
    mu1.c=matrix(0,cat1,c1)
  if(is.null(Omega1.c))
    Omega1.c=diag(preci,c1)
  for (j in 1:p3)
  { mu_mc.1=matrix(1,ncol=cat2[j],nrow=nrow(M1))
    mu_mc.2=mu_mc.1
    for (k in 2:cat2[j])
    {temp.alpha=glm(M1[,catm1[j,1]+k-2]~.-1,data=data.frame(x1,mcov[,mind[catm[j],]]),family=binomial(link = "logit"))
    mu_mc.1[,k]=predict(temp.alpha,type="response")
    mu_mc.2[,k]=predict(temp.alpha,newdata=data.frame(x1+deltax,mcov[,mind[catm[j],]]),type="response")
    }
    temp.1=(diag(1/apply(mu_mc.1,1,sum)))%*%mu_mc.1
    temp.2=(diag(1/apply(mu_mc.2,1,sum)))%*%mu_mc.2
    para.a[catm1[j,1]:catm1[j,2]]=apply((temp.2-temp.1)[,-1],2,mean)
  }
}

para1=abs(para1*para.a)



if(y.type==4)
 {zero=rep(0,nrow(y))#is.cen=ifelse(y[,2]==1,0,1)   #censored or not
  if(is.null(tmax))
   tmax=max(y[,1], na.rm=T)+100
  if(is.null(multi))
    multi=TRUE}             #maximum time to event
 #Cen=ifelse(y[,2]==1,tmax,y[,1])  #censoring time
 #y=ifelse(y[,2]==1,y[,1],NA)}

  data0.1<- list (N=N,x=x,x1=x1,y=y,M1=M1,M2=M2,P=P,c1=c1,c2=c2,para0=para0,r1=r1,
                  nmc=nmc,mcov=mcov,mind=mind,deltax=deltax,para1=para1,para2=para2)
  para=c("beta0","c","beta","lambda1","temp1") #"alpha",

if(y.type==3)
  {data0.1[['caty']]=caty}

if(y.type==4)
  {data0.1[['zero']]=zero
  #data0.1[['is.cen']]=is.cen
  data0.1[['tmax']]=tmax
  para=c(para,"lambda","r")}

  if(p1>0)
  {data0.1[['contm']]=contm
   data0.1[['contm1']]=contm1
   data0.1[['p1']]=p1
   data0.1[['mu0.1']]=mu0.1
   data0.1[['mu1.1']]=mu1.1
   data0.1[['Omega0.1']]=Omega0.1
   data0.1[['Omega1.1']]=Omega1.1
   data0.1[['mu0.a']]=mu0.a
   data0.1[['mu1.a']]=mu1.a
   data0.1[['Omega0.a']]=Omega0.a
   data0.1[['Omega1.a']]=Omega1.a
   para=c(para,"alpha0","alpha1","alpha0.a","alpha1.a")
  }

  if(p2>0)
    {data0.1[['binm']]=binm
     data0.1[['binm1']]=binm1
     data0.1[['p2']]=p2
     data0.1[['mu0.b']]=mu0.b
     data0.1[['mu1.b']]=mu1.b
     data0.1[['Omega0.b']]=Omega0.b
     data0.1[['Omega1.b']]=Omega1.b
     para=c(para,"alpha0.b","alpha1.b")}

  if(p3>0)
    {data0.1[['cat1']]=cat1
     data0.1[['cat2']]=cat2
     data0.1[['catm']]=catm
     data0.1[['catm1']]=catm1
     data0.1[['p3']]=p3
     data0.1[['mu0.c']]=mu0.c
     data0.1[['mu1.c']]=mu1.c
     data0.1[['Omega0.c']]=Omega0.c
     data0.1[['Omega1.c']]=Omega1.c
     para=c(para,"alpha0.c","alpha1.c")}

  if(!is.null(cova)){
    cv1=ncol(cova)
    if(is.null(mucv))
      mucv=rep(0,cv1)
    if(is.null(Omegacv))
      Omegacv=diag(preci,cv1)
    data0.1[['cv1']]=cv1
    data0.1[['mucv']]=mucv
    data0.1[['Omegacv']]=Omegacv
    data0.1[['cova']]=cova
    para=c(para,"eta")
  }

  if(is.null(inits))
    inits<- function(){list()}

  if (y.type==1 & is.null(cova))
      ytype="contnc"
  else if(y.type==1 & !is.null(cova))
    ytype="contc"
  else if(y.type==2 & is.null(cova))
    ytype="binnc"
  else if(y.type==2 & !is.null(cova))
    ytype="binc"
  else if(y.type==3 & is.null(cova))
    ytype="catnc"
  else if(y.type==3 & !is.null(cova))
    ytype="catc"
  else if(y.type==4 & is.null(cova))
    ytype="survnc"
  else if(y.type==4 & !is.null(cova))
    ytype="survc"

  if (is.null(filename))
  {filename=paste(tempdir(),"model.txt",sep="/")
   writeLines(jags_model(contm=ifelse(p1==0,'ncontm', 'ycontm'),
                        binm=ifelse(p2==0, 'nbinm','ybinm'),
                        catm=ifelse(p3==0, 'ncatm','ycatm'),
                        cova=ifelse(is.null(cova),'ncova','ycova'),
                        ytype=ytype), filename)
  }

  med0<- jags(data0.1, inits,
              model.file = filename, #"O:/My Documents/My Research/Research/Bayesian Mediation Analysis/codes/promis/bx_cy.txt",
              parameters.to.save = para,
              n.chains = n.chains, n.iter = n.iter, n.burnin=n.burnin, n.thin = n.thin)

  #check the results
  #calculate the mediation effects
  N1=(n.iter-n.burnin)/n.thin  #10000
  #m.mcov=apply(mcov,2,mean,na.rm=T)  #effects are calculated at the mean of mcov
  #med0$BUGSoutput$sims.list$r=1/med0$BUGSoutput$sims.list$r
  #attach(med0$BUGSoutput$sims.list)
#browser()
  if(y.type==3){
    aie1=array(0,c(N1,p1+p2+p3,c1,caty-1))
    ie1=array(0,dim=c(N1,N,p1+p2+p3,c1,caty-1))
    tt2=1
    de1=array(0,c(N1,c1,caty-1))
    n.cont=1

    aie2=array(0,c(N1,p1+p2+p3,c1,caty-1))
    ie2=array(0,dim=c(N1,N,p1+p2+p3,c1,caty-1))
    de2=array(0,c(N1,c1,caty-1))

    te3=array(0,c(N1,N,caty-1))
    ate3=array(0,c(N1,c1,caty-1))
    omu3=ate3
    ade3=ate3
    aie3=array(0,c(N1,p1+p2+p3,c1,caty-1))
    ie3=array(0,dim=c(N1,N,p1+p2+p3,c1,caty-1))
    de3=array(0,c(N1,N,caty-1))
    mu_M2=array(0,c(dim(M1),N1))
    mu_M3=mu_M2

    te4=array(0,c(N1,N,caty-1))
    de4=te4
    ade4=array(0,c(N1,c1,caty-1))
    ie4=array(0,dim=c(N1,N,p1+p2+p3,c1,caty-1))
    mu.M0=array(0,c(dim(M1),N1))
    mu.M1=mu.M0
    ate4=ade4
    omu4=ate4
    aie4=array(0,c(N1,p1+p2+p3,c1,caty-1))

    for (l in 1:c1){
      x1.temp=x1
      x3.temp=x1
      if(l%in%data0$contpred1.0){
        x3.temp[,l]=x1[,l]+deltap[l]
      }
      else if(l%in%data0$binpred1.0)
      {x1.temp[,l]=0
      x3.temp[,l]=1
      deltap[l]=1}
      else{ #categorical predictor
        for (i in 1:nrow(data0$catpred1.0))
          if(l%in%(data0$catpred1.0[i,1]:data0$catpred1.0[i,2]))
          {x1.temp[,data0$catpred1.0[i,1]:data0$catpred1.0[i,2]]=0
          x3.temp[,data0$catpred1.0[i,1]:data0$catpred1.0[i,2]]=0
          x3.temp[,l]=1}
        deltap[l]=1
      }

      #method 1: the same for binary or continuous predictors
      if(p1>0){
        for(k in 1:(caty-1))
          if (p1+p2+p3==1 & c1==1 & contm1[1,1]==contm1[1,2])
            aie1[,contm[1],l,k]=med0$BUGSoutput$sims.list$alpha1[,1]*med0$BUGSoutput$sims.list$beta[,k,contm1[1,1]]#*mean(data0$m.cont.der[,contm3[1,1]]) since alpha1 is the coefficient of x on the change of f(M), beta is the coefficient of f(M) on the change of y
        else
          for (j in 1: p1)
            aie1[,contm[j],l,k]=apply(as.matrix(med0$BUGSoutput$sims.list$alpha1[,contm1[j,1]:contm1[j,2],l])*as.matrix(med0$BUGSoutput$sims.list$beta[,k,contm1[j,1]:contm1[j,2]]),1,sum)#*mean(data0$m.cont.der[,contm3[j,1]])
      }


      if(p2>0){
        if(p2==1 & c1==1) #since c1=1, the predictor can only be binary or continuous
        {if (nmc==1)
        {temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,1]+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x1.temp)
        temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,1]+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x3.temp)}
          else
          {temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,1,mind[binm[1],]]%*%t(mcov[,mind[binm[1],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x1.temp)
          temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,1,mind[binm[1],]]%*%t(mcov[,mind[binm[1],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x3.temp)}
          for(j in 1:(caty-1)){
            if(!1%in%data0$contpred1.0) #if the predictor is binary
              ie1[,,binm[1],1,j]=med0$BUGSoutput$sims.list$beta[,j,binm1[1]]*(exp(temp.x3)/(1+exp(temp.x3))-exp(temp.x1)/(1+exp(temp.x1)))
            else #if the predictor is continuous
              ie1[,,binm[1],1,j]=med0$BUGSoutput$sims.list$alpha1.b[,1]*med0$BUGSoutput$sims.list$beta[,j,binm1[1]]*exp(temp.x1)/(1+exp(temp.x1))^2
            aie1[,binm[1],l,j] <-apply(ie1[,,binm[1],1,j],1,mean)}
        }
        else
          for (k in 1:p2){
            if (nmc==1 & p2==1)
            {temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,k]+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x1.temp)
            temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,k]+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x3.temp)
            }
            else
            {temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,k,mind[binm[k],]]%*%t(mcov[,mind[binm[k],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x1.temp)
            temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,k,mind[binm[k],]]%*%t(mcov[,mind[binm[k],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x3.temp)
            }
            for(j in 1:(caty-1)){
              if(!l%in%data0$contpred1.0) #if the predictor is binary or categorical
                ie1[,,binm[k],l,j]=med0$BUGSoutput$sims.list$beta[,j,binm1[k]]*(exp(temp.x3)/(1+exp(temp.x3))-exp(temp.x1)/(1+exp(temp.x1)))
              else if(data0$contpred1.0) #for continuous predictor
                ie1[,,binm[1],l,j]=med0$BUGSoutput$sims.list$alpha1.b[,k,l]*med0$BUGSoutput$sims.list$beta[,j,binm1[k]]*exp(temp.x1)/(1+exp(temp.x1))^2
              aie1[,binm[k],l,j] <- apply(ie1[,,binm[k],l,j],1,mean)}
          }
      }

      if(p3>0){
        for (j in 1:p3){
          mu_Mc1<-array(0,c(N1,N,cat2[j]-1))
          mu_Mc0<-array(0,c(N1,N,cat2[j]-1))
          for (k in 2:cat2[j]){
            mu_Mc0[,,k-1] <- exp(as.matrix(med0$BUGSoutput$sims.list$alpha0.c[,j,k-1,mind[catm[j],]])%*%t(mcov[,mind[catm[j],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.c[,j,k-1,])%*%t(x1.temp))
            mu_Mc1[,,k-1] <- exp(as.matrix(med0$BUGSoutput$sims.list$alpha0.c[,j,k-1,mind[catm[j],]])%*%t(mcov[,mind[catm[j],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.c[,j,k-1,])%*%t(x3.temp))
          }
          sum_Mc1 <-apply(mu_Mc1,c(1,2),sum)+1
          sum_Mc0 <-apply(mu_Mc0,c(1,2),sum)+1
          if(l%in%data0$contpred1.0) #for continuous predictor
          {tt.0=rep(0,N1)
          for (k in 2:cat2[j])
            tt.0=mu_Mc0[,,k-1]/sum_Mc0*med0$BUGSoutput$sims.list$alpha1.c[,j,k-1,l]
          for (q1 in 1:(caty-1))
            for (k in 2:cat2[j])
              ie1[,,catm[j],l,q1]=ie1[,,catm[j],l,q1]+(mu_Mc0[,,k-1]/sum_Mc0)*med0$BUGSoutput$sims.list$beta[,q1,catm1[j,1]+k-2]*(med0$BUGSoutput$sims.list$alpha1.c[,j,k-1,l]-tt.0)
          }
          else #for binary or categorical predictor
            for (q1 in 1:(caty-1))
              for (k in 2:cat2[j])
                ie1[,,catm[j],l,q1]=ie1[,,catm[j],l]+(mu_Mc1[,,k-1]/sum_Mc1-mu_Mc0[,,k-1]/sum_Mc0)*med0$BUGSoutput$sims.list$beta[,q1,catm1[j,1]+k-2]
          for (q1 in 1:(caty-1))
            aie1[,catm[j],l,q1]<-apply(ie1[,,catm[j],l,q1],1,mean)
        }}

      if(l%in%data0$contpred1)
      {tt1=match(l,data0$contpred1)
      for (q1 in 1:(caty-1))
        de1[,l,q1]=as.matrix(med0$BUGSoutput$sims.list$c[,q1,data0$contpred3[tt1,1]:data0$contpred3[tt1,2]])%*%
          apply(as.matrix(data0$pred.cont.der[,data0$contpred3[tt1,1]:data0$contpred3[tt1,2]]),2,mean)
      tt2=tt2+data0$contpred3[tt1,2]-data0$contpred3[tt1,1]+1
      }
      else
      {for (q1 in 1:(caty-1))
        de1[,l,q1]=med0$BUGSoutput$sims.list$c[,q1,tt2]
      tt2=tt2+1}

      #method2: the same for binary and continuous predictors
      if(p1>0){
        if(p1==1 & c1==1 & contm1[1,1]==contm1[1,2])
        {temp.M3=M1
        temp.M3[,contm1[1,1]]=M3[,contm3[1,1]]
        temp.mu1<-array(0,c(N1,N,caty))
        temp.mu3<-temp.mu1
        for(q1 in 2:caty)
        {temp.mu1[,,q1]=med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x)+med0$BUGSoutput$sims.list$beta[,q1-1,]%*%t(M1)+matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,nrow(x))
        temp.mu3[,,q1]=med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x)+med0$BUGSoutput$sims.list$beta[,q1-1,]%*%t(temp.M3)+matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,nrow(x))
        if(!is.null(cova))
        {temp.mu1[,,q1]=temp.mu1[,,q1]+med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)
        temp.mu3[,,q1]=temp.mu3[,,q1]+med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)}}
        temp.mu1=exp(temp.mu1)
        temp.mu3=exp(temp.mu3)
        temp.mu1.sum=apply(temp.mu1,c(1,2),sum)
        temp.mu3.sum=apply(temp.mu3,c(1,2),sum)
        for(q1 in 1:(caty-1))
        {ie2[,,contm[1],l,q1]=(med0$BUGSoutput$sims.list$alpha1.a[,1]/deltam[contm[1]])*(temp.mu3[,,q1+1]/temp.mu3.sum-temp.mu1[,,q1+1]/temp.mu1.sum)
        aie2[,contm[1],l,q1]=apply(ie2[,,contm[1],l,q1],1,mean,na.rm=T)}
        }
        else
          for (j in 1:p1){
            temp.M3=M1
            temp.M3[,contm1[j,1]:contm1[j,2]]=M3[,contm3[j,1]:contm3[j,2]]
            temp.mu1<-array(0,c(N1,N,caty))
            temp.mu3<-temp.mu1
            for(q1 in 2:caty)
            {temp.mu1[,,q1]=med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x)+med0$BUGSoutput$sims.list$beta[,q1-1,]%*%t(M1)+matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,nrow(x))
            temp.mu3[,,q1]=med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x)+med0$BUGSoutput$sims.list$beta[,q1-1,]%*%t(temp.M3)+matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,nrow(x))
            if(!is.null(cova))
            {temp.mu1[,,q1]=temp.mu1[,,q1]+med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)
            temp.mu3[,,q1]=temp.mu3[,,q1]+med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)}}
            temp.mu1=exp(temp.mu1)
            temp.mu3=exp(temp.mu3)
            temp.mu1.sum=apply(temp.mu1,c(1,2),sum)
            temp.mu3.sum=apply(temp.mu3,c(1,2),sum)
            for(q1 in 1:(caty-1))
            {ie2[,,contm[j],l,q1]=(med0$BUGSoutput$sims.list$alpha1.a[,j,l]/deltam[contm[j]])*(temp.mu3[,,q1+1]/temp.mu3.sum-temp.mu1[,,q1+1]/temp.mu1.sum)
            aie2[,contm[j],l,q1]=apply(ie2[,,contm[j],l,q1],1,mean,na.rm=T)}
          }}


      #for binary and categorical mediators, method 2 and method 1 are not the same
      if(p2>0){
        if(p2==1 & c1==1){
          temp.M3=M1
          temp.M1=M1
          temp.M3[,binm1[1]]=1
          temp.M1[,binm1[1]]=0
          temp.mu1<-array(0,c(N1,N,caty))
          temp.mu3<-temp.mu1
          for(q1 in 2:caty)
          {temp.mu1[,,q1]=med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x)+med0$BUGSoutput$sims.list$beta[,q1-1,]%*%t(temp.M1)+matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,nrow(x))
          temp.mu3[,,q1]=med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x)+med0$BUGSoutput$sims.list$beta[,q1-1,]%*%t(temp.M3)+matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,nrow(x))
          if(!is.null(cova))
          {temp.mu1[,,q1]=temp.mu1[,,q1]+med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)
          temp.mu3[,,q1]=temp.mu3[,,q1]+med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)}}
          temp.mu1=exp(temp.mu1)
          temp.mu3=exp(temp.mu3)
          temp.mu1.sum=apply(temp.mu1,c(1,2),sum)
          temp.mu3.sum=apply(temp.mu3,c(1,2),sum)
          bpart=array(0,c(N1,N,caty))
          for(q1 in 1:(caty-1))
            bpart[,,q1]=temp.mu3[,,q1+1]/temp.mu3.sum-temp.mu1[,,q1+1]/temp.mu1.sum
          if (nmc==1){
            temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,1]+matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x1.temp)
            temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,1]+matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x3.temp)
          }
          else
          {temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,1,mind[binm[1],]]%*%t(mcov[,mind[binm[1],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x1.temp)
          temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,1,mind[binm[1],]]%*%t(mcov[,mind[binm[1],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x3.temp)}
          for(q1 in 1:(caty-1)){
            ie2[,,binm[1],1,q1]=bpart[,,q1]*(exp(temp.x3)/(1+exp(temp.x3))-exp(temp.x1)/(1+exp(temp.x1)))/deltap[l]
            aie2[,binm[1],1,q1] <-apply(ie2[,,binm[1],1,q1],1,mean,na.rm=T)}}
        else
          for (k in 1:p2){
            temp.M3=M1
            temp.M1=M1
            temp.M3[,binm1[k]]=1
            temp.M1[,binm1[k]]=0
            temp.mu1<-array(0,c(N1,N,caty))
            temp.mu3<-temp.mu1
            for(q1 in 2:caty)
            {temp.mu1[,,q1]=med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x)+med0$BUGSoutput$sims.list$beta[,q1-1,]%*%t(temp.M1)+matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,nrow(x))
            temp.mu3[,,q1]=med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x)+med0$BUGSoutput$sims.list$beta[,q1-1,]%*%t(temp.M3)+matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,nrow(x))
            if(!is.null(cova))
            {temp.mu1[,,q1]=temp.mu1[,,q1]+med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)
            temp.mu3[,,q1]=temp.mu3[,,q1]+med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)}}
            temp.mu1=exp(temp.mu1)
            temp.mu3=exp(temp.mu3)
            temp.mu1.sum=apply(temp.mu1,c(1,2),sum)
            temp.mu3.sum=apply(temp.mu3,c(1,2),sum)
            bpart=array(0,c(N1,N,caty-1))
            for(q1 in 1:(caty-1))
              bpart[,,q1]=temp.mu3[,,q1+1]/temp.mu3.sum-temp.mu1[,,q1+1]/temp.mu1.sum

            if(nmc==1 & p2==1){
              temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,k]+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x1.temp)
              temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,k]+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x3.temp)
            }
            else{
              temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,k,mind[binm[k],]]%*%t(mcov[,mind[binm[k],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x1.temp)
              temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,k,mind[binm[k],]]%*%t(mcov[,mind[binm[k],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x3.temp)
            }
            for(q1 in 1:(caty-1)){
              ie2[,,binm[k],l,q1]=bpart[,,q1]*(exp(temp.x3)/(1+exp(temp.x3))-exp(temp.x1)/(1+exp(temp.x1)))/deltap[l]
              aie2[,binm[k],l,q1] <- apply(ie2[,,binm[k],l,q1],1,mean,na.rm=T)}
          }
      }

      if(p3>0){
        for (j in 1:p3){
          bpart=array(0,c(N1,N,cat2[j]-1,caty-1))
          M1.temp=M1
          M1.temp[,catm1[j,1]:catm1[j,2]]=0
          for (k in catm1[j,1]:catm1[j,2]){
            temp.M3=M1.temp
            temp.M1=M1.temp
            temp.M3[,k]=1
            temp.mu1<-array(0,c(N1,N,caty))
            temp.mu3<-temp.mu1
            for(q1 in 2:caty)
            {temp.mu1[,,q1]=med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x)+med0$BUGSoutput$sims.list$beta[,q1-1,]%*%t(temp.M1)+matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,nrow(x))
            temp.mu3[,,q1]=med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x)+med0$BUGSoutput$sims.list$beta[,q1-1,]%*%t(temp.M3)+matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,nrow(x))
            if(!is.null(cova))
            {temp.mu1[,,q1]=temp.mu1[,,q1]+med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)
            temp.mu3[,,q1]=temp.mu3[,,q1]+med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)}}
            temp.mu1=exp(temp.mu1)
            temp.mu3=exp(temp.mu3)
            temp.mu1.sum=apply(temp.mu1,c(1,2),sum)
            temp.mu3.sum=apply(temp.mu3,c(1,2),sum)
            for(q1 in 1:(caty-1))
              bpart[,,k-catm1[j,1]+1,q1]=exp(temp.mu3)/(1+exp(temp.mu3))-exp(temp.mu1)/(1+exp(temp.mu1))
          }

          mu_Mc1<-array(0,c(N1,N,cat2[j]-1))
          mu_Mc0<-array(0,c(N1,N,cat2[j]-1))
          for (k in 2:cat2[j]){
            mu_Mc1[,,k-1] <- exp(as.matrix(med0$BUGSoutput$sims.list$alpha0.c[,j,k-1,mind[catm[j],]])%*%t(mcov[,mind[catm[j],]])+med0$BUGSoutput$sims.list$alpha1.c[,j,k-1,]%*%t(x3.temp))
            mu_Mc0[,,k-1] <- exp(as.matrix(med0$BUGSoutput$sims.list$alpha0.c[,j,k-1,mind[catm[j],]])%*%t(mcov[,mind[catm[j],]])+med0$BUGSoutput$sims.list$alpha1.c[,j,k-1,]%*%t(x1.temp))
          }
          sum_Mc1 <-apply(mu_Mc1,c(1,2),sum)+1
          sum_Mc0 <-apply(mu_Mc0,c(1,2),sum)+1
          for(q1 in 1:(caty-1))
          {for (k in 2:cat2[j])
            ie2[,,catm[j],l,q1]=ie2[,,catm[j],l,q1]+(mu_Mc1[,,k-1]/sum_Mc1-mu_Mc0[,,k-1]/sum_Mc0)*bpart[,,k-1,q1]
          aie2[,catm[j],l,q1]<-apply(ie2[,,catm[j],l,q1],1,mean,na.rm=T)}
        }
      }

      temp.x1=x
      temp.x3=x
      if(l%in%data0$contpred1.0)
      {tt1=match(l,data0$contpred1.0)
      temp.x3[,data0$contpred2[tt1,1]:data0$contpred2[tt1,2]]=data0$pred3[,data0$contpred3[tt1,1]:data0$contpred3[tt1,2]]}
      else if(l%in%data0$binpred1.0)
      {tt1=match(l,data0$binpred1)
      temp.x1[,data0$binpred2[tt1]]=0
      temp.x3[,data0$binpred2[tt1]]=1
      deltap[l]=1}
      else
      {for (i in 1:nrow(data0$catpred1.0))
        if(l%in%(data0$catpred1.0[i,1]:data0$catpred1.0[i,2]))
          tt1=i
      temp.x1[,data0$catpred2[tt1,1]:data0$catpred2[tt1,2]]=0
      temp.x3[,data0$catpred2[tt1,1]:data0$catpred2[tt1,2]]=0
      temp.x3[,l]=1
      deltap[l]=1}
      temp.mu1<-array(0,c(N1,N,caty))
      temp.mu3<-temp.mu1
      for(q1 in 2:caty)
      {temp.mu1[,,q1]=med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(temp.x1)+med0$BUGSoutput$sims.list$beta[,q1-1,]%*%t(M1)+matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,nrow(x))
      temp.mu3[,,q1]=med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(temp.x3)+med0$BUGSoutput$sims.list$beta[,q1-1,]%*%t(M1)+matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,nrow(x))
      if(!is.null(cova))
      {temp.mu1[,,q1]=temp.mu1[,,q1]+med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)
      temp.mu3[,,q1]=temp.mu3[,,q1]+med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)}}
      temp.mu1=exp(temp.mu1)
      temp.mu3=exp(temp.mu3)
      temp.mu1.sum=apply(temp.mu1,c(1,2),sum)
      temp.mu3.sum=apply(temp.mu3,c(1,2),sum)
      for(q1 in 1:(caty-1))
      {de2.1=(temp.mu3[,,q1+1]/temp.mu3.sum-temp.mu1[,,q1+1]/temp.mu1.sum)/deltap[l]
      de2[,l,q1]=apply(de2.1,1,mean,na.rm=T)}

      #method 3:parametric
      #3.1. get M1(x) and M1(x+dx) for y
      if(p1>0){
        if(c1==1 & p1+p2+p3==1 & contm1[1,1]==contm1[1,2]){
          if(nmc==1)
          {mu_M2[,contm1[1,1],] <- med0$BUGSoutput$sims.list$alpha0[,contm1[1,1]]+as.matrix(med0$BUGSoutput$sims.list$alpha1[,contm1[1,1]])%*%t(x1.temp)
          mu_M3[,contm1[1,1],] <- med0$BUGSoutput$sims.list$alpha0[,contm1[1,1]]+as.matrix(med0$BUGSoutput$sims.list$alpha1[,contm1[1,1]])%*%t(x3.temp)}
          else
          {mu_M2[,contm1[1,1],] <- med0$BUGSoutput$sims.list$alpha0[,contm1[1,1],mind[contm[1],]]%*%t(mcov[,mind[contm[1],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1[,contm1[1,1]])%*%t(x1.temp)
          mu_M3[,contm1[1,1],] <- med0$BUGSoutput$sims.list$alpha0[,contm1[1,1],mind[contm[1],]]%*%t(mcov[,mind[contm[1],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1[,contm1[1,1]])%*%t(x3.temp)}
        }
        else
          for (j in 1:p1){
            for (k in contm1[j,1]:contm1[j,2]){
              if(nmc==1 & p1+p2+p3==1)
              {mu_M2[,k,] <- med0$BUGSoutput$sims.list$alpha0[,k]+med0$BUGSoutput$sims.list$alpha1[,k,l]%*%t(x1.temp)
              mu_M3[,k,] <- med0$BUGSoutput$sims.list$alpha0[,k]+med0$BUGSoutput$sims.list$alpha1[,k,l]%*%t(x3.temp)}
              else
              {mu_M2[,k,] <- med0$BUGSoutput$sims.list$alpha0[,k,mind[contm[j],]]%*%t(mcov[,mind[contm[j],]])+med0$BUGSoutput$sims.list$alpha1[,k,]%*%t(x1.temp)
              mu_M3[,k,] <- med0$BUGSoutput$sims.list$alpha0[,k,mind[contm[j],]]%*%t(mcov[,mind[contm[j],]])+med0$BUGSoutput$sims.list$alpha1[,k,]%*%t(x3.temp)}
            }}}

      if(p2>0){
        if(p2==1 & c1==1)
        {if(nmc==1)
        {temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,1]+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x1.temp)
        temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,1]+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x3.temp)}
          else
          {temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,1,mind[binm[1],]]%*%t(mcov[,mind[binm[1],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x1.temp)
          temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,1,mind[binm[1],]]%*%t(mcov[,mind[binm[1],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x3.temp)}
          mu_M2[,binm[1],] <- exp(temp.x1)/(1+exp(temp.x1))
          mu_M3[,binm[1],] <- exp(temp.x3)/(1+exp(temp.x3))
        }
        else{
          for (k in 1:p2){
            if(nmc==1 & p2==1)
            {temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,k]+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x1.temp)
            temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,k]+med0$BUGSoutput$sims.list$alpha1.b[,k,]%*%t(x3.temp)}
            else
            {temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,k,mind[binm[k],]]%*%t(mcov[,mind[binm[k],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x1.temp)
            temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,k,mind[binm[k],]]%*%t(mcov[,mind[binm[k],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x3.temp)}
            mu_M2[,binm1[k],] <- exp(temp.x1)/(1+exp(temp.x1))
            mu_M3[,binm1[k],] <- exp(temp.x3)/(1+exp(temp.x3))}
        }}

      if(p3>0){
        for (j in 1:p3){
          mu_Mc1<-array(0,c(N1,N,cat2[j]-1))
          mu_Mc0<-array(0,c(N1,N,cat2[j]-1))
          for (k in 2:cat2[j]){
            mu_Mc1[,,k-1] <- exp(as.matrix(med0$BUGSoutput$sims.list$alpha0.c[,j,k-1,mind[catm[j],]])%*%t(mcov[,mind[catm[j],]])+med0$BUGSoutput$sims.list$alpha1.c[,j,k-1,]%*%t(x3.temp))
            mu_Mc0[,,k-1] <- exp(as.matrix(med0$BUGSoutput$sims.list$alpha0.c[,j,k-1,mind[catm[j],]])%*%t(mcov[,mind[catm[j],]])+med0$BUGSoutput$sims.list$alpha1.c[,j,k-1,]%*%t(x1.temp))
          }
          sum_Mc1 <-apply(mu_Mc1,c(1,2),sum)+1
          sum_Mc0 <-apply(mu_Mc0,c(1,2),sum)+1
          for (k in 2:cat2[j])
          {mu_M2[,catm1[j,1]+k-2,]=mu_Mc0[,,k-1]/sum_Mc0
          mu_M3[,catm1[j,1]+k-2,]=mu_Mc1[,,k-1]/sum_Mc1}
        }}

      #3.2. get x and dx for y
      x1.temp1=x
      x3.temp1=x
      if(l%in%as.vector(data0$contpred1.0)){ #need the continuous predictor in its original format
        i=match(l,data0$contpred1.0)
        x3.temp1[,data0$contpred2[i,1]:data0$contpred2[i,2]]=data0$pred3[,data0$contpred3[i,1]:data0$contpred3[i,2]]
      }
      else if(l%in%data0$binpred1.0)
      {i=match(l,data0$binpred1.0)
      x1.temp1[,data0$binpred2[i]]=0
      x3.temp1[,data0$binpred2[i]]=1
      deltap[l]=1}
      else{ #categorical predictor
        for (i in 1:nrow(data0$catpred1.0))
          if(l%in%(data0$catpred1.0[i,1]:data0$catpred1.0[i,2]))
          {x1.temp1[,data0$catpred2[i,1]:data0$catpred2[i,2]]=0
          x3.temp1[,data0$catpred2[i,1]:data0$catpred2[i,2]]=0
          di=match(l,data0$catpred1.0[i,1]:data0$catpred1.0[i,2])
          x3.temp1[,data0$catpred2[i,1]+di-1]=1}
        deltap[l]=1
      }

      #3.3. get the total effect
      mu_y0<-array(0,c(N1,N,caty))
      mu_y1<-mu_y0
      for(q1 in 2:caty)
      {temp1=array(0,c(nrow(M1)+1,ncol(M1),N1))
      temp1[2:(nrow(M1)+1),,]=mu_M2
      temp1[1,,]=t(med0$BUGSoutput$sims.list$beta[,q1-1,])
      temp2=temp1
      temp2[2:(nrow(M1)+1),,]=mu_M3
      if(is.null(cova))
      {mu_y0[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x1.temp1) + t(apply(temp1,3,matrix.prod))
      mu_y1[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x3.temp1) + t(apply(temp2,3,matrix.prod))}
      else
      {mu_y0[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x1.temp1) + med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(as.matrix(cova))+t(apply(temp1,3,matrix.prod))
      mu_y1[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x3.temp1) + med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(as.matrix(cova))+t(apply(temp2,3,matrix.prod))
      }} #get the linear part

      mu_y0=exp(mu_y0)
      mu_y1=exp(mu_y1)
      mu_y0.sum=apply(mu_y0,c(1,2),sum)
      mu_y1.sum=apply(mu_y1,c(1,2),sum)
      for(q1 in 1:(caty-1)){
        te3[,,q1]=(mu_y1[,,q1+1]/mu_y1.sum-mu_y0[,,q1+1]/mu_y0.sum)/deltap[l]
        ate3[,l,q1]=apply(te3[,,q1],1,mean,na.rm=T)}

      #3.4. calculate the ie
      j1=sample(1:N,size=N*N1,replace=T)
      j2=sample(1:N,size=N*N1,replace=T)

      #3.4.1. continuous mediators
      if(p1>0){
        for (j in 1:p1){
          mu_y0.2<-array(0,c(N1,N,caty))
          mu_y1.2<-mu_y0.2

          for(q1 in 2:caty)
          {temp1.1=temp1
          temp1.2=temp2
          temp1.1[1,,]=t(med0$BUGSoutput$sims.list$beta[,q1-1,])
          temp1.2[1,,]=t(med0$BUGSoutput$sims.list$beta[,q1-1,])

          for (i in contm1[j,1]:contm1[j,2])
          {temp1.1[2:(nrow(M1)+1),i,]=matrix(M1[j1,i],N,N1)
          temp1.2[2:(nrow(M1)+1),i,]=matrix(M1[j1,i],N,N1)} #j2/j1
          if(is.null(cova))
          {mu_y0.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x1.temp1) +t(apply(temp1.1,3,matrix.prod))
          mu_y1.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x3.temp1) +t(apply(temp1.2,3,matrix.prod)) }
          else{
            mu_y0.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x1.temp1) + med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)+t(apply(temp1.1,3,matrix.prod))
            mu_y1.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x3.temp1) + med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)+t(apply(temp1.2,3,matrix.prod))}
          }
          mu_y0.2=exp(mu_y0.2)
          mu_y1.2=exp(mu_y1.2)
          mu_y0.2.sum=apply(mu_y0.2,c(1,2),sum)
          mu_y1.2.sum=apply(mu_y1.2,c(1,2),sum)

          for(q1 in 1:(caty-1))
            ie3[,,contm[j],l,q1]<-te3[,,q1]-(mu_y1.2[,,q1+1]/mu_y1.2.sum-mu_y0.2[,,q1+1]/mu_y0.2.sum)/deltap[l]
        }}

      #3.4.2. binary mediators
      if(p2>0){
        for (k in 1:p2){
          mu_y0.2<-array(0,c(N1,N,caty))
          mu_y1.2<-mu_y0.2

          for(q1 in 2:caty)
          {temp1.1=temp1
          temp1.2=temp2
          temp1.1[1,,]=t(med0$BUGSoutput$sims.list$beta[,q1-1,])
          temp1.2[1,,]=t(med0$BUGSoutput$sims.list$beta[,q1-1,])
          temp1.1[2:(nrow(M1)+1),binm1[k],]=matrix(M1[j1,binm1[k]],N,N1)
          temp1.2[2:(nrow(M1)+1),binm1[k],]=matrix(M1[j1,binm1[k]],N,N1) #j2/j1
          if (is.null(cova)){
            mu_y0.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x1.temp1) + t(apply(temp1.1,3,matrix.prod))
            mu_y1.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x3.temp1) + t(apply(temp1.2,3,matrix.prod))}
          else{
            mu_y0.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x1.temp1) + med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)+t(apply(temp1.1,3,matrix.prod))
            mu_y1.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x3.temp1) + med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)+t(apply(temp1.2,3,matrix.prod))}
          }
          mu_y0.2=exp(mu_y0.2)
          mu_y1.2=exp(mu_y1.2)
          mu_y0.2.sum=apply(mu_y0.2,c(1,2),sum)
          mu_y1.2.sum=apply(mu_y1.2,c(1,2),sum)

          for(q1 in 1:(caty-1))
            ie3[,,binm[k],l,q1]<-te3[,,q1]-(mu_y1.2[,,q1+1]/mu_y1.2.sum-mu_y0.2[,,q1+1]/mu_y0.2.sum)/deltap[l]
        }}


      if(p3>0){
        for (j in 1:p3){
          mu_y0.2<-array(0,c(N1,N,caty))
          mu_y1.2<-mu_y0.2

          for(q1 in 2:caty)
          {temp1.1=temp1
          temp1.2=temp2
          temp1.1[1,,]=t(med0$BUGSoutput$sims.list$beta[,q1-1,])
          temp1.2[1,,]=t(med0$BUGSoutput$sims.list$beta[,q1-1,])
          for (i in catm1[j,1]:catm1[j,2])
          {temp1.1[2:(nrow(M1)+1),i,]=matrix(M1[j1,i],N,N1)
          temp1.2[2:(nrow(M1)+1),i,]=matrix(M1[j1,i],N,N1)} #j2/j1
          if(is.null(cova)){
            mu_y0.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x1.temp1) + t(apply(temp1.1,3,matrix.prod))
            mu_y1.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x3.temp1) + t(apply(temp1.2,3,matrix.prod))}
          else{
            mu_y0.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x1.temp1) + med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)+t(apply(temp1.1,3,matrix.prod))
            mu_y1.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x3.temp1) + med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)+t(apply(temp1.2,3,matrix.prod))}
          }

          mu_y0.2=exp(mu_y0.2)
          mu_y1.2=exp(mu_y1.2)
          mu_y0.2.sum=apply(mu_y0.2,c(1,2),sum)
          mu_y1.2.sum=apply(mu_y1.2,c(1,2),sum)

          for(q1 in 1:(caty-1))
            ie3[,,catm[j],l,q1]<-te3[,,q1]-(mu_y1.2[,,q1+1]/mu_y1.2.sum-mu_y0.2[,,q1+1]/mu_y0.2.sum)/deltap[l]
        }}

      aie3[,,l,]<-apply(array(ie3[,,,l,],c(N1,N,p1+p2+p3,caty-1)),c(1,3,4),mean,na.rm=T)

      #3.5. Calculate the de
      mu_y0.2<-array(0,c(N1,N,caty))
      mu_y1.2<-mu_y0.2

      for(q1 in 2:caty)
      {temp1.1=temp1
      temp1.2=temp2
      temp1.1[1,,]=t(med0$BUGSoutput$sims.list$beta[,q1-1,])
      temp1.2[1,,]=t(med0$BUGSoutput$sims.list$beta[,q1-1,])
      for (i in 1:ncol(M1))
      {temp1.1[2:(nrow(M1)+1),i,]=matrix(M1[j1,i],N,N1)
      temp1.2[2:(nrow(M1)+1),i,]=matrix(M1[j1,i],N,N1)} #j2/j1
      if(is.null(cova)){
        mu_y0.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x1.temp1) + t(apply(temp1.1,3,matrix.prod))
        mu_y1.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x3.temp1) + t(apply(temp1.2,3,matrix.prod))}
      else{
        mu_y0.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x1.temp1) + med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)+t(apply(temp1.1,3,matrix.prod))
        mu_y1.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x3.temp1) + med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)+t(apply(temp1.2,3,matrix.prod))}
      }

      mu_y0.2=exp(mu_y0.2)
      mu_y1.2=exp(mu_y1.2)
      mu_y0.2.sum=apply(mu_y0.2,c(1,2),sum)
      mu_y1.2.sum=apply(mu_y1.2,c(1,2),sum)

      for(q1 in 1:(caty-1))
      {de3[,,q1]<-(mu_y1.2[,,q1+1]/mu_y1.2.sum-mu_y0.2[,,q1+1]/mu_y0.2.sum)/deltap[l]
      ade3[,l,q1]=apply(de3[,,q1],1,mean,na.rm=T)}

      #method3: semi-parametric for binary or categorical predictors
      if(!l%in%data0$contpred1.0){
        if(!is.null(data0$binpred1.0))
        {if (l%*%data0$binpred1.0)
        {M.0=data.frame(M1[x1[,l]==0,])
        y.0=y[x1[,l]==0]}
          else
          {for(i in 1:nrow(data0$catpred1.0))
            if(l%in%(data0$catpred1.0[i,1]:data0$catpred1.0[i,2]))
              tt1=i
          M.0=data.frame(M1[apply(x1[,data0$catpred1.0[tt1,1]:data0$catpred1.0[[tt1,2]]]==1,1,sum)==0,])
          y.0=y[apply(x1[,data0$catpred1.0[tt1,1]:data0$catpred1.0[[tt1,2]]]==1,1,sum)==0]}}
        else
        {for(i in 1:nrow(data0$catpred1.0))
          if(l%in%(data0$catpred1.0[i,1]:data0$catpred1.0[i,2]))
            tt1=i
        M.0=data.frame(M1[apply(x1[,data0$catpred1.0[tt1,1]:data0$catpred1.0[[tt1,2]]]==1,1,sum)==0,])
        y.0=y[apply(x1[,data0$catpred1.0[tt1,1]:data0$catpred1.0[[tt1,2]]]==1,1,sum)==0]}

        y.1=y[x1[,l]==1]
        M.1=data.frame(M1[x1[,l]==1,])

        j1=sample(1:N,size=N*N1,replace=T)
        j2=sample(1:N,size=N*N1,replace=T)
        n3=nrow(M.0)
        n4=nrow(M.1)
        j3=sample(1:n3,size=N*N1,replace = T)
        j4=sample(1:n4,size=N*N1,replace = T)

        for (i in 1:ncol(M1))
        {mu.M0[,i,]=matrix(M.0[j3,i],N,N1)
        mu.M1[,i,]=matrix(M.1[j4,i],N,N1)}

        #4.1. get the total effect
        mu_y0=matrix(y.0[j3],N1,N)
        mu_y1=matrix(y.1[j4],N1,N)

        temp1=array(0,c(nrow(M1)+1,ncol(M1),N1))
        temp1[2:(nrow(M1)+1),,]=mu.M0
        temp2=temp1
        temp2[2:(nrow(M1)+1),,]=mu.M1

        temp.levels=levels(y)[-1]
        for(q1 in 1:length(temp.levels))
        {te4[,,q1]<- (mu_y1==temp.levels[q1])-(mu_y0==temp.levels[q1])
        ate4[,l,q1]=apply(te4[,,q1],1,mean,na.rm=T)}

        #4.2. Get the ies
        #ie for continuous mediators
        if(p1>0){
          for (j in 1:p1){
            mu_y0.2<-array(0,c(N1,N,caty))
            mu_y1.2<-mu_y0.2

            for(q1 in 2:caty)
            {temp1.1=temp1
            temp1.2=temp2
            temp1.1[1,,]=t(med0$BUGSoutput$sims.list$beta[,q1-1,])
            temp1.2[1,,]=t(med0$BUGSoutput$sims.list$beta[,q1-1,])

            for (i in contm1[j,1]:contm1[j,2])
            {temp1.1[2:(nrow(M1)+1),i,]=matrix(M1[j1,i],N,N1)
            temp1.2[2:(nrow(M1)+1),i,]=matrix(M1[j1,i],N,N1)} #j2/j1
            if(is.null(cova))
            {mu_y0.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x1.temp1) +t(apply(temp1.1,3,matrix.prod))
            mu_y1.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x3.temp1) +t(apply(temp1.2,3,matrix.prod)) }
            else{
              mu_y0.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x1.temp1) + med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)+t(apply(temp1.1,3,matrix.prod))
              mu_y1.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x3.temp1) + med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)+t(apply(temp1.2,3,matrix.prod))}

            }
            mu_y0.2=exp(mu_y0.2)
            mu_y1.2=exp(mu_y1.2)
            mu_y0.2.sum=apply(mu_y0.2,c(1,2),sum)
            mu_y1.2.sum=apply(mu_y1.2,c(1,2),sum)

            for(q1 in 1:(caty-1))
              ie4[,,contm[j],l,q1]<-te4[,,q1]-(mu_y1.2[,,q1+1]/mu_y1.2.sum-mu_y0.2[,,q1+1]/mu_y0.2.sum)/deltap[l]
          }}

        #ie for binary mediators
        if(p2>0){
          for (k in 1:p2){
            mu_y0.2<-array(0,c(N1,N,caty))
            mu_y1.2<-mu_y0.2

            for(q1 in 2:caty)
            {temp1.1=temp1
            temp1.2=temp2
            temp1.1[2:(nrow(M1)+1),binm1[k],]=matrix(M1[j1,binm1[k]],N,N1)
            temp1.2[2:(nrow(M1)+1),binm1[k],]=matrix(M1[j1,binm1[k]],N,N1) #j2/j1
            temp1.1[1,,]=t(med0$BUGSoutput$sims.list$beta[,q1-1,])
            temp1.2[1,,]=t(med0$BUGSoutput$sims.list$beta[,q1-1,])

            if (is.null(cova)){
              mu_y0.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x1.temp1) + t(apply(temp1.1,3,matrix.prod))
              mu_y1.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x3.temp1) + t(apply(temp1.2,3,matrix.prod))}
            else{
              mu_y0.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x1.temp1) + med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)+t(apply(temp1.1,3,matrix.prod))
              mu_y1.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x3.temp1) + med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)+t(apply(temp1.2,3,matrix.prod))}
            }
            mu_y0.2=exp(mu_y0.2)
            mu_y1.2=exp(mu_y1.2)
            mu_y0.2.sum=apply(mu_y0.2,c(1,2),sum)
            mu_y1.2.sum=apply(mu_y1.2,c(1,2),sum)

            for(q1 in 1:(caty-1))
              ie4[,,binm[k],l,q1]<-te4[,,q1]-(mu_y1.2[,,q1+1]/mu_y1.2.sum-mu_y0.2[,,q1+1]/mu_y0.2.sum)
          }}

        #ie for categorical mediators
        if(p3>0){
          for (j in 1:p3){
            mu_y0.2<-array(0,c(N1,N,caty))
            mu_y1.2<-mu_y0.2

            for(q1 in 2:caty)
            {temp1.1=temp1
            temp1.2=temp2
            temp1.1[1,,]=t(med0$BUGSoutput$sims.list$beta[,q1-1,])
            temp1.2[1,,]=t(med0$BUGSoutput$sims.list$beta[,q1-1,])
            for (i in catm1[j,1]:catm1[j,2])
            {temp1.1[2:(nrow(M1)+1),i,]=matrix(M1[j1,i],N,N1)
            temp1.2[2:(nrow(M1)+1),i,]=matrix(M1[j1,i],N,N1)} #j2/j1
            if(is.null(cova)){
              mu_y0.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x1.temp1) + t(apply(temp1.1,3,matrix.prod))
              mu_y1.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x3.temp1) + t(apply(temp1.2,3,matrix.prod))}
            else{
              mu_y0.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x1.temp1) + med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)+t(apply(temp1.1,3,matrix.prod))
              mu_y1.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x3.temp1) + med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)+t(apply(temp1.2,3,matrix.prod))}
            }
            mu_y0.2=exp(mu_y0.2)
            mu_y1.2=exp(mu_y1.2)
            mu_y0.2.sum=apply(mu_y0.2,c(1,2),sum)
            mu_y1.2.sum=apply(mu_y1.2,c(1,2),sum)

            for(q1 in 1:(caty-1))
              ie4[,,catm[j],l,q1]<-te4[,,q1]-(mu_y1.2[,,q1+1]/mu_y1.2.sum-mu_y0.2[,,q1+1]/mu_y0.2.sum)
          }}

        aie4[,,l,]<-apply(array(ie4[,,,l,],c(N1,N,p1+p2+p3,caty-1)),c(1,3,4),mean,na.rm=T)

        #4.3. Calculate the de
        mu_y0.2<-array(0,c(N1,N,caty))
        mu_y1.2<-mu_y0.2

        for(q1 in 2:caty)
        {temp1.1=temp1
        temp1.2=temp2
        for (i in 1:ncol(M1))
        {temp1.1[2:(nrow(M1)+1),i,]=matrix(M1[j1,i],N,N1)
        temp1.2[2:(nrow(M1)+1),i,]=matrix(M1[j1,i],N,N1)} #j2/j1
        temp1.1[1,,]=t(med0$BUGSoutput$sims.list$beta[,q1-1,])
        temp1.2[1,,]=t(med0$BUGSoutput$sims.list$beta[,q1-1,])
        if(is.null(cova)){
          mu_y0.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x1.temp1) + t(apply(temp1.1,3,matrix.prod))
          mu_y1.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x3.temp1) + t(apply(temp1.2,3,matrix.prod))}
        else{
          mu_y0.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x1.temp1) + med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)+t(apply(temp1.1,3,matrix.prod))
          mu_y1.2[,,q1]<- matrix(med0$BUGSoutput$sims.list$beta0[,q1-1],N1,N) + med0$BUGSoutput$sims.list$c[,q1-1,]%*%t(x3.temp1) + med0$BUGSoutput$sims.list$eta[,q1-1,]%*%t(cova)+t(apply(temp1.2,3,matrix.prod))}
        }
        mu_y0.2=exp(mu_y0.2)
        mu_y1.2=exp(mu_y1.2)
        mu_y0.2.sum=apply(mu_y0.2,c(1,2),sum)
        mu_y1.2.sum=apply(mu_y1.2,c(1,2),sum)

      }
    }

    ate1=apply(aie1,c(1,3,4),sum)+de1
    ate2=apply(aie2,c(1,3,4),sum)+de2

    colnames(aie4)=colnames(M2)
    colnames(aie1)=colnames(M2)
    colnames(aie2)=colnames(M2)
    colnames(aie3)=colnames(M2)
  }
  else{
    aie1=array(0,c(N1,p1+p2+p3,c1))
    ie1=array(0,dim=c(N1,N,p1+p2+p3,c1))
    tt2=1
    de1=NULL
    n.cont=1

    aie2=array(0,c(N1,p1+p2+p3,c1))
    ie2=array(0,dim=c(N1,N,p1+p2+p3,c1))
    de2=NULL

    ate3=matrix(0,N1,c1)
    omu3=ate3
    ade3=matrix(0,N1,c1)
    aie3=array(0,c(N1,p1+p2+p3,c1))
    ie3=array(0,dim=c(N1,N,p1+p2+p3,c1))
    mu_M2=array(0,c(dim(M1),N1))
    mu_M3=mu_M2

    ade4=matrix(0,N1,c1)
    ie4=array(0,dim=c(N1,N,p1+p2+p3,c1))
    mu.M0=array(0,c(dim(M1),N1))
    mu.M1=mu.M0
    ate4=matrix(0,N1,c1)
    omu4=ate4
    aie4=array(0,c(N1,p1+p2+p3,c1))

    for (l in 1:c1){
      x1.temp=x1
      x3.temp=x1
      if(l%in%data0$contpred1.0){
        x3.temp[,l]=x1[,l]+deltap[l]
      }
      else if(l%in%data0$binpred1.0)
      {x1.temp[,l]=0
      x3.temp[,l]=1
      deltap[l]=1}
      else{ #categorical predictor
        for (i in 1:nrow(data0$catpred1.0))
          if(l%in%(data0$catpred1.0[i,1]:data0$catpred1.0[i,2]))
          {x1.temp[,data0$catpred1.0[i,1]:data0$catpred1.0[i,2]]=0
          x3.temp[,data0$catpred1.0[i,1]:data0$catpred1.0[i,2]]=0
          x3.temp[,l]=1}
        deltap[l]=1
      }

      #method 1: the same for binary or continuous predictors
      if(p1>0){
        if (p1+p2+p3==1 & c1==1 & contm1[1,1]==contm1[1,2])
          aie1[,contm[1],l]=med0$BUGSoutput$sims.list$alpha1[,1]*med0$BUGSoutput$sims.list$beta[,contm1[1,1]]#*mean(data0$m.cont.der[,contm3[1,1]]) since alpha1 is the coefficient of x on the change of f(M), beta is the coefficient of f(M) on the change of y
        else
          for (j in 1: p1)
            aie1[,contm[j],l]=apply(as.matrix(med0$BUGSoutput$sims.list$alpha1[,contm1[j,1]:contm1[j,2],l])*as.matrix(med0$BUGSoutput$sims.list$beta[,contm1[j,1]:contm1[j,2]]),1,sum)#*mean(data0$m.cont.der[,contm3[j,1]])
      }


      if(p2>0){
        if(p2==1 & c1==1) #since c1=1, the predictor can only be binary or continuous
        {if (nmc==1)
        {temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,1]+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x1.temp)
        temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,1]+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x3.temp)}
          else
          {temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,1,mind[binm[1],]]%*%t(mcov[,mind[binm[1],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x1.temp)
          temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,1,mind[binm[1],]]%*%t(mcov[,mind[binm[1],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x3.temp)}
          if(!1%in%data0$contpred1.0) #if the predictor is binary
            ie1[,,binm[1],1]=med0$BUGSoutput$sims.list$beta[,binm1[1]]*(exp(temp.x3)/(1+exp(temp.x3))-exp(temp.x1)/(1+exp(temp.x1)))
          else #if the predictor is continuous
            ie1[,,binm[1],1]=med0$BUGSoutput$sims.list$alpha1.b[,1]*med0$BUGSoutput$sims.list$beta[,binm1[1]]*exp(temp.x1)/(1+exp(temp.x1))^2
          aie1[,binm[1],l] <-apply(ie1[,,binm[1],1],1,mean)
        }
        else
          for (k in 1:p2){
            if (nmc==1 & p2==1)
            {temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,k]+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x1.temp)
            temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,k]+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x3.temp)
            }
            else
            {temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,k,mind[binm[k],]]%*%t(mcov[,mind[binm[k],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x1.temp)
            temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,k,mind[binm[k],]]%*%t(mcov[,mind[binm[k],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x3.temp)
            }
            if(!l%in%data0$contpred1.0) #if the predictor is binary or categorical
              ie1[,,binm[k],l]=med0$BUGSoutput$sims.list$beta[,binm1[k]]*(exp(temp.x3)/(1+exp(temp.x3))-exp(temp.x1)/(1+exp(temp.x1)))
            else if(data0$contpred1.0) #for continuous predictor
              ie1[,,binm[1],l]=med0$BUGSoutput$sims.list$alpha1.b[,k,l]*med0$BUGSoutput$sims.list$beta[,binm1[k]]*exp(temp.x1)/(1+exp(temp.x1))^2
            aie1[,binm[k],l] <- apply(ie1[,,binm[k],l],1,mean)
          }
      }

      if(p3>0){
        for (j in 1:p3){
          mu_Mc1<-array(0,c(N1,N,cat2[j]-1))
          mu_Mc0<-array(0,c(N1,N,cat2[j]-1))
          for (k in 2:cat2[j]){
            mu_Mc0[,,k-1] <- exp(as.matrix(med0$BUGSoutput$sims.list$alpha0.c[,j,k-1,mind[catm[j],]])%*%t(mcov[,mind[catm[j],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.c[,j,k-1,])%*%t(x1.temp))
            mu_Mc1[,,k-1] <- exp(as.matrix(med0$BUGSoutput$sims.list$alpha0.c[,j,k-1,mind[catm[j],]])%*%t(mcov[,mind[catm[j],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.c[,j,k-1,])%*%t(x3.temp))
          }
          sum_Mc1 <-apply(mu_Mc1,c(1,2),sum)+1
          sum_Mc0 <-apply(mu_Mc0,c(1,2),sum)+1
          if(l%in%data0$contpred1.0) #for continuous predictor
          {tt.0=rep(0,N1)
          for (k in 2:cat2[j])
            tt.0=mu_Mc0[,,k-1]/sum_Mc0*med0$BUGSoutput$sims.list$alpha1.c[,j,k-1,l]
          for (k in 2:cat2[j])
            ie1[,,catm[j],l]=ie1[,,catm[j],l]+(mu_Mc0[,,k-1]/sum_Mc0)*med0$BUGSoutput$sims.list$beta[,catm1[j,1]+k-2]*(med0$BUGSoutput$sims.list$alpha1.c[,j,k-1,l]-tt.0)
          }
          else #for binary or categorical predictor
            for (k in 2:cat2[j])
              ie1[,,catm[j],l]=ie1[,,catm[j],l]+(mu_Mc1[,,k-1]/sum_Mc1-mu_Mc0[,,k-1]/sum_Mc0)*med0$BUGSoutput$sims.list$beta[,catm1[j,1]+k-2]
          aie1[,catm[j],l]<-apply(ie1[,,catm[j],l],1,mean)
        }}

      if(l%in%data0$contpred1)
      {tt1=match(l,data0$contpred1)
      de1=cbind(de1,as.matrix(med0$BUGSoutput$sims.list$c[,data0$contpred3[tt1,1]:data0$contpred3[tt1,2]])%*%
                  apply(as.matrix(data0$pred.cont.der[,data0$contpred3[tt1,1]:data0$contpred3[tt1,2]]),2,mean))
      tt2=tt2+data0$contpred3[tt1,2]-data0$contpred3[tt1,1]+1
      }
      else
      {de1=cbind(de1,med0$BUGSoutput$sims.list$c[,tt2])
      tt2=tt2+1}

      #method2: the same for binary and continuous predictors
      if(p1>0){
        if(y.type==1){
          if(p1==1 & c1==1 & contm1[1,1]==contm1[1,2])
          { ie2[,,contm[1],l]=(med0$BUGSoutput$sims.list$alpha1.a[,1]/deltam[1])*   #
            (as.matrix(med0$BUGSoutput$sims.list$beta[,contm1[1,1]])%*%(M3[,contm3[1,1]]-M1[,contm1[1,1]]))
          aie2[,contm[1],l]=apply(ie2[,,contm[1],l],1,mean,na.rm=T)
          }
          else
            for (j in 1:p1){
              ie2[,,contm[j],l]=(med0$BUGSoutput$sims.list$alpha1.a[,j,l]/deltam[j])*   #a unit change in x, result in the unit change in m
                (med0$BUGSoutput$sims.list$beta[,contm1[j,1]:contm1[j,2]]%*%t(M3[,contm3[j,1]:contm3[j,2]]-M1[,contm1[j,1]:contm1[j,2]]))
              aie2[,contm[j],l]=apply(ie2[,,contm[j],l],1,mean,na.rm=T)  #unit change in m, result in the unit change in y
            }}
        else if(y.type==2)
        {if(p1==1 & c1==1 & contm1[1,1]==contm1[1,2])
        { temp.M3=M1
        temp.M3[,contm1[1,1]]=M3[,contm3[1,1]]
        temp.mu1=med0$BUGSoutput$sims.list$c%*%t(x)+med0$BUGSoutput$sims.list$beta%*%t(M1)+matrix(med0$BUGSoutput$sims.list$beta0,N1,nrow(x))
        temp.mu3=med0$BUGSoutput$sims.list$c%*%t(x)+med0$BUGSoutput$sims.list$beta%*%t(temp.M3)+matrix(med0$BUGSoutput$sims.list$beta0,N1,nrow(x))
        if(!is.null(cova))
        {temp.mu1=temp.mu1+med0$BUGSoutput$sims.list$eta%*%t(cova)
        temp.mu3=temp.mu3+med0$BUGSoutput$sims.list$eta%*%t(cova)}
        ie2[,,contm[1],l]=(med0$BUGSoutput$sims.list$alpha1.a[,1]/deltam[1])*(exp(temp.mu3)/(1+exp(temp.mu3))-exp(temp.mu1)/(1+exp(temp.mu1)))
        aie2[,contm[1],l]=apply(ie2[,,contm[1],l],1,mean,na.rm=T)
        }
          else
            for (j in 1:p1){
              temp.M3=M1
              temp.M3[,contm1[j,1]:contm1[j,2]]=M3[,contm3[j,1]:contm3[j,2]]
              temp.mu1=med0$BUGSoutput$sims.list$c%*%t(x)+med0$BUGSoutput$sims.list$beta%*%t(M1)+matrix(med0$BUGSoutput$sims.list$beta0,N1,nrow(x))
              temp.mu3=med0$BUGSoutput$sims.list$c%*%t(x)+med0$BUGSoutput$sims.list$beta%*%t(temp.M3)+matrix(med0$BUGSoutput$sims.list$beta0,N1,nrow(x))
              if(!is.null(cova))
              {temp.mu1=temp.mu1+med0$BUGSoutput$sims.list$eta%*%t(cova)
              temp.mu3=temp.mu3+med0$BUGSoutput$sims.list$eta%*%t(cova)}
              ie2[,,contm[j],l]=(med0$BUGSoutput$sims.list$alpha1.a[,j,l]/deltam[j])*(exp(temp.mu3)/(1+exp(temp.mu3))-exp(temp.mu1)/(1+exp(temp.mu1)))
              aie2[,contm[j],l]=apply(ie2[,,contm[j],l],1,mean,na.rm=T)
            }}
        else if(y.type==4)
        {if(p1==1 & c1==1 & contm1[1,1]==contm1[1,2])
        {temp.M3=M1
        temp.M3[,contm1[1,1]]=M3[,contm3[1,1]]
        temp.mu1=med0$BUGSoutput$sims.list$c%*%t(x)+med0$BUGSoutput$sims.list$beta%*%t(M1)
        temp.mu3=med0$BUGSoutput$sims.list$c%*%t(x)+med0$BUGSoutput$sims.list$beta%*%t(temp.M3)
        if(!is.null(cova))
        {temp.mu1=temp.mu1+med0$BUGSoutput$sims.list$eta%*%t(cova)
        temp.mu3=temp.mu3+med0$BUGSoutput$sims.list$eta%*%t(cova)}
        vec=cbind(med0$BUGSoutput$sims.list$r,
                  as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(temp.mu3),2,expp, 1/med0$BUGSoutput$sims.list$r)))
        tmean3=t(apply(vec,1,weib.trunc.mean,right=tmax))
        #if(multi)
        #  tmean3=log(tmean3)
        vec=cbind(med0$BUGSoutput$sims.list$r,
                  as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(temp.mu1),2,expp, 1/med0$BUGSoutput$sims.list$r)))
        tmean1=t(apply(vec,1,weib.trunc.mean,right=tmax))
        #if(multi)
        #  tmean1=log(tmean1)
        #tmean3=ifelse(tmean3==Inf,tmax-runif(1,0,0.1),tmean3)
        #tmean1=ifelse(tmean1==Inf,tmax-runif(1,0,0.1),tmean1)
        ie2[,,contm[1],l]=as.vector(med0$BUGSoutput$sims.list$alpha1.a[,1]/deltam[1])*(tmean3-tmean1)
        aie2[,contm[1],l]=apply(ie2[,,contm[1],l],1,mean,na.rm=T)
        }
          else
            for (j in 1:p1){
              temp.M3=M1
              temp.M3[,contm1[j,1]:contm1[j,2]]=M3[,contm3[j,1]:contm3[j,2]]
              temp.mu1=med0$BUGSoutput$sims.list$c%*%t(x)+med0$BUGSoutput$sims.list$beta%*%t(M1)
              temp.mu3=med0$BUGSoutput$sims.list$c%*%t(x)+med0$BUGSoutput$sims.list$beta%*%t(temp.M3)
              if(!is.null(cova))
              {temp.mu1=temp.mu1+med0$BUGSoutput$sims.list$eta%*%t(cova)
              temp.mu3=temp.mu3+med0$BUGSoutput$sims.list$eta%*%t(cova)}
              vec=cbind(med0$BUGSoutput$sims.list$r,
                        as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(temp.mu3),2,expp, 1/med0$BUGSoutput$sims.list$r)))
              tmean3=t(apply(vec,1,weib.trunc.mean,right=tmax))
              #if(multi)
              #  tmean3=log(tmean3)
              vec=cbind(med0$BUGSoutput$sims.list$r,
                        as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(temp.mu1),2,expp, 1/med0$BUGSoutput$sims.list$r)))
              tmean1=t(apply(vec,1,weib.trunc.mean,right=tmax))
              #if(multi)
              #  tmean1=log(tmean1)
              #tmean3=ifelse(tmean3==Inf,tmax-runif(1,0,0.1),tmean3)
              #tmean1=ifelse(tmean1==Inf,tmax-runif(1,0,0.1),tmean1)
              ie2[,,contm[j],l]=as.vector(med0$BUGSoutput$sims.list$alpha1.a[,j,l]/deltam[j])*(tmean3-tmean1)
              aie2[,contm[j],l]=apply(ie2[,,contm[j],l],1,mean)
            }}
      }

      #for binary and categorical mediators, method 2 and method 1 are not the same
      if(p2>0){
        if(y.type==1){
          if(p2==1 & c1==1){
            if (nmc==1){
              temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,1]+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x1.temp)
              temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,1]+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x3.temp)
            }
            else
            {temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,1,mind[binm[1],]]%*%t(mcov[,mind[binm[1],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x1.temp)
            temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,1,mind[binm[1],]]%*%t(mcov[,mind[binm[1],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x3.temp)}
            ie2[,,binm[1],1]=med0$BUGSoutput$sims.list$beta[,binm1[1]]*(exp(temp.x3)/(1+exp(temp.x3))-exp(temp.x1)/(1+exp(temp.x1)))/deltap[l]
            aie2[,binm[1],1] <-apply(ie2[,,binm[1],1],1,mean,na.rm=T)}
          else
            for (k in 1:p2){
              if(nmc==1 & p2==1){
                temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,k]+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x1.temp)
                temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,k]+med0$BUGSoutput$sims.list$alpha1.b[,k,]%*%t(x3.temp)
              }
              else{
                temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,k,mind[binm[k],]]%*%t(mcov[,mind[binm[k],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x1.temp)
                temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,k,mind[binm[k],]]%*%t(mcov[,mind[binm[k],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x3.temp)
              }
              ie2[,,binm[k],l]=med0$BUGSoutput$sims.list$beta[,binm1[k]]*(exp(temp.x3)/(1+exp(temp.x3))-exp(temp.x1)/(1+exp(temp.x1)))/deltap[l]
              aie2[,binm[k],l] <- apply(ie2[,,binm[k],l],1,mean,na.rm=T)
            }}
        else if(y.type==2){
          if(p2==1 & c1==1){
            temp.M3=M1
            temp.M1=M1
            temp.M3[,binm1[1]]=1
            temp.M1[,binm1[1]]=0
            temp.mu1=med0$BUGSoutput$sims.list$c%*%t(x)+med0$BUGSoutput$sims.list$beta%*%t(temp.M1)+matrix(med0$BUGSoutput$sims.list$beta0,N1,nrow(x))
            temp.mu3=med0$BUGSoutput$sims.list$c%*%t(x)+med0$BUGSoutput$sims.list$beta%*%t(temp.M3)+matrix(med0$BUGSoutput$sims.list$beta0,N1,nrow(x))
            if(!is.null(cova))
            {temp.mu1=temp.mu1+med0$BUGSoutput$sims.list$eta%*%t(cova)
            temp.mu3=temp.mu3+med0$BUGSoutput$sims.list$eta%*%t(cova)}
            bpart=exp(temp.mu3)/(1+exp(temp.mu3))-exp(temp.mu1)/(1+exp(temp.mu1))
            if (nmc==1){
              temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,1]+matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x1.temp)
              temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,1]+matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x3.temp)
            }
            else
            {temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,1,mind[binm[1],]]%*%t(mcov[,mind[binm[1],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x1.temp)
            temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,1,mind[binm[1],]]%*%t(mcov[,mind[binm[1],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x3.temp)}
            ie2[,,binm[1],1]=bpart*(exp(temp.x3)/(1+exp(temp.x3))-exp(temp.x1)/(1+exp(temp.x1)))/deltap[l]
            aie2[,binm[1],1] <-apply(ie2[,,binm[1],1],1,mean,na.rm=T)}
          else
            for (k in 1:p2){
              temp.M3=M1
              temp.M1=M1
              temp.M3[,binm1[k]]=1
              temp.M1[,binm1[k]]=0
              temp.mu1=med0$BUGSoutput$sims.list$c%*%t(x)+med0$BUGSoutput$sims.list$beta%*%t(temp.M1)+matrix(med0$BUGSoutput$sims.list$beta0,N1,nrow(x))
              temp.mu3=med0$BUGSoutput$sims.list$c%*%t(x)+med0$BUGSoutput$sims.list$beta%*%t(temp.M3)+matrix(med0$BUGSoutput$sims.list$beta0,N1,nrow(x))
              if(!is.null(cova))
              {temp.mu1=temp.mu1+med0$BUGSoutput$sims.list$eta%*%t(cova)
              temp.mu3=temp.mu3+med0$BUGSoutput$sims.list$eta%*%t(cova)}
              bpart=exp(temp.mu3)/(1+exp(temp.mu3))-exp(temp.mu1)/(1+exp(temp.mu1))

              if(nmc==1 & p2==1){
                temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,k]+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x1.temp)
                temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,k]+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x3.temp)
              }
              else{
                temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,k,mind[binm[k],]]%*%t(mcov[,mind[binm[k],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x1.temp)
                temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,k,mind[binm[k],]]%*%t(mcov[,mind[binm[k],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x3.temp)
              }
              ie2[,,binm[k],l]=bpart*(exp(temp.x3)/(1+exp(temp.x3))-exp(temp.x1)/(1+exp(temp.x1)))/deltap[l]
              aie2[,binm[k],l] <- apply(ie2[,,binm[k],l],1,mean,na.rm=T)
            }}
        else if(y.type==4){
          if(p2==1 & c1==1){
            temp.M3=M1
            temp.M1=M1
            temp.M3[,binm1[1]]=1
            temp.M1[,binm1[1]]=0
            temp.mu1=med0$BUGSoutput$sims.list$c%*%t(x)+med0$BUGSoutput$sims.list$beta%*%t(temp.M1)
            temp.mu3=med0$BUGSoutput$sims.list$c%*%t(x)+med0$BUGSoutput$sims.list$beta%*%t(temp.M3)
            if(!is.null(cova))
            {temp.mu1=temp.mu1+med0$BUGSoutput$sims.list$eta%*%t(cova)
            temp.mu3=temp.mu3+med0$BUGSoutput$sims.list$eta%*%t(cova)}
            vec=cbind(med0$BUGSoutput$sims.list$r,
                      as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(temp.mu3),2,expp, 1/med0$BUGSoutput$sims.list$r)))
            tmean3=t(apply(vec,1,weib.trunc.mean,right=tmax))
            #if(multi)
            #  tmean3=log(tmean3)
            vec=cbind(med0$BUGSoutput$sims.list$r,
                      as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(temp.mu1),2,expp, 1/med0$BUGSoutput$sims.list$r)))
            tmean1=t(apply(vec,1,weib.trunc.mean,right=tmax))
            #if(multi)
            #  tmean1=log(tmean1)
            #tmean3=ifelse(tmean3==Inf,tmax-runif(1,0,0.1),tmean3)
            #tmean1=ifelse(tmean1==Inf,tmax-runif(1,0,0.1),tmean1)
            bpart=tmean3-tmean1
            if (nmc==1){
              temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,1]+matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x1.temp)
              temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,1]+matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x3.temp)
            }
            else
            {temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,1,mind[binm[1],]]%*%t(mcov[,mind[binm[1],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x1.temp)
            temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,1,mind[binm[1],]]%*%t(mcov[,mind[binm[1],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x3.temp)}
            ie2[,,binm[1],1]=bpart*(exp(temp.x3)/(1+exp(temp.x3))-exp(temp.x1)/(1+exp(temp.x1)))/deltap[l]
            aie2[,binm[1],1] <-apply(ie2[,,binm[1],1],1,mean,na.rm=T)}
          else
            for (k in 1:p2){
              temp.M3=M1
              temp.M1=M1
              temp.M3[,binm1[k]]=1
              temp.M1[,binm1[k]]=0
              temp.mu1=med0$BUGSoutput$sims.list$c%*%t(x)+med0$BUGSoutput$sims.list$beta%*%t(temp.M1)
              temp.mu3=med0$BUGSoutput$sims.list$c%*%t(x)+med0$BUGSoutput$sims.list$beta%*%t(temp.M3)
              if(!is.null(cova))
              {temp.mu1=temp.mu1+med0$BUGSoutput$sims.list$eta%*%t(cova)
              temp.mu3=temp.mu3+med0$BUGSoutput$sims.list$eta%*%t(cova)}
              vec=cbind(med0$BUGSoutput$sims.list$r,
                        as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(temp.mu3),2,expp, 1/med0$BUGSoutput$sims.list$r)))
              tmean3=t(apply(vec,1,weib.trunc.mean,right=tmax))
              #if(multi)
              #  tmean3=log(tmean3)
              vec=cbind(med0$BUGSoutput$sims.list$r,
                        as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(temp.mu1),2,expp, 1/med0$BUGSoutput$sims.list$r)))
              tmean1=t(apply(vec,1,weib.trunc.mean,right=tmax))
              #if(multi)
              #  tmean1=log(tmean1)
              #tmean3=ifelse(tmean3==Inf,tmax-runif(1,0,0.1),tmean3)
              #tmean1=ifelse(tmean1==Inf,tmax-runif(1,0,0.1),tmean1)
              bpart=tmean3-tmean1

              if(nmc==1 & p2==1){
                temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,k]+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x1.temp)
                temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,k]+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x3.temp)
              }
              else{
                temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,k,mind[binm[k],]]%*%t(mcov[,mind[binm[k],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x1.temp)
                temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,k,mind[binm[k],]]%*%t(mcov[,mind[binm[k],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x3.temp)
              }
              ie2[,,binm[k],l]=bpart*(exp(temp.x3)/(1+exp(temp.x3))-exp(temp.x1)/(1+exp(temp.x1)))/deltap[l]
              aie2[,binm[k],l] <- apply(ie2[,,binm[k],l],1,mean,na.rm=T)
            }}
      }

      if(p3>0){
        if(y.type==1)
          for (j in 1:p3){
            mu_Mc1<-array(0,c(N1,N,cat2[j]-1))
            mu_Mc0<-array(0,c(N1,N,cat2[j]-1))
            for (k in 2:cat2[j]){
              mu_Mc1[,,k-1] <- exp(as.matrix(med0$BUGSoutput$sims.list$alpha0.c[,j,k-1,mind[catm[j],]])%*%t(mcov[,mind[catm[j],]])+med0$BUGSoutput$sims.list$alpha1.c[,j,k-1,]%*%t(x3.temp))
              mu_Mc0[,,k-1] <- exp(as.matrix(med0$BUGSoutput$sims.list$alpha0.c[,j,k-1,mind[catm[j],]])%*%t(mcov[,mind[catm[j],]])+med0$BUGSoutput$sims.list$alpha1.c[,j,k-1,]%*%t(x1.temp))
            }
            sum_Mc1 <-apply(mu_Mc1,c(1,2),sum)+1
            sum_Mc0 <-apply(mu_Mc0,c(1,2),sum)+1
            for (k in 2:cat2[j])
              ie2[,,catm[j],l]=ie2[,,catm[j],l]+(mu_Mc1[,,k-1]/sum_Mc1-mu_Mc0[,,k-1]/sum_Mc0)*med0$BUGSoutput$sims.list$beta[,catm1[j,1]+k-2]/deltap[l]
            aie2[,catm[j],l]<-apply(ie2[,,catm[j],l],1,mean,na.rm=T)
          }
        else if(y.type==2)
          for (j in 1:p3){
            bpart=array(0,c(N1,N,cat2[j]-1))
            M1.temp=M1
            M1.temp[,catm1[j,1]:catm1[j,2]]=0
            for (k in catm1[j,1]:catm1[j,2]){
              temp.M3=M1.temp
              temp.M1=M1.temp
              temp.M3[,k]=1
              temp.mu1=med0$BUGSoutput$sims.list$c%*%t(x)+med0$BUGSoutput$sims.list$beta%*%t(temp.M1)+matrix(med0$BUGSoutput$sims.list$beta0,N1,nrow(x))
              temp.mu3=med0$BUGSoutput$sims.list$c%*%t(x)+med0$BUGSoutput$sims.list$beta%*%t(temp.M3)+matrix(med0$BUGSoutput$sims.list$beta0,N1,nrow(x))
              if(!is.null(cova))
              {temp.mu1=temp.mu1+med0$BUGSoutput$sims.list$eta%*%t(cova)
              temp.mu3=temp.mu3+med0$BUGSoutput$sims.list$eta%*%t(cova)}
              bpart[,,k-catm1[j,1]+1]=exp(temp.mu3)/(1+exp(temp.mu3))-exp(temp.mu1)/(1+exp(temp.mu1))
            }

            mu_Mc1<-array(0,c(N1,N,cat2[j]-1))
            mu_Mc0<-array(0,c(N1,N,cat2[j]-1))
            for (k in 2:cat2[j]){
              mu_Mc1[,,k-1] <- exp(as.matrix(med0$BUGSoutput$sims.list$alpha0.c[,j,k-1,mind[catm[j],]])%*%t(mcov[,mind[catm[j],]])+med0$BUGSoutput$sims.list$alpha1.c[,j,k-1,]%*%t(x3.temp))
              mu_Mc0[,,k-1] <- exp(as.matrix(med0$BUGSoutput$sims.list$alpha0.c[,j,k-1,mind[catm[j],]])%*%t(mcov[,mind[catm[j],]])+med0$BUGSoutput$sims.list$alpha1.c[,j,k-1,]%*%t(x1.temp))
            }
            sum_Mc1 <-apply(mu_Mc1,c(1,2),sum)+1
            sum_Mc0 <-apply(mu_Mc0,c(1,2),sum)+1
            for (k in 2:cat2[j])
              ie2[,,catm[j],l]=ie2[,,catm[j],l]+(mu_Mc1[,,k-1]/sum_Mc1-mu_Mc0[,,k-1]/sum_Mc0)*bpart[,,k-1]

            # for (k in 1:N)
            #    ie2[,k,catm[j],l]=diag((diag(1/sum_Mc1[,k])%*%mu_Mc1[,k,]-diag(1/sum_Mc0[,k])%*%mu_Mc0[,k,])%*%t(bpart[,k,]))
            aie2[,catm[j],l]<-apply(ie2[,,catm[j],l],1,mean,na.rm=T)
          }
        else if(y.type==4)
          for (j in 1:p3){
            bpart=array(0,c(N1,N,cat2[j]-1))
            M1.temp=M1
            M1.temp[,catm1[j,1]:catm1[j,2]]=0
            for (k in catm1[j,1]:catm1[j,2]){
              temp.M3=M1.temp
              temp.M1=M1.temp
              temp.M3[,k]=1
              temp.mu1=med0$BUGSoutput$sims.list$c%*%t(x)+med0$BUGSoutput$sims.list$beta%*%t(temp.M1)
              temp.mu3=med0$BUGSoutput$sims.list$c%*%t(x)+med0$BUGSoutput$sims.list$beta%*%t(temp.M3)
              if(!is.null(cova))
              {temp.mu1=temp.mu1+med0$BUGSoutput$sims.list$eta%*%t(cova)
              temp.mu3=temp.mu3+med0$BUGSoutput$sims.list$eta%*%t(cova)}
              vec=cbind(med0$BUGSoutput$sims.list$r,
                        as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(temp.mu3),2,expp, 1/med0$BUGSoutput$sims.list$r)))
              tmean3=t(apply(vec,1,weib.trunc.mean,right=tmax))
              #if(multi)
              #  tmean3=log(tmean3)
              vec=cbind(med0$BUGSoutput$sims.list$r,
                        as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(temp.mu1),2,expp, 1/med0$BUGSoutput$sims.list$r)))
              tmean1=t(apply(vec,1,weib.trunc.mean,right=tmax))
              #if(multi)
              #  tmean1=log(tmean1)
              #tmean3=ifelse(tmean3==Inf,tmax-runif(1,0,0.1),tmean3)
              #tmean1=ifelse(tmean1==Inf,tmax-runif(1,0,0.1),tmean1)
              bpart[,,k-catm1[j,1]+1]=tmean3-tmean1
              #          bpart[,,k-catm1[j,1]+1]=as.vector(gamma(1+1/med0$BUGSoutput$sims.list$r)/med0$BUGSoutput$sims.list$lambda^(1/med0$BUGSoutput$sims.list$r))*
              #            (1/apply(exp(temp.mu3),2,expp, 1/med0$BUGSoutput$sims.list$r)-1/apply(exp(temp.mu1),2,expp,1/med0$BUGSoutput$sims.list$r))
            }

            mu_Mc1<-array(0,c(N1,N,cat2[j]-1))
            mu_Mc0<-array(0,c(N1,N,cat2[j]-1))
            for (k in 2:cat2[j]){
              mu_Mc1[,,k-1] <- exp(as.matrix(med0$BUGSoutput$sims.list$alpha0.c[,j,k-1,mind[catm[j],]])%*%t(mcov[,mind[catm[j],]])+med0$BUGSoutput$sims.list$alpha1.c[,j,k-1,]%*%t(x3.temp))
              mu_Mc0[,,k-1] <- exp(as.matrix(med0$BUGSoutput$sims.list$alpha0.c[,j,k-1,mind[catm[j],]])%*%t(mcov[,mind[catm[j],]])+med0$BUGSoutput$sims.list$alpha1.c[,j,k-1,]%*%t(x1.temp))
            }
            sum_Mc1 <-apply(mu_Mc1,c(1,2),sum)+1
            sum_Mc0 <-apply(mu_Mc0,c(1,2),sum)+1
            for (k in 2:cat2[j])
              ie2[,,catm[j],l]=ie2[,,catm[j],l]+(mu_Mc1[,,k-1]/sum_Mc1-mu_Mc0[,,k-1]/sum_Mc0)*bpart[,,k-1]
            aie2[,catm[j],l]<-apply(ie2[,,catm[j],l],1,mean,na.rm=T)
          }
      }

      if(y.type==2)
      {temp.x1=x
      temp.x3=x
      if(l%in%data0$contpred1.0)
      {tt1=match(l,data0$contpred1.0)
      temp.x3[,data0$contpred2[tt1,1]:data0$contpred2[tt1,2]]=data0$pred3[,data0$contpred3[tt1,1]:data0$contpred3[tt1,2]]}
      else if(l%in%data0$binpred1.0)
      {tt1=match(l,data0$binpred1)
      temp.x1[,data0$binpred2[tt1]]=0
      temp.x3[,data0$binpred2[tt1]]=1
      deltap[l]=1}
      else
      {for (i in 1:nrow(data0$catpred1.0))
        if(l%in%(data0$catpred1.0[i,1]:data0$catpred1.0[i,2]))
          tt1=i
      temp.x1[,data0$catpred2[tt1,1]:data0$catpred2[tt1,2]]=0
      temp.x3[,data0$catpred2[tt1,1]:data0$catpred2[tt1,2]]=0
      temp.x3[,l]=1
      deltap[l]=1}
      temp.mu1=med0$BUGSoutput$sims.list$c%*%t(temp.x1)+med0$BUGSoutput$sims.list$beta%*%t(M1)+matrix(med0$BUGSoutput$sims.list$beta0,N1,nrow(x))
      temp.mu3=med0$BUGSoutput$sims.list$c%*%t(temp.x3)+med0$BUGSoutput$sims.list$beta%*%t(M1)+matrix(med0$BUGSoutput$sims.list$beta0,N1,nrow(x))
      if(!is.null(cova))
      {temp.mu1=temp.mu1+med0$BUGSoutput$sims.list$eta%*%t(cova)
      temp.mu3=temp.mu3+med0$BUGSoutput$sims.list$eta%*%t(cova)}
      de2.1=(exp(temp.mu3)/(1+exp(temp.mu3))-exp(temp.mu1)/(1+exp(temp.mu1)))/deltap[l]
      de2=cbind(de2,apply(de2.1,1,mean,na.rm=T))}
      else if(y.type==4)
      {temp.x1=x
      temp.x3=x
      if(l%in%data0$contpred1.0)
      {tt1=match(l,data0$contpred1.0)
      temp.x3[,data0$contpred2[tt1,1]:data0$contpred2[tt1,2]]=data0$pred3[,data0$contpred3[tt1,1]:data0$contpred3[tt1,2]]}
      else if(l%in%data0$binpred1.0)
      {tt1=match(l,data0$binpred1)
      temp.x1[,data0$binpred2[tt1]]=0
      temp.x3[,data0$binpred2[tt1]]=1
      deltap[l]=1}
      else
      {for (i in 1:nrow(data0$catpred1.0))
        if(l%in%(data0$catpred1.0[i,1]:data0$catpred1.0[i,2]))
          tt1=i
      temp.x1[,data0$catpred2[tt1,1]:data0$catpred2[tt1,2]]=0
      temp.x3[,data0$catpred2[tt1,1]:data0$catpred2[tt1,2]]=0
      temp.x3[,l]=1
      deltap[l]=1}
      temp.mu1=med0$BUGSoutput$sims.list$c%*%t(temp.x1)+med0$BUGSoutput$sims.list$beta%*%t(M1)
      temp.mu3=med0$BUGSoutput$sims.list$c%*%t(temp.x3)+med0$BUGSoutput$sims.list$beta%*%t(M1)
      if(!is.null(cova))
      {temp.mu1=temp.mu1+med0$BUGSoutput$sims.list$eta%*%t(cova)
      temp.mu3=temp.mu3+med0$BUGSoutput$sims.list$eta%*%t(cova)}
      vec=cbind(med0$BUGSoutput$sims.list$r,
                as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(temp.mu3),2,expp, 1/med0$BUGSoutput$sims.list$r)))
      tmean3=t(apply(vec,1,weib.trunc.mean,right=tmax))
      #if(multi)
      #  tmean3=log(tmean3)
      vec=cbind(med0$BUGSoutput$sims.list$r,
                as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(temp.mu1),2,expp, 1/med0$BUGSoutput$sims.list$r)))
      tmean1=t(apply(vec,1,weib.trunc.mean,right=tmax))
      #if(multi)
      #  tmean1=log(tmean1)
      de2.1=(tmean3-tmean1)/deltap[l]
      de2=cbind(de2,apply(de2.1,1,mean,na.rm=T))}

      #method 3:parametric
      #3.1. get M1(x) and M1(x+dx) for y
      if(p1>0){
        if(c1==1 & p1+p2+p3==1 & contm1[1,1]==contm1[1,2]){
          if(nmc==1)
          {mu_M2[,contm1[1,1],] <- med0$BUGSoutput$sims.list$alpha0[,contm1[1,1]]+as.matrix(med0$BUGSoutput$sims.list$alpha1[,contm1[1,1]])%*%t(x1.temp)
          mu_M3[,contm1[1,1],] <- med0$BUGSoutput$sims.list$alpha0[,contm1[1,1]]+as.matrix(med0$BUGSoutput$sims.list$alpha1[,contm1[1,1]])%*%t(x3.temp)}
          else
          {mu_M2[,contm1[1,1],] <- med0$BUGSoutput$sims.list$alpha0[,contm1[1,1],mind[contm[1],]]%*%t(mcov[,mind[contm[1],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1[,contm1[1,1]])%*%t(x1.temp)
          mu_M3[,contm1[1,1],] <- med0$BUGSoutput$sims.list$alpha0[,contm1[1,1],mind[contm[1],]]%*%t(mcov[,mind[contm[1],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1[,contm1[1,1]])%*%t(x3.temp)}
        }
        else
          for (j in 1:p1){
            for (k in contm1[j,1]:contm1[j,2]){
              if(nmc==1 & p1+p2+p3==1)
              {mu_M2[,k,] <- med0$BUGSoutput$sims.list$alpha0[,k]+med0$BUGSoutput$sims.list$alpha1[,k,l]%*%t(x1.temp)
              mu_M3[,k,] <- med0$BUGSoutput$sims.list$alpha0[,k]+med0$BUGSoutput$sims.list$alpha1[,k,l]%*%t(x3.temp)}
              else
              {mu_M2[,k,] <- med0$BUGSoutput$sims.list$alpha0[,k,mind[contm[j],]]%*%t(mcov[,mind[contm[j],]])+med0$BUGSoutput$sims.list$alpha1[,k,]%*%t(x1.temp)
              mu_M3[,k,] <- med0$BUGSoutput$sims.list$alpha0[,k,mind[contm[j],]]%*%t(mcov[,mind[contm[j],]])+med0$BUGSoutput$sims.list$alpha1[,k,]%*%t(x3.temp)}
            }}}

      if(p2>0){
        if(p2==1 & c1==1)
        {if(nmc==1)
        {temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,1]+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x1.temp)
        temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,1]+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x3.temp)}
          else
          {temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,1,mind[binm[1],]]%*%t(mcov[,mind[binm[1],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x1.temp)
          temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,1,mind[binm[1],]]%*%t(mcov[,mind[binm[1],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,1])%*%t(x3.temp)}
          mu_M2[,binm[1],] <- exp(temp.x1)/(1+exp(temp.x1))
          mu_M3[,binm[1],] <- exp(temp.x3)/(1+exp(temp.x3))
        }
        else{
          for (k in 1:p2){
            if(nmc==1 & p2==1)
            {temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,k]+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x1.temp)
            temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,k]+med0$BUGSoutput$sims.list$alpha1.b[,k,]%*%t(x3.temp)}
            else
            {temp.x1=med0$BUGSoutput$sims.list$alpha0.b[,k,mind[binm[k],]]%*%t(mcov[,mind[binm[k],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x1.temp)
            temp.x3=med0$BUGSoutput$sims.list$alpha0.b[,k,mind[binm[k],]]%*%t(mcov[,mind[binm[k],]])+as.matrix(med0$BUGSoutput$sims.list$alpha1.b[,k,])%*%t(x3.temp)}
            mu_M2[,binm1[k],] <- exp(temp.x1)/(1+exp(temp.x1))
            mu_M3[,binm1[k],] <- exp(temp.x3)/(1+exp(temp.x3))}
        }}

      if(p3>0){
        for (j in 1:p3){
          mu_Mc1<-array(0,c(N1,N,cat2[j]-1))
          mu_Mc0<-array(0,c(N1,N,cat2[j]-1))
          for (k in 2:cat2[j]){
            mu_Mc1[,,k-1] <- exp(as.matrix(med0$BUGSoutput$sims.list$alpha0.c[,j,k-1,mind[catm[j],]])%*%t(mcov[,mind[catm[j],]])+med0$BUGSoutput$sims.list$alpha1.c[,j,k-1,]%*%t(x3.temp))
            mu_Mc0[,,k-1] <- exp(as.matrix(med0$BUGSoutput$sims.list$alpha0.c[,j,k-1,mind[catm[j],]])%*%t(mcov[,mind[catm[j],]])+med0$BUGSoutput$sims.list$alpha1.c[,j,k-1,]%*%t(x1.temp))
          }
          sum_Mc1 <-apply(mu_Mc1,c(1,2),sum)+1
          sum_Mc0 <-apply(mu_Mc0,c(1,2),sum)+1
          for (k in 2:cat2[j])
          {mu_M2[,catm1[j,1]+k-2,]=mu_Mc0[,,k-1]/sum_Mc0
          mu_M3[,catm1[j,1]+k-2,]=mu_Mc1[,,k-1]/sum_Mc1}
        }}

      #3.2. get x and dx for y
      x1.temp1=x
      x3.temp1=x
      if(l%in%as.vector(data0$contpred1.0)){ #need the continuous predictor in its original format
        i=match(l,data0$contpred1.0)
        x3.temp1[,data0$contpred2[i,1]:data0$contpred2[i,2]]=data0$pred3[,data0$contpred3[i,1]:data0$contpred3[i,2]]
      }
      else if(l%in%data0$binpred1.0)
      {i=match(l,data0$binpred1.0)
      x1.temp1[,data0$binpred2[i]]=0
      x3.temp1[,data0$binpred2[i]]=1
      deltap[l]=1}
      else{ #categorical predictor
        for (i in 1:nrow(data0$catpred1.0))
          if(l%in%(data0$catpred1.0[i,1]:data0$catpred1.0[i,2]))
          {x1.temp1[,data0$catpred2[i,1]:data0$catpred2[i,2]]=0
          x3.temp1[,data0$catpred2[i,1]:data0$catpred2[i,2]]=0
          di=match(l,data0$catpred1.0[i,1]:data0$catpred1.0[i,2])
          x3.temp1[,data0$catpred2[i,1]+di-1]=1}
        deltap[l]=1
      }

      #3.3. get the total effect
      temp1=array(0,c(nrow(M1)+1,ncol(M1),N1))
      temp1[2:(nrow(M1)+1),,]=mu_M2
      temp1[1,,]=t(med0$BUGSoutput$sims.list$beta)
      temp2=temp1
      temp2[2:(nrow(M1)+1),,]=mu_M3
      if(is.null(cova))
      {mu_y0<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x1.temp1) + t(apply(temp1,3,matrix.prod))
      mu_y1<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x3.temp1) + t(apply(temp2,3,matrix.prod))
      }
      else
      {mu_y0<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x1.temp1) + med0$BUGSoutput$sims.list$eta%*%t(as.matrix(cova))+t(apply(temp1,3,matrix.prod))
      mu_y1<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x3.temp1) + med0$BUGSoutput$sims.list$eta%*%t(as.matrix(cova))+t(apply(temp2,3,matrix.prod))
      } #get the linear part

      if(y.type==2)
        te3=(exp(mu_y1)/(1+exp(mu_y1))-exp(mu_y0)/(1+exp(mu_y0)))/deltap[l]
      else if(y.type==1)
        te3<- (mu_y1-mu_y0)/deltap[l]
      else if (y.type==4)
      {mu_y0<- mu_y0-matrix(med0$BUGSoutput$sims.list$beta0,N1,N)
      mu_y1<- mu_y1-matrix(med0$BUGSoutput$sims.list$beta0,N1,N)
      vec=cbind(med0$BUGSoutput$sims.list$r,
                as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(mu_y1),2,expp, 1/med0$BUGSoutput$sims.list$r)))
      tmean3=t(apply(vec,1,weib.trunc.mean,right=tmax))
      vec=cbind(med0$BUGSoutput$sims.list$r,
                as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(mu_y0),2,expp, 1/med0$BUGSoutput$sims.list$r)))
      tmean1=t(apply(vec,1,weib.trunc.mean,right=tmax))
      if(multi)
      {tmean1=log(tmean1)
      tmean3=log(tmean3)}

      te3=(tmean3-tmean1)/deltap[l]
      }

      ate3[,l]=apply(te3,1,mean,na.rm=T)
      if(y.type==4){
        if(multi)
          omu3[,l]=apply(exp(tmean1),1,mean,na.rm=T)
        else
          omu3[,l]=apply(tmean1,1,mean,na.rm=T)}

      #3.4. calculate the ie
      j1=sample(1:N,size=N*N1,replace=T)
      j2=sample(1:N,size=N*N1,replace=T)

      #3.4.1. continuous mediators
      if(p1>0){
        for (j in 1:p1){
          temp1.1=temp1
          temp1.2=temp2
          for (i in contm1[j,1]:contm1[j,2])
          {temp1.1[2:(nrow(M1)+1),i,]=matrix(M1[j1,i],N,N1)
          temp1.2[2:(nrow(M1)+1),i,]=matrix(M1[j1,i],N,N1)} #j2/j1
          if(is.null(cova))
          {mu_y0.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x1.temp1) +t(apply(temp1.1,3,matrix.prod))
          mu_y1.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x3.temp1) +t(apply(temp1.2,3,matrix.prod)) }
          else{
            mu_y0.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x1.temp1) + med0$BUGSoutput$sims.list$eta%*%t(cova)+t(apply(temp1.1,3,matrix.prod))
            mu_y1.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x3.temp1) + med0$BUGSoutput$sims.list$eta%*%t(cova)+t(apply(temp1.2,3,matrix.prod))}
          if(y.type==1)
            ie3[,,contm[j],l]=te3-(mu_y1.2-mu_y0.2)/deltap[l]
          else if(y.type==2)
            ie3[,,contm[j],l]=te3-(exp(mu_y1.2)/(1+exp(mu_y1.2))-exp(mu_y0.2)/(1+exp(mu_y0.2)))/deltap[l]
          else if (y.type==4)
          {mu_y0.2<- mu_y0.2-matrix(med0$BUGSoutput$sims.list$beta0,N1,N)
          mu_y1.2<- mu_y1.2-matrix(med0$BUGSoutput$sims.list$beta0,N1,N)
          vec=cbind(med0$BUGSoutput$sims.list$r,
                    as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(mu_y1.2),2,expp, 1/med0$BUGSoutput$sims.list$r)))
          tmean3=t(apply(vec,1,weib.trunc.mean,right=tmax))
          if(multi)
            tmean3=log(tmean3)
          vec=cbind(med0$BUGSoutput$sims.list$r,
                    as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(mu_y0.2),2,expp, 1/med0$BUGSoutput$sims.list$r)))
          tmean1=t(apply(vec,1,weib.trunc.mean,right=tmax))
          if(multi)
            tmean1=log(tmean1)
          ie3[,,contm[j],l]<-te3-(tmean3-tmean1)/deltap[l]
          }
        }}

      #3.4.2. binary mediators
      if(p2>0){
        for (k in 1:p2){
          temp1.1=temp1
          temp1.2=temp2
          temp1.1[2:(nrow(M1)+1),binm1[k],]=matrix(M1[j1,binm1[k]],N,N1)
          temp1.2[2:(nrow(M1)+1),binm1[k],]=matrix(M1[j1,binm1[k]],N,N1) #j2/j1
          if (is.null(cova)){
            mu_y0.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x1.temp1) + t(apply(temp1.1,3,matrix.prod))
            mu_y1.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x3.temp1) + t(apply(temp1.2,3,matrix.prod))}
          else{
            mu_y0.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x1.temp1) + med0$BUGSoutput$sims.list$eta%*%t(cova)+t(apply(temp1.1,3,matrix.prod))
            mu_y1.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x3.temp1) + med0$BUGSoutput$sims.list$eta%*%t(cova)+t(apply(temp1.2,3,matrix.prod))}
          if(y.type==1)
            ie3[,,binm[k],l]=te3-(mu_y1.2-mu_y0.2)/deltap[l]
          else if(y.type==2)
            ie3[,,binm[k],l]=te3-(exp(mu_y1.2)/(1+exp(mu_y1.2))-exp(mu_y0.2)/(1+exp(mu_y0.2)))/deltap[l]
          else if (y.type==4)
          {mu_y0.2<- mu_y0.2-matrix(med0$BUGSoutput$sims.list$beta0,N1,N)
          mu_y1.2<- mu_y1.2-matrix(med0$BUGSoutput$sims.list$beta0,N1,N)
          vec=cbind(med0$BUGSoutput$sims.list$r,
                    as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(mu_y1.2),2,expp, 1/med0$BUGSoutput$sims.list$r)))
          tmean3=t(apply(vec,1,weib.trunc.mean,right=tmax))
          if(multi)
            tmean3=log(tmean3)
          vec=cbind(med0$BUGSoutput$sims.list$r,
                    as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(mu_y0.2),2,expp, 1/med0$BUGSoutput$sims.list$r)))
          tmean1=t(apply(vec,1,weib.trunc.mean,right=tmax))
          if(multi)
            tmean1=log(tmean1)
          ie3[,,binm[k],l]<-te3-(tmean3-tmean1)/deltap[l]
          }
        }}

      if(p3>0){
        for (j in 1:p3){
          temp1.1=temp1
          temp1.2=temp2
          for (i in catm1[j,1]:catm1[j,2])
          {temp1.1[2:(nrow(M1)+1),i,]=matrix(M1[j1,i],N,N1)
          temp1.2[2:(nrow(M1)+1),i,]=matrix(M1[j1,i],N,N1)} #j2/j1
          if(is.null(cova)){
            mu_y0.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x1.temp1) + t(apply(temp1.1,3,matrix.prod))
            mu_y1.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x3.temp1) + t(apply(temp1.2,3,matrix.prod))}
          else{
            mu_y0.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x1.temp1) + med0$BUGSoutput$sims.list$eta%*%t(cova)+t(apply(temp1.1,3,matrix.prod))
            mu_y1.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x3.temp1) + med0$BUGSoutput$sims.list$eta%*%t(cova)+t(apply(temp1.2,3,matrix.prod))}

          if(y.type==1)
            ie3[,,catm[j],l]=te3-(mu_y1.2-mu_y0.2)/deltap[l]
          else if(y.type==2)
            ie3[,,catm[j],l]=te3-(exp(mu_y1.2)/(1+exp(mu_y1.2))-exp(mu_y0.2)/(1+exp(mu_y0.2)))/deltap[l]
          else if (y.type==4)
          {mu_y0.2<- mu_y0.2-matrix(med0$BUGSoutput$sims.list$beta0,N1,N)
          mu_y1.2<- mu_y1.2-matrix(med0$BUGSoutput$sims.list$beta0,N1,N)
          vec=cbind(med0$BUGSoutput$sims.list$r,
                    as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(mu_y1.2),2,expp, 1/med0$BUGSoutput$sims.list$r)))
          tmean3=t(apply(vec,1,weib.trunc.mean,right=tmax))
          if(multi)
            tmean3=log(tmean3)
          vec=cbind(med0$BUGSoutput$sims.list$r,
                    as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(mu_y0.2),2,expp, 1/med0$BUGSoutput$sims.list$r)))
          tmean1=t(apply(vec,1,weib.trunc.mean,right=tmax))
          if(multi)
            tmean1=log(tmean1)
          ie3[,,catm[j],l]<-te3-(tmean3-tmean1)/deltap[l]
          }
        }}

      aie3[,,l]<-apply(array(ie3[,,,l],c(N1,N,p1+p2+p3)),c(1,3),mean,na.rm=T)

      #3.5. Calculate the de
      temp1.1=temp1
      temp1.2=temp2
      for (i in 1:ncol(M1))
      {temp1.1[2:(nrow(M1)+1),i,]=matrix(M1[j1,i],N,N1)
      temp1.2[2:(nrow(M1)+1),i,]=matrix(M1[j1,i],N,N1)} #j2/j1
      if(is.null(cova)){
        mu_y0.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x1.temp1) + t(apply(temp1.1,3,matrix.prod))
        mu_y1.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x3.temp1) + t(apply(temp1.2,3,matrix.prod))}
      else{
        mu_y0.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x1.temp1) + med0$BUGSoutput$sims.list$eta%*%t(cova)+t(apply(temp1.1,3,matrix.prod))
        mu_y1.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x3.temp1) + med0$BUGSoutput$sims.list$eta%*%t(cova)+t(apply(temp1.2,3,matrix.prod))}

      if(y.type==1)
        de3=(mu_y1.2-mu_y0.2)/deltap[l]
      else if(y.type==2)
        de3=(exp(mu_y1.2)/(1+exp(mu_y1.2))-exp(mu_y0.2)/(1+exp(mu_y0.2)))/deltap[l]
      else if (y.type==4)
      {mu_y0.2<- mu_y0.2-matrix(med0$BUGSoutput$sims.list$beta0,N1,N)
      mu_y1.2<- mu_y1.2-matrix(med0$BUGSoutput$sims.list$beta0,N1,N)
      vec=cbind(med0$BUGSoutput$sims.list$r,
                as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(mu_y1.2),2,expp, 1/med0$BUGSoutput$sims.list$r)))
      tmean3=t(apply(vec,1,weib.trunc.mean,right=tmax))
      if(multi)
        tmean3=log(tmean3)
      vec=cbind(med0$BUGSoutput$sims.list$r,
                as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(mu_y0.2),2,expp, 1/med0$BUGSoutput$sims.list$r)))
      tmean1=t(apply(vec,1,weib.trunc.mean,right=tmax))
      if(multi)
        tmean1=log(tmean1)
      de3=(tmean3-tmean1)/deltap[l]
      }

      ade3[,l]=apply(de3,1,mean,na.rm=T)

      #method3: semi-parametric for binary or categorical predictors
      if(!l%in%data0$contpred1.0){
        if(!is.null(data0$binpred1.0))
        {if (l%*%data0$binpred1.0)
        {M.0=data.frame(M1[x1[,l]==0,])
        y.0=y[x1[,l]==0]}
          else
          {for(i in 1:nrow(data0$catpred1.0))
            if(l%in%(data0$catpred1.0[i,1]:data0$catpred1.0[i,2]))
              tt1=i
          M.0=data.frame(M1[apply(x1[,data0$catpred1.0[tt1,1]:data0$catpred1.0[[tt1,2]]]==1,1,sum)==0,])
          y.0=y[apply(x1[,data0$catpred1.0[tt1,1]:data0$catpred1.0[[tt1,2]]]==1,1,sum)==0]}}
        else
        {for(i in 1:nrow(data0$catpred1.0))
          if(l%in%(data0$catpred1.0[i,1]:data0$catpred1.0[i,2]))
            tt1=i
        M.0=data.frame(M1[apply(x1[,data0$catpred1.0[tt1,1]:data0$catpred1.0[[tt1,2]]]==1,1,sum)==0,])
        y.0=y[apply(x1[,data0$catpred1.0[tt1,1]:data0$catpred1.0[[tt1,2]]]==1,1,sum)==0]}

        y.1=y[x1[,l]==1]
        M.1=data.frame(M1[x1[,l]==1,])

        j1=sample(1:N,size=N*N1,replace=T)
        j2=sample(1:N,size=N*N1,replace=T)
        n3=nrow(M.0)
        n4=nrow(M.1)
        j3=sample(1:n3,size=N*N1,replace = T)
        j4=sample(1:n4,size=N*N1,replace = T)

        for (i in 1:ncol(M1))
        {mu.M0[,i,]=matrix(M.0[j3,i],N,N1)
        mu.M1[,i,]=matrix(M.1[j4,i],N,N1)}

        #4.1. get the total effect
        mu_y0=matrix(y.0[j3],N1,N)
        mu_y1=matrix(y.1[j4],N1,N)

        temp1=array(0,c(nrow(M1)+1,ncol(M1),N1))
        temp1[2:(nrow(M1)+1),,]=mu.M0
        temp1[1,,]=t(med0$BUGSoutput$sims.list$beta)
        temp2=temp1
        temp2[2:(nrow(M1)+1),,]=mu.M1

        if(is.null(cova))
        {mu_y0.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x1.temp1) +t(apply(temp1,3,matrix.prod))
        mu_y1.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x3.temp1) +t(apply(temp2,3,matrix.prod)) }
        else{
          mu_y0.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x1.temp1) + med0$BUGSoutput$sims.list$eta%*%t(cova)+t(apply(temp1,3,matrix.prod))
          mu_y1.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x3.temp1) + med0$BUGSoutput$sims.list$eta%*%t(cova)+t(apply(temp2,3,matrix.prod))}
        if(y.type==1)
          te4<- (mu_y1.2-mu_y0.2)/deltap[l]
        else if(y.type==2)
          te4=(exp(mu_y1.2)/(1+exp(mu_y1.2))-exp(mu_y0.2)/(1+exp(mu_y0.2)))/deltap[l]
        else if (y.type==4)
        {mu_y0.2<- mu_y0.2-matrix(med0$BUGSoutput$sims.list$beta0,N1,N)
        mu_y1.2<- mu_y1.2-matrix(med0$BUGSoutput$sims.list$beta0,N1,N)
        vec=cbind(med0$BUGSoutput$sims.list$r,
                  as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(mu_y1.2),2,expp, 1/med0$BUGSoutput$sims.list$r)))
        tmean3=t(apply(vec,1,weib.trunc.mean,right=tmax))
        vec=cbind(med0$BUGSoutput$sims.list$r,
                  as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(mu_y0.2),2,expp, 1/med0$BUGSoutput$sims.list$r)))
        tmean1=t(apply(vec,1,weib.trunc.mean,right=tmax))
        if(multi)
        {tmean1=log(tmean1)
        tmean3=log(tmean3)}
        te4=(tmean3-tmean1)/deltap[l]
        }

        ate4[,l]=apply(te4,1,mean,na.rm=T)
        if(y.type==4){
          if(multi)
            omu4[,l]=apply(exp(tmean1),1,mean,na.rm=T)
          else
            omu4[,l]=apply(tmean1,1,mean,na.rm=T)}


        #4.2. Get the ies
        #ie for continuous mediators
        if(p1>0){
          for (j in 1:p1){
            temp1.1=temp1
            temp1.2=temp2
            for (i in contm1[j,1]:contm1[j,2])
            {temp1.1[2:(nrow(M1)+1),i,]=matrix(M1[j1,i],N,N1)
            temp1.2[2:(nrow(M1)+1),i,]=matrix(M1[j1,i],N,N1)} #j2/j1
            if(is.null(cova))
            {mu_y0.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x1.temp1) +t(apply(temp1.1,3,matrix.prod))
            mu_y1.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x3.temp1) +t(apply(temp1.2,3,matrix.prod)) }
            else{
              mu_y0.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x1.temp1) + med0$BUGSoutput$sims.list$eta%*%t(cova)+t(apply(temp1.1,3,matrix.prod))
              mu_y1.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x3.temp1) + med0$BUGSoutput$sims.list$eta%*%t(cova)+t(apply(temp1.2,3,matrix.prod))}

            if(y.type==1)
              ie4[,,contm[j],l]=te4-(mu_y1.2-mu_y0.2)/deltap[l]
            else if(y.type==2)
              ie4[,,contm[j],l]=te4-(exp(mu_y1.2)/(1+exp(mu_y1.2))-exp(mu_y0.2)/(1+exp(mu_y0.2)))/deltap[l]
            else if (y.type==4)
            {mu_y0.2<- mu_y0.2-matrix(med0$BUGSoutput$sims.list$beta0,N1,N)
            mu_y1.2<- mu_y1.2-matrix(med0$BUGSoutput$sims.list$beta0,N1,N)
            vec=cbind(med0$BUGSoutput$sims.list$r,
                      as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(mu_y1.2),2,expp, 1/med0$BUGSoutput$sims.list$r)))
            tmean3=t(apply(vec,1,weib.trunc.mean,right=tmax))
            if(multi)
              tmean3=log(tmean3)
            vec=cbind(med0$BUGSoutput$sims.list$r,
                      as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(mu_y0.2),2,expp, 1/med0$BUGSoutput$sims.list$r)))
            tmean1=t(apply(vec,1,weib.trunc.mean,right=tmax))
            if(multi)
              tmean1=log(tmean1)
            ie4[,,contm[j],l]=te4-(tmean3-tmean1)/deltap[l]
            }
          }}

        #ie for binary mediators
        if(p2>0){
          for (k in 1:p2){
            temp1.1=temp1
            temp1.2=temp2
            temp1.1[2:(nrow(M1)+1),binm1[k],]=matrix(M1[j1,binm1[k]],N,N1)
            temp1.2[2:(nrow(M1)+1),binm1[k],]=matrix(M1[j1,binm1[k]],N,N1) #j2/j1
            if (is.null(cova)){
              mu_y0.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x1.temp1) + t(apply(temp1.1,3,matrix.prod))
              mu_y1.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x3.temp1) + t(apply(temp1.2,3,matrix.prod))}
            else{
              mu_y0.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x1.temp1) + med0$BUGSoutput$sims.list$eta%*%t(cova)+t(apply(temp1.1,3,matrix.prod))
              mu_y1.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x3.temp1) + med0$BUGSoutput$sims.list$eta%*%t(cova)+t(apply(temp1.2,3,matrix.prod))}
            if(y.type==1)
              ie4[,,binm[k],l]=te4-(mu_y1.2-mu_y0.2)/deltap[l]
            else if(y.type==2)
              ie4[,,binm[k],l]=te4-(exp(mu_y1.2)/(1+exp(mu_y1.2))-exp(mu_y0.2)/(1+exp(mu_y0.2)))/deltap[l]
            else if (y.type==4)
            {mu_y0.2<- mu_y0.2-matrix(med0$BUGSoutput$sims.list$beta0,N1,N)
            mu_y1.2<- mu_y1.2-matrix(med0$BUGSoutput$sims.list$beta0,N1,N)
            vec=cbind(med0$BUGSoutput$sims.list$r,
                      as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(mu_y1.2),2,expp, 1/med0$BUGSoutput$sims.list$r)))
            tmean3=t(apply(vec,1,weib.trunc.mean,right=tmax))
            if(multi)
              tmean3=log(tmean3)
            vec=cbind(med0$BUGSoutput$sims.list$r,
                      as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(mu_y0.2),2,expp, 1/med0$BUGSoutput$sims.list$r)))
            tmean1=t(apply(vec,1,weib.trunc.mean,right=tmax))
            if(multi)
              tmean1=log(tmean1)
            ie4[,,binm[k],l]=te4-(tmean3-tmean1)/deltap[l]
            }
          }}

        #ie for categorical mediators
        if(p3>0){
          for (j in 1:p3){
            temp1.1=temp1
            temp1.2=temp2
            for (i in catm1[j,1]:catm1[j,2])
            {temp1.1[2:(nrow(M1)+1),i,]=matrix(M1[j1,i],N,N1)
            temp1.2[2:(nrow(M1)+1),i,]=matrix(M1[j1,i],N,N1)} #j2/j1
            if(is.null(cova)){
              mu_y0.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x1.temp1) + t(apply(temp1.1,3,matrix.prod))
              mu_y1.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x3.temp1) + t(apply(temp1.2,3,matrix.prod))}
            else{
              mu_y0.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x1.temp1) + med0$BUGSoutput$sims.list$eta%*%t(cova)+t(apply(temp1.1,3,matrix.prod))
              mu_y1.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x3.temp1) + med0$BUGSoutput$sims.list$eta%*%t(cova)+t(apply(temp1.2,3,matrix.prod))}

            if(y.type==1)
              ie4[,,catm[j],l]=te4-(mu_y1.2-mu_y0.2)/deltap[l]
            else if(y.type==2)
              ie4[,,catm[j],l]=te4-(exp(mu_y1.2)/(1+exp(mu_y1.2))-exp(mu_y0.2)/(1+exp(mu_y0.2)))/deltap[l]
            else if (y.type==4)
            {mu_y0.2<- mu_y0.2-matrix(med0$BUGSoutput$sims.list$beta0,N1,N)
            mu_y1.2<- mu_y1.2-matrix(med0$BUGSoutput$sims.list$beta0,N1,N)
            vec=cbind(med0$BUGSoutput$sims.list$r,
                      as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(mu_y1.2),2,expp, 1/med0$BUGSoutput$sims.list$r)))
            tmean3=t(apply(vec,1,weib.trunc.mean,right=tmax))
            if(multi)
              tmean3=log(tmean3)
            vec=cbind(med0$BUGSoutput$sims.list$r,
                      as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(mu_y0.2),2,expp, 1/med0$BUGSoutput$sims.list$r)))
            tmean1=t(apply(vec,1,weib.trunc.mean,right=tmax))
            if(multi)
              tmean1=log(tmean1)
            ie4[,,catm[j],l]=te4-(tmean3-tmean1)/deltap[l]
            }
          }}

        aie4[,,l]<-apply(array(ie4[,,,l],c(N1,N,p1+p2+p3)),c(1,3),mean,na.rm=T)

        #4.3. Calculate the de
        temp1.1=temp1
        temp1.2=temp2
        for (i in 1:ncol(M1))
        {temp1.1[2:(nrow(M1)+1),i,]=matrix(M1[j1,i],N,N1)
        temp1.2[2:(nrow(M1)+1),i,]=matrix(M1[j1,i],N,N1)} #j2/j1
        if(is.null(cova)){
          mu_y0.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x1.temp1) + t(apply(temp1.1,3,matrix.prod))
          mu_y1.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x3.temp1) + t(apply(temp1.2,3,matrix.prod))}
        else{
          mu_y0.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x1.temp1) + med0$BUGSoutput$sims.list$eta%*%t(cova)+t(apply(temp1.1,3,matrix.prod))
          mu_y1.2<- matrix(med0$BUGSoutput$sims.list$beta0,N1,N) + med0$BUGSoutput$sims.list$c%*%t(x3.temp1) + med0$BUGSoutput$sims.list$eta%*%t(cova)+t(apply(temp1.2,3,matrix.prod))}

        if(y.type==1)
          de4=(mu_y1.2-mu_y0.2)/deltap[l]
        else if(y.type==2)
          de4=(exp(mu_y1.2)/(1+exp(mu_y1.2))-exp(mu_y0.2)/(1+exp(mu_y0.2)))/deltap[l]
        else if (y.type==4)
        {mu_y0.2<- mu_y0.2-matrix(med0$BUGSoutput$sims.list$beta0,N1,N)
        mu_y1.2<- mu_y1.2-matrix(med0$BUGSoutput$sims.list$beta0,N1,N)
        vec=cbind(med0$BUGSoutput$sims.list$r,
                  as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(mu_y1.2),2,expp, 1/med0$BUGSoutput$sims.list$r)))
        tmean3=t(apply(vec,1,weib.trunc.mean,right=tmax))
        if(multi)
          tmean3=log(tmean3)
        vec=cbind(med0$BUGSoutput$sims.list$r,
                  as.vector(med0$BUGSoutput$sims.list$lambda^(-1/med0$BUGSoutput$sims.list$r))*(1/apply(exp(mu_y0.2),2,expp, 1/med0$BUGSoutput$sims.list$r)))
        tmean1=t(apply(vec,1,weib.trunc.mean,right=tmax))
        if(multi)
          tmean1=log(tmean1)
        de4=(tmean3-tmean1)/deltap[l]
        }
        ade4[,l]=apply(de4,1,mean,na.rm=T)
      }
    }

    if(y.type==1)
      de2=de1
    ate1=apply(aie1,c(1,3),sum)+de1
    ate2=apply(aie2,c(1,3),sum)+de2

    colnames(aie4)=colnames(M2)
    colnames(aie1)=colnames(M2)
    colnames(aie2)=colnames(M2)
    colnames(aie3)=colnames(M2)
  }
  #detach(med0$BUGSoutput$sims.list)

  result=list(aie1=aie1, ade1=de1, ate1=ate1, aie2=aie2, ade2=de2,  ate2=ate2,
              aie3=aie3, ade3=de3, ate3=ate3, aie4=aie4, ade4=ade4, ate4=ate4,
              sims.list=med0,data0=data0,omu3=omu3,omu4=omu4) #ie2=ie2, med0$BUGSoutput$sims.list
#  result=list(aie1=aie1, ade1=de1, ate1=ate1, aie2=aie2, ade2=de2,  ate2=ate2,
#              sims.list=med0,data0=data0)
  class(result)='hdbma'
  return(result)
}

summary.hdbma<-function(object, ..., plot= TRUE, RE=TRUE, quant=c(0.025, 0.25, 0.5, 0.75,0.975),digit=4,method=1)
{
  summary.med<-function(vec,qua=quant, digit=digit)
  {c(mean=mean(vec,na.rm=T),sd=sd(vec,na.rm=T),quantile(vec,qua,na.rm=T))
  }

  summary.med.re<-function(vec,vec1,qua=quant, digit=digit)
  {vec=vec/vec1
  c(mean=mean(vec,na.rm=T),sd=sd(vec,na.rm=T),quantile(vec,qua,na.rm=T))}

  if(object$data0$y_type!=3){
    result1<-array(0,c(7,2+ncol(object$aie1),ncol(object$ade1)))
    result1.re<-array(0,c(7,1+ncol(object$aie1),ncol(object$ade1)))
    result2<-result1
    result2.re<-result1.re
    result3<-result1
    result3.re<-result1.re
    result4<-result1
    result4.re<-result1.re
    for(j in 1:ncol(object$ade1)){
      result1[,,j]<-apply(cbind(TE=object$ate1[,j],DE=object$ade1[,j],object$aie1[,,j]),2,summary.med)
      result1.re[,,j]<-apply(cbind(DE=object$ade1[,j],object$aie1[,,j]),2,summary.med.re,object$ate1[,j])

      result2[,,j]<-apply(cbind(TE=object$ate2[,j],DE=object$ade2[,j],object$aie2[,,j]),2,summary.med)
      result2.re[,,j]<-apply(cbind(DE=object$ade2[,j],object$aie2[,,j]),2,summary.med.re,object$ate2[,j])

      result3[,,j]<-apply(cbind(TE=object$ate3[,j],DE=object$ade3[,j],object$aie3[,,j]),2,summary.med)
      result3.re[,,j]<-apply(cbind(DE=object$ade3[,j],object$aie3[,,j]),2,summary.med.re,object$ate3[,j])

      result4[,,j]<-apply(cbind(TE=object$ate4[,j],DE=object$ade4[,j],object$aie4[,,j]),2,summary.med)
      result4.re[,,j]<-apply(cbind(DE=object$ade4[,j],object$aie4[,,j]),2,summary.med.re,object$ate4[,j])}}
  else{
    result1<-array(0,c(7,2+ncol(object$aie1),ncol(object$ade1),dim(object$ade1)[3]))
    result1.re<-array(0,c(7,1+ncol(object$aie1),ncol(object$ade1),dim(object$ade1)[3]))
    result2<-result1
    result2.re<-result1.re
    result3<-result1
    result3.re<-result1.re
    result4<-result1
    result4.re<-result1.re
    for(j in 1:ncol(object$ade1)){
      for(q1 in 1:dim(object$ade1)[3]){
        result1[,,j,q1]<-apply(cbind(TE=object$ate1[,j,q1],DE=object$ade1[,j,q1],object$aie1[,,j,q1]),2,summary.med)
        result1.re[,,j,q1]<-apply(cbind(DE=object$ade1[,j,q1],object$aie1[,,j,q1]),2,summary.med.re,object$ate1[,j,q1])

        result2[,,j,q1]<-apply(cbind(TE=object$ate2[,j,q1],DE=object$ade2[,j,q1],object$aie2[,,j,q1]),2,summary.med)
        result2.re[,,j,q1]<-apply(cbind(DE=object$ade2[,j,q1],object$aie2[,,j,q1]),2,summary.med.re,object$ate2[,j,q1])

        result3[,,j,q1]<-apply(cbind(TE=object$ate3[,j,q1],DE=object$ade3[,j,q1],object$aie3[,,j,q1]),2,summary.med)
        result3.re[,,j,q1]<-apply(cbind(DE=object$ade3[,j,q1],object$aie3[,,j,q1]),2,summary.med.re,object$ate3[,j,q1])

        result4[,,j,q1]<-apply(cbind(TE=object$ate4[,j,q1],DE=object$ade4[,j,q1],object$aie4[,,j,q1]),2,summary.med)
        result4.re[,,j,q1]<-apply(cbind(DE=object$ade4[,j,q1],object$aie4[,,j,q1]),2,summary.med.re,object$ate4[,j,q1])}
    }
  }
  c.names=colnames(object$aie1)
  colnames(result1)=c("TE","DE",c.names)
  colnames(result2)=c("TE","DE",c.names)
  colnames(result3)=c("TE","DE",c.names)
  colnames(result4)=c("TE","DE",c.names)
  colnames(result1.re)=c("DE",c.names)
  colnames(result2.re)=c("DE",c.names)
  colnames(result3.re)=c("DE",c.names)
  colnames(result4.re)=c("DE",c.names)
  rownames(result1)=c("mean","sd",paste("q",quant,sep="_"))
  rownames(result2)=c("mean","sd",paste("q",quant,sep="_"))
  rownames(result3)=c("mean","sd",paste("q",quant,sep="_"))
  rownames(result4)=c("mean","sd",paste("q",quant,sep="_"))
  rownames(result1.re)=c("mean","sd",paste("q",quant,sep="_"))
  rownames(result2.re)=c("mean","sd",paste("q",quant,sep="_"))
  rownames(result3.re)=c("mean","sd",paste("q",quant,sep="_"))
  rownames(result4.re)=c("mean","sd",paste("q",quant,sep="_"))

  result=list(result1=result1,result1.re=result1.re,result2=result2,result2.re=result2.re,
              result3=result3,result3.re=result3.re,result4=result4,result4.re=result4.re,method=method,
              digit=digit,plot=plot,RE=RE,y.type=object$data0$y_type)
  class(result)="summary.hdbma"
  result
}

print.summary.hdbma<-function (x, ..., digit = x$digit, method=x$method, RE=x$RE)
{plot.sum<-function(obj1,main1="Estimated Effects")
{oldpar <- par(no.readonly = TRUE)
on.exit(par(oldpar))
re <- obj1[1,]
upper <- obj1[7,]
lower <- obj1[3,]
name1 <- colnames(obj1)
par(mfrow = c(1, 1), mar = c(1, 6, 1, 1), oma = c(3, 2, 2, 4))
bp <- barplot2(re, horiz = TRUE, main = main1,
               names.arg = name1, plot.ci = TRUE, ci.u = upper,
               ci.l = lower, cex.names = 0.9, beside = FALSE,
               cex.axis = 0.9, las = 1, xlim = range(c(upper,lower), na.rm = TRUE,finite=T),
               col = rainbow(length(re), start = 3/6, end = 4/6))}
if(x$y.type!=3){
  for (j in 1:dim(x$result1)[3])
  {cat('\n For Predictor', j, ':\n')
    if (method==1)
    {if(!RE)
    {cat('Estimated Effects for Method 1:\n')
      print(round(x$result1[,,j],digits = digit))
      if(x$plot)
        plot.sum(x$result1[,,j],main1=paste("Estimated Effects Using Method", method, "(predictor",j,")"))}
      else{
        cat('Estimated Relative Effects for Method 1:\n')
        print(round(x$result1.re[,,j],digits = digit))
        if(x$plot)
          plot.sum(x$result1.re[,,j],main1=paste("Estimated Relative Effects Using Method", method, "(predictor",j,")"))}}
    else if(method==2){
      if(!RE)
      {cat('Estimated Effects for Method 2:\n')
        print(round(x$result2[,,j],digits = digit))
        if(x$plot)
          plot.sum(x$result2[,,j],main1=paste("Estimated Effects Using Method", method, "(predictor",j,")"))}
      else{
        cat('Estimated Relative Effects for Method 2:\n')
        print(round(x$result2.re[,,j],digits = digit))
        if(x$plot)
          plot.sum(x$result2.re[,,j],main1=paste("Estimated Relative Effects Using Method", method, "(predictor",j,")"))}
    }
    else if(method==3){
      if(!RE)
      {cat('Estimated Effects for Method 3:\n')
        print(round(x$result3[,,j],digits = digit))
        if(x$plot)
          plot.sum(x$result3[,,j],main1=paste("Estimated Effects Using Method", method, "(predictor",j,")"))}
      else{
        cat('Estimated Relative Effects for Method 3:\n')
        print(round(x$result3.re[,,j],digits = digit))
        if(x$plot)
          plot.sum(x$result3.re[,,j],main1=paste("Estimated Relative Effects Using Method", method, "(predictor",j,")"))}
    }
    else if(method==4){
      if(!RE)
      {cat('Estimated Effects for Method 4:\n')
        print(round(x$result4[,,j],digits = digit))
        if(x$plot)
          plot.sum(x$result4[,,j],main1=paste("Estimated Effects Using Method", method, "(predictor",j,")"))}
      else{
        cat('Estimated Relative Effects for Method 4:\n')
        print(round(x$result4.re[,,j],digits = digit))
        if(x$plot)
          plot.sum(x$result4.re[,,j],main1=paste("Estimated Relative Effects Using Method", method, "(predictor",j,")"))}
    }
  }}
else{
  for (j in 1:dim(x$result1)[3])
    for (q1 in 1:dim(x$result1)[4])
    {cat('\n For Predictor', j, 'outcome',q1, ':\n')
      if (method==1)
      {if(!RE)
      {cat('Estimated Effects for Method 1:\n')
        print(round(x$result1[,,j,q1],digits = digit))
        if(x$plot)
          plot.sum(x$result1[,,j,q1],main1=paste("Estimated Effects Using Method", method, "(predictor",j,"outcome",q1,")"))}
        else{
          cat('Estimated Relative Effects for Method 1:\n')
          print(round(x$result1.re[,,j,q1],digits = digit))
          if(x$plot)
            plot.sum(x$result1.re[,,j,q1],main1=paste("Estimated Relative Effects Using Method", method, "(predictor",j,"outcome",q1,")"))}}
      else if(method==2){
        if(!RE)
        {cat('Estimated Effects for Method 2:\n')
          print(round(x$result2[,,j,q1],digits = digit))
          if(x$plot)
            plot.sum(x$result2[,,j,q1],main1=paste("Estimated Effects Using Method", method, "(predictor",j,"outcome",q1,")"))}
        else{
          cat('Estimated Relative Effects for Method 2:\n')
          print(round(x$result2.re[,,j,q1],digits = digit))
          if(x$plot)
            plot.sum(x$result2.re[,,j,q1],main1=paste("Estimated Relative Effects Using Method", method, "(predictor",j,"outcome",q1,")"))}
      }
      else if(method==3){
        if(!RE)
        {cat('Estimated Effects for Method 3:\n')
          print(round(x$result3[,,j,q1],digits = digit))
          if(x$plot)
            plot.sum(x$result3[,,j,q1],main1=paste("Estimated Effects Using Method", method, "(predictor",j,"outcome",q1,")"))}
        else{
          cat('Estimated Relative Effects for Method 3:\n')
          print(round(x$result3.re[,,j,q1],digits = digit))
          if(x$plot)
            plot.sum(x$result3.re[,,j,q1],main1=paste("Estimated Relative Effects Using Method", method, "(predictor",j,"outcome",q1,")"))}
      }
      else if(method==4){
        if(!RE)
        {cat('Estimated Effects for Method 4:\n')
          print(round(x$result4[,,j,q1],digits = digit))
          if(x$plot)
            plot.sum(x$result4[,,j,q1],main1=paste("Estimated Effects Using Method", method, "(predictor",j,"outcome",q1,")"))}
        else{
          cat('Estimated Relative Effects for Method 4:\n')
          print(round(x$result4.re[,,j,q1],digits = digit))
          if(x$plot)
            plot.sum(x$result4.re[,,j,q1],main1=paste("Estimated Relative Effects Using Method", method, "(predictor",j,"outcome",q1,")"))}
      }
    }
}
}
