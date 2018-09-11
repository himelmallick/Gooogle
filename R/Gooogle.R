#' PD approximation of a matrix.
#'
#' @description The function uses eigen values and eigen vectors to recreate a pd matrix.
#'
#' @param mat A matrix input.
#'
#' @return The function returns a modified matrix which is pd.
#'
makePD = function(mat){
  N = nrow(mat)
  HC = mat
  D = eigen(mat)
  E = D$values
  U = D$vectors
  v = as.numeric(E < 0)
  idx = which(v==1)
  m = sum(v) # number of negative values
  if(m > 0){
    S = sum(v*E)*2
    W = (S*S*100)+1
    P = min(abs(E[idx])) # smallest positive value
    for(i in idx){
      C = E[i]
      E[i] = P * (S-C)*(S-C)/W
    }
  }
  return(E)
  }


#' Function for the final model
#'
#' @description Obtain the final model according to a given criterion (AIC or BIC).
#'
#' @param fit The fitted models from grpreg.
#' @param crit The criteria for model selection.
#'
#' @return A list of the following objects
#' \item{coefficients}{A list containing coefficients of the count part and zero inflation part.}
#' \item{aic}{The AIC for the final model.}
#' \item{bic}{The BIC for the final model.}
#' \item{loglik}{The loglikelihood for the final model.}
#'
#' @export
#'
fit.final.func<-function(fit,crit)
{
  fit.final.idx<-which.min(get(crit)(fit))
  fit.final.aic<-AIC(fit)[fit.final.idx]
  fit.final.bic<-BIC(fit)[fit.final.idx]
  fit.final.loglik<-logLik(fit)[fit.final.idx]
  fit.final.coeff<-fit$beta[,fit.final.idx]

  fit.final.coeff.count<-fit.final.coeff[str_sub(names(fit.final.coeff),start = -5)=="count"]
  names(fit.final.coeff.count)<-substr(names(fit.final.coeff.count),1,nchar(names(fit.final.coeff.count))-6)
  fit.final.coeff.zero<-fit.final.coeff[str_sub(names(fit.final.coeff),start = -4)=="zero"]
  names(fit.final.coeff.zero)<-substr(names(fit.final.coeff.zero),1,nchar(names(fit.final.coeff.zero))-5)

  return(list(coefficients=list(count=fit.final.coeff.count,zero=fit.final.coeff.zero),aic=fit.final.aic, bic=fit.final.bic, loglik=fit.final.loglik))
}




#' A group regularized fit to the zero inflated count data.
#'
#' @description Fit zero inflated count data with a group regularization algorithm.
#'
#' @usage gooogle(data,xvars,zvars,yvar,group=1:ncol(data),samegrp.overlap=T,penalty=c("grLasso", "grMCP", "grSCAD", "gBridge"),dist=c("poisson","negbin"), nlambda=100, lambda,lambda.min=ifelse((nrow(data[,unique(c(xvars,zvars))])>ncol(data[,unique(c(xvars,zvars))])),1e-4,.05),lambda.max, crit="BIC",alpha=1, eps=.001, max.iter=1000, gmax=length(unique(group)),gamma=ifelse(penalty=="gBridge",0.5,ifelse(penalty == "grSCAD", 4, 3)), warn=TRUE)
#'
#'
#' @param data The data frame or matrix consisting of outcome and predictors.
#' @param xvars The vector of variable names to be included in count model.
#' @param zvars The vector of variable names for excess zero model.
#' @param yvar The outcome variable name.
#' @param group The vector of integers describing the grouping of the coefficients. For greatest efficiency and least ambiguity, it is best if group is a vector of consecutive integers. If there are coefficientss to be included in the model without being penalized, assign them to group 0 (or "0").
#' @param samegrp.overlap A logical argument. If TRUE (default) same grouping indices will be assigned to shared predictors in the count and degenerate distribution.
#' @param penalty The penalty to be applied in the model. For group level selection, one of "grLasso", "grMCP" or "grSCAD". For bi-level selection "gBridge" can be specified.
#' @param dist The distribution for count model - "poisson" for poisson or "negbin" for negative binomial.
#' @param nlambda The number of lambda values. Default is 100.
#' @param lambda A user specified sequence of lambda values.
#' @param lambda.min The smallest value for lambda, as a fraction of lambda.max. Default is .0001 if the number of observations is larger than the number of covariates and .05 otherwise.
#' @param lambda.max The maximum value for lambda (only needed for gBridge penalty).
#' @param crit The selection criteria for the best model. It can either be "AIC" or \code{BIC} (default).
#' @param alpha The tuning parameter for the balance between the group penalty and the L2 penalty, as in grpreg. Default value is 1.
#' @param max.iter Maximum number of iterations allowed.
#' @param gmax Maximum number of non-zero groups allowed.
#' @param gamma Tuning parameter of group MCP/SCAD. Default is 3 for MCP and 4 for SCAD.
#' @param warn A logical argument indicating whether this function gives warning in case of convergence issue.
#' @param eps The convergence threshhold, as in grpreg.
#'
#' @details The algorithm fits zero inflated count data to conduct variable selection in the presence of intrinsic grouping structure in the predictor set. Group wise penalties are considered for both count and zero abundance part of the mixture model where the likelihood is optimized using group level or bi-level co-ordinate descent algorithms.
#'
#' @return A list containing the following components is returned
#' \item{coefficients}{A list with two sets of coefficients corresponding to  count and zero inflation parts of the mixture model.}
#' \item{aic}{The AIC of the selected model.}
#' \item{bic}{The BIC of the selected model.}
#' \item{loglik}{The log-likelihood of the selected model.}
#'
#' @export
#' @examples
#' \dontrun{
#' ## Auto Insurance Claim Data
#' library(HDtweedie)
#' data("auto")
#' y<-auto$y
#' y<-round(y)
#' x<-auto$x
#' data<-cbind.data.frame(y,x)
#' group=c(rep(1,5),rep(2,7),rep(3,4),rep(4:14,each=3),15:21)
#' yvar<-names(data)[1]
#' xvars<-names(data)[-1]
#' zvars<-xvars
#'
#' ## ZIP regression
#' fit.poisson<-gooogle(data=data,yvar=yvar,xvars=xvars,zvars=zvars,group=group,samegrp.overlap=T,dist="poisson",penalty="gBridge")
#' fit.poisson$aic
#'
#' ## ZINB regression
#' fit.negbin<-gooogle(data=data,yvar=yvar,xvars=xvars,zvars=zvars,group=group,samegrp.overlap=T,dist="negbin",penalty="gBridge")
#' fit.negbin$aic
#' }
#'
#' @importFrom pscl zeroinfl
#' @importFrom grpreg grpreg gBridge
#' @importFrom stringr str_sub
#' @importFrom gamlss.dist dZIP dZINBI
#'

gooogle<-function(data,xvars,zvars,yvar,group=1:ncol(data),
                  samegrp.overlap=T,
                  penalty=c("grLasso", "grMCP", "grSCAD", "gBridge"),
                  dist=c("poisson","negbin"), nlambda=100, lambda,
                  lambda.min=ifelse((nrow(data[,unique(c(xvars,zvars))])>ncol(data[,unique(c(xvars,zvars))])),1e-4,.05), lambda.max, crit="BIC",
                  alpha=1, eps=.001, max.iter=1000, gmax=length(unique(group)),
                  gamma=ifelse(penalty=="gBridge",0.5,ifelse(penalty == "grSCAD", 4, 3)), warn=TRUE)
{
  ll.func<-function(beta.count,beta.zero,y,X,Z,dist)
  {
    y<-as.numeric(y[,1])

    if (dist=="negbin") theta<-a

    if (is.null(Z))
    {
      zgam<-rep(beta.zero,length(y))
    } else {
      zgam<-as.matrix(cbind(1,Z))%*%beta.zero
      zgam<-as.numeric(zgam[,1])
    }
    pzero<-exp(zgam)/(1+exp(zgam))

    xbet<-as.matrix(cbind(1,X))%*%beta.count
    xbet<-as.numeric(xbet[,1])
    mu<-exp(xbet)

    if (dist=="poisson"){
      ll<-try(sum(dZIP(y,mu=mu,sigma=pzero,log=T)),silent = T)
    } else {
      ll<-try(sum(dZINBI(y,mu=mu,sigma=1/theta,nu=pzero,log=T)),silent = T)
    }
    if (class(ll)=="try-error")
    {
      ll<-NA
    }
    return(ll)
  }

  y<-data.frame(data[,yvar])
  X<-data.frame(data[,xvars])
  pred.names<-union(xvars,zvars)
  names(X)<-paste(names(X),".count",sep="")
  names(y)<-yvar
  group.x<-group[which(pred.names %in% xvars)] #group of xvars
  n<-nrow(data)
  if (is.null(zvars)) #if there is no covariate in the zero model
  {
    Z<-NULL
    group.z<-NULL
    data<-cbind.data.frame(y,X)
    xvars<-names(X)
    fit.formula<-as.formula(paste(yvar,"~",paste(paste(xvars,collapse="+"),"|1"),sep=""))
  } else {
    Z<-data.frame(data[,zvars])
    names(Z)<-paste(names(Z),".zero",sep="")

    if (samegrp.overlap) # if X and Z assign same groups for shared covariates
    {
      group.z<-group[which(pred.names %in% zvars)]
    } else {
      group.z<-max(group.x)+group[which(pred.names %in% zvars)]
    }
    data<-cbind.data.frame(y,X,Z)
    xvars<-names(X)
    zvars<-names(Z)
    fit.formula<-as.formula(paste(yvar,"~",paste(paste(xvars,collapse="+"),"|",paste(zvars,collapse="+")),sep=""))
  }
  fit.zero <- zeroinfl(fit.formula, dist=dist,data = data)
  b2.mle<-c(fit.zero$coefficients$count[-1],fit.zero$coefficients$zero[-1])

  vcov<-fit.zero$vcov
  p<-length(xvars)

  if (dist=="poisson")
  {
    sigma.11<-vcov[c(1,p+2),c(1,p+2)]
    sigma.12<-vcov[c(1,p+2),-c(1,p+2)]
    sigma.22<-vcov[-c(1,p+2),-c(1,p+2)]

    vcov.bar<-sigma.22-t(sigma.12)%*%ginv(sigma.11)%*%sigma.12

    e<-eigen(vcov.bar)
    if(det(vcov.bar)>0) # in case vcov is not pd add small values to the diagonal
    {
      cov.star=e$vectors%*%diag(1/sqrt(abs(e$values)))%*%t(e$vectors)
    } else {
      cov.star=e$vectors%*%diag(1/sqrt(makePD(vcov.bar)))%*%t(e$vectors)
    }

    y.star<-cov.star%*%b2.mle # transformed y
    cov.star<-data.frame(cov.star) # scaled x matrix
    names(cov.star)<-c(xvars,zvars)

    ## reorder the group index and the covariates
    group<-c(group.x,group.z)
    u<-unique(group)
    cov.star.reordered<-rep(0,dim(cov.star)[1])
    for(i in 1:length(u))
    {
      index<-which(group==u[i])
      cov.star.reordered<-cbind(cov.star.reordered,cov.star[,index])
    }
    cov.star.reordered<-cov.star.reordered[,-1]
    group<-sort(group)

    if (penalty=="gBridge")
    {
      fit<-gBridge(X=cov.star.reordered, y=y.star, group=group, gamma=gamma)

    } else {
      fit=grpreg(cov.star.reordered,y.star,group=group,penalty=penalty,family="gaussian", nlambda=nlambda, lambda, lambda.min=lambda.min,  alpha=alpha, eps=eps, max.iter=max.iter,gmax=gmax, gamma=gamma, warn=warn)
    }
  } else if(dist=="negbin") {
    sigma.11<-vcov[c(1,p+2),c(1,p+2)]
    sigma.12<-vcov[c(1,p+2),-c(1,p+2)]
    sigma.22<-vcov[-c(1,p+2),-c(1,p+2)]

    vcov.bar<-try(sigma.22-t(sigma.12)%*%ginv(sigma.11)%*%sigma.12,silent=T)

    a<-1/fit.zero$theta
    var.a<-(fit.zero$SE.logtheta)^2*exp(2*log(a))

    e<-eigen(vcov.bar)

    if(round(det(vcov.bar),4)>0) # transform vcov to pd if it's not
    {
      cov.star=e$vectors%*%diag(1/sqrt(abs(e$values)))%*%t(e$vectors)
    } else {
      cov.star=e$vectors%*%diag(1/sqrt(makePD(vcov.bar)))%*%t(e$vectors)
    }
    y.star<-cov.star%*%b2.mle
    cov.star<-data.frame(cov.star) # scaled x matrix
    names(cov.star)<-c(xvars,zvars)

    group<-c(group.x,group.z) # 0 indicator correspond to disperison and intercept terms
    u<-unique(group)
    ## reorder the group index and the covariates
    cov.star.reordered<-rep(0,dim(cov.star)[1])
    for(i in 1:length(u))
    {
      index<-which(group==u[i])
      cov.star.reordered<-cbind(cov.star.reordered,cov.star[,index])
    }
    cov.star.reordered<-cov.star.reordered[,-1]
    group<-sort(group)

    if (penalty=="gBridge")
    {
      fit<-gBridge(X=cov.star.reordered, y=y.star, group=group, gamma=gamma)

    } else {
      fit=grpreg(cov.star.reordered,y.star,group=group,penalty=penalty,family="gaussian", nlambda=nlambda, lambda, lambda.min=lambda.min,  alpha=alpha, eps=eps, max.iter=max.iter, gmax=gmax, gamma=gamma, warn=warn)
    }
  } else {
    stop("Gooogle only works with Poisson or Negative Binomial distribution for count data")
  }

    fit.final<-fit.final.func(fit=fit,crit=crit)
    b2.final<-c(fit.final$coefficients$count,fit.final$coefficients$zero)
    b1.final<-sigma.12%*%solve(sigma.22)%*%(b2.mle-b2.final)
    fit.final.coeff.count<-c(intercept=b1.final[1],fit.final$coefficients$count)
    fit.final.coeff.zero<-c(intercept=b1.final[2],fit.final$coefficients$zero)

    dfc<-sum(fit.final.coeff.count!=0)
    dfz<-sum(fit.final.coeff.count!=0)
    ll<-ll.func(fit.final.coeff.count,fit.final.coeff.zero,y,X,Z,dist=dist)
    aic<--2*ll+2*(dfc+dfz)
    bic<--2*ll+log(n)*(dfc+dfz)

    if (dist=="poisson"){
      coeff.final=list(count=fit.final.coeff.count,zero=fit.final.coeff.zero)
    } else {
      coeff.final=list(count=fit.final.coeff.count,zero=fit.final.coeff.zero,dispersion=a)
    }

    return(list(coefficients=coeff.final,aic=aic, bic=bic, loglik=ll))
}
