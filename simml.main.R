
#' Single-index models with multiple-links (main function)
#'
#' \code{simml} is the wrapper function for Single-index models with multiple-links (SIMML).
#' The function estimates a linear combination (a single-index) of covariates X, and models the treatment-specific outcome y, via treatment-specific nonparametrically-defined link functions.
#'
#'
#' SIMML captures the effect of covariates via a single-index and their interaction with the treatment via nonparametric link functions.
#' Interaction effects are determined by distinct shapes of the link functions.
#' The estimated single-index is useful for comparing differential treatment efficacy.
#' The resulting \code{simml} object can be used to estimate an optimal treatment decision rule
#' for a new patient with pretreatment clinical information.
#'
#' @param y  a n-by-1 vector of treatment outcomes; y is a member of the exponential family; any distribution supported by \code{mgcv::gam}; y can also be an ordinal categorial response with \code{R} categories taking a value from 1 to \code{R}.
#' @param A  a n-by-1 vector of treatment variable; each element is assumed to take a value in a finite discrete space.
#' @param X  a n-by-p matrix of baseline covarates.
#' @param Xm  a n-by-q design matrix associated with an X main effect model; the defult is \code{NULL} and it is taken as a vector of zeros
#' @param family  specifies the distribution of y; e.g., "gaussian", "binomial", "poisson"; can be any family supported by \code{mgcv::gam}; can also be "ordinal", for an ordinal categorical response y.
#' @param R   the number of response categories for the case of family = "ordinal".
#' @param bs basis type for the treatment (A) and single-index joint effect; the defult is "ps" (p-splines); any basis supported by \code{mgcv::gam} can be used, e.g., "cr" (cubic regression splines); see \code{mgcv::s} for detail.
#' @param k  basis dimension for the spline-type-represented treatment-specific link functions.
#' @param sp  smoothing paramter for the treatment-specific link functions; if \code{NULL}, then estimated from the data.
#' @param linear.link  if \code{TRUE}, the link function is restricted to be linear.
#' @param method  the smoothing parameter estimation method; "GCV.Cp" to use GCV for unknown scale parameter and Mallows' Cp/UBRE/AIC for known scale; any method supported by \code{mgcv::gam} can be used.
#' @param gamma  increase this beyond 1 to produce smoother models. \code{gamma} multiplies the effective degrees of freedom in the GCV or UBRE/AIC (see \code{mgcv::gam} for detail); the default is 1.
#' @param aug a n-by-1 additional augmentation vector associated with the X main effect; the default is \code{NULL} and it is taken as a vector of zeros
#' @param rho    a tuning parameter associated with the additional augmentation vector \code{aug}; the default is 0.
#' @param beta.ini  an initial value for \code{beta.coef}; a p-by-1 vector; the defult is \code{NULL}, in which case a linear model estimate is used.
#' @param ind.to.be.positive  for identifiability of the solution \code{beta.coef}, the user can restrict the jth (e.g., j=1) component of \code{beta.coef} to be positive; by default, we match the "overall" sign of \code{beta.coef} with that of the linear estimate (i.e., the initial estimate), by restricting the inner product between the two to be positive.
#' @param scale.si.01 if \code{TRUE}, re-scale the index coefficients to restrict the index to the interval [0,1]; in such a case, an intercept term is induced.
#' @param max.iter  an integer specifying the maximum number of iterations for \code{beta.coef} update.
#' @param eps.iter a value specifying the convergence criterion of algorithm.
#' @param trace.iter if \code{TRUE}, trace the estimation process and print the differences in \code{beta.coef}.
#' @param pen.order 0 indicates the ridge penalty; 1 indicates the 1st difference penalty; 2 indicates the 2nd difference penalty, used in a penalized least squares (LS) estimation of \code{beta.coef}.
#' @param lambda  a regularization parameter associated with the penalized LS for \code{beta.coef} update; the default is 0, and the index coefficients are not penalized.
#' @param center.X   if \code{TRUE}, center X to have zero mean.
#' @param scale.X    if \code{TRUE}, scale X to have unit variance.
#' @param ortho.constr  separates the interaction effects from the main effect (without this, the interaction effect can be confounded by the main effect; the default is \code{TRUE}.
#' @param si.main.effect  if \code{TRUE}, once the convergence in the estimates of \code{beta.coef} is reached, include the main effect associated with the fitted single-index (beta.coef'X) to the final fit; the default is \code{FALSE}.
#' @param random.effect  if \code{TRUE}, as part of the main effects, the user can incorporate z-specific random intercepts.
#' @param z  a factor that specifies the random intercepts when \code{random.effect = TRUE}.
#' @param plots if \code{TRUE}, produce a plot for the estimated effect contrast (for binary treatment cases) (on a linear predictor scale).
#' @param bootstrap if \code{TRUE}, compute bootstrap confidence intervals for the single-index coefficients, \code{beta.coef}; the default is \code{FALSE}.
#' @param boot.conf  a value specifying the confidence level of the bootstrap confidence intervals; the defult is \code{boot.conf = 0.95}.
#' @param nboot  when \code{bootstrap=TRUE}, a value specifying the number of bootstrap replications.
#' @param seed  when  \code{bootstrap=TRUE}, randomization seed used in bootstrap resampling.
#'
#' @return a list of information of the fitted SIMML including
#'  \item{beta.coef}{ the estimated single-index coefficients.} \item{g.fit}{a \code{mgcv:gam} object containing information about the estimated treatment-specific link functions.} \item{beta.ini}{the initial value used in the estimation of \code{beta.coef}} \item{beta.path}{solution path of \code{beta.coef} over the iterations} \item{d.beta}{records the change in \code{beta.coef} over the solution path, \code{beta.path}} \item{scale.X}{sd of pretreatment covariates X} \item{center.X}{mean of pretreatment covariates X} \item{L}{number of different treatment options} \item{p}{number of pretreatment covariates X} \item{n}{number of subjects} \item{boot.ci}{(1-boot.alpha/2) percentile bootstrap CIs (LB, UB) associated with \code{beta.coef}}
#'
#' @author Park, Petkova, Tarpey, Ogden
#' @import mgcv stats graphics
#' @seealso \code{pred.simml},  \code{fit.simml}
#' @export
#'
#' @examples
#'
#'
#' family <- "gaussian"   #"poisson"
#' delta = 1              # moderate main effect
#' s=2                    # if s=2 (s=1), a nonlinear (linear) contrast function
#' n=500                  # number of subjects
#' p=10                   # number of pretreatment covariates
#'
#' # generate training data
#' data <- generate.data(n= n, p=p, delta = delta, s= s, family = family)
#' data$SNR  # the ratio of interactions("signal") vs. main effects("noise")
#' A <- data$A
#' y <- data$y
#' X <- data$X
#'
#' # generate testing data
#' data.test <- generate.data(n=10^5, p=p, delta = delta,  s= s, family = family)
#' A.test <- data.test$A
#' y.test <- data.test$y
#' X.test <- data.test$X
#' data.test$value.opt     # the optimal "value"
#'
#'
#' # fit SIMML
#' #1) SIMML without X main effect
#' simml.obj1 <- simml(y, A, X, family = family)
#'
#' #2) SIMML with X main effect (estimation efficiency for the g term of SIMML can be improved)
#' simml.obj2 <- simml(y, A, X, Xm = X, family = family)
#'
#'
#' # apply the estimated SIMML to the testing set and obtain treatment assignment rules.
#' simml.trt.rule1 <- pred.simml(simml.obj1, newX= X.test)$trt.rule
#' # "value" estimation (estimated by IPWE)
#' simml.value1 <-  mean(y.test[simml.trt.rule1 == A.test])
#' simml.value1
#'
#' simml.trt.rule2 <- pred.simml(simml.obj2, newX= X.test)$trt.rule
#' simml.value2 <-  mean(y.test[simml.trt.rule2 == A.test])
#' simml.value2
#'
#' # compare these to the optimal "value"
#' data.test$value.opt
#'
#'
#'
#' # fit MC (modified covariates) model of Tien et al 2014
#' n.A <- summary(as.factor(A)); pi.A <- n.A/sum(n.A)
#' mc  <- (as.numeric(A) + pi.A[1] -2) *cbind(1, X)  # 0.5*(-1)^as.numeric(A) *cbind(1, X)
#' mc.coef  <-  coef(glm(y ~ mc, family =  family))
#' mc.trt.rule <- (cbind(1, X.test) %*% mc.coef[-1] > 0) +1
#' # "value" estimation (estimated by IPWE)
#' mc.value  <-  mean(y.test[mc.trt.rule == A.test])
#' mc.value
#'
#'
#' # visualization of the estimated link functions of SIMML
#' simml.obj1$beta.coef        # estimated single-index coefficients
#' g.fit <- simml.obj1$g.fit   # estimated trt-specific link functions; "g.fit" is a mgcv::gam object.
#' #plot(g.fit)
#'
#'\donttest{
#' # can improve visualization by using the package "mgcViz"
#' #install.packages("mgcViz")
#' # mgcViz depends on "rgl". "rgl" depends on XQuartz, which you can download from xquartz.org
#' #library(mgcViz)
#' # transform the "mgcv::gam" object to a "mgcViz" object (to improve visualization)
#' g.fit <- getViz(g.fit)
#'
#' plot1  <- plot( sm(g.fit,1) )  # for treatment group 1
#' plot1 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
#'   l_ciLine(mul = 5, colour = "blue", linetype = 2) +
#'   l_points(shape = 19, size = 1, alpha = 0.1) +
#'   xlab(expression(paste("z = ", alpha*minute, "x")))  +  ylab("y") +
#'   ggtitle("Treatment group 1 (Trt =1)") +  theme_classic()
#'
#' plot2 <- plot( sm(g.fit,2) )   # for treatment group 2
#' plot2 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
#'   l_ciLine(mul = 5, colour = "blue", linetype = 2) +
#'   l_points(shape = 19, size = 1, alpha = 0.1) +
#'   xlab(expression(paste("z = ", alpha*minute, "x"))) +ylab("y") +
#'   ggtitle("Treatment group 2 (Trt =2)") + theme_classic()
#'
#'
#' trans = function(x) x + g.fit$coefficients[2]
#' plotDiff(s1 = sm(g.fit, 2), s2 = sm(g.fit, 1), trans=trans) +  l_ciPoly() +
#'   l_fitLine() + geom_hline(yintercept = 0, linetype = 2) +
#'   xlab(expression(paste("z = ", alpha*minute, "x")) ) +
#'   ylab("(Treatment 2 effect) - (Treatment 1 effect)") +
#'   ggtitle("Contrast between two treatment effects") +
#'   theme_classic()
#'
#'
#' # yet another way of visualization, using ggplot2
#' #library(ggplot2)
#' dat  <- data.frame(y= simml.obj1$g.fit$model$y,
#'                    x= simml.obj1$g.fit$model$single.index,
#'                    Treatment= simml.obj1$g.fit$model$A)
#' g.plot<- ggplot(dat, aes(x=x,y=y,color=Treatment,shape=Treatment,linetype=Treatment))+
#'    geom_point(aes(color=Treatment, shape=Treatment), size=1, fill="white") +
#'    scale_colour_brewer(palette="Set1", direction=-1) +
#'    xlab(expression(paste(beta*minute,"x"))) + ylab("y")
#' g.plot + geom_smooth(method=gam, formula= y~ s(x, bs=simml.obj1$bs, k=simml.obj1$k),
#'                      se=TRUE, fullrange=TRUE, alpha = 0.35)
#'}
#'
#'\donttest{
#' # can obtain bootstrap CIs for beta.coef.
#' simml.obj <- simml(y,A,X,Xm=X, family=family,bootstrap=TRUE,nboot=15)  #nboot=500.
#' simml.obj$beta.coef
#' round(simml.obj$boot.ci,3)
#'
#' # compare the estimates to the true beta.coef.
#' data$true.beta
#'}
#'
#'
#'# an application to data with ordinal categorical response
#'dat <- ordinal.data(n=500, p=5, R = 11,  # 11 response levels
#'                    s = "nonlinear",     # nonlinear interactions
#'                    delta = 1)
#'dat$SNR
#'y <- dat$y  # ordinal response
#'X <- dat$X  # X matrix
#'A <- dat$A  # treatment
#'dat$true.beta  # the "true" single-index coefficient
#'
#'\donttest{
#'# 1) fit a cumulative logit simml, with a flexible link function
#'res <-  simml(y,A,X, family="ordinal", R=11)
#'res$beta.coef  # single-index coefficients.
#'res$g.fit$family$getTheta(TRUE)  # the estimated R-1 threshold values.
#'}
#'# 2) fit a cumulative logit simml, with a linear link function
#'res2 <-  simml(y,A,X, family="ordinal", R=11, linear.link = TRUE)
#'res2$beta.coef  # single-index coefficients.
#'
#'\donttest{
#'family = mgcv::ocat(R=11)  # ocat: ordered categorical response family, with R categories.
#'# the treatment A's effect.
#'tmp <- mgcv::gam(y ~ A, family =family)
#'exp(coef(tmp)[2])  #odds ratio (OR) comparing treatment A=2 vs. A=1.
#'
#'ind2 <- pred.simml(res)$trt.rule ==2  # subgroup recommended with A=2 under SIMML ITR
#'tmp2 <- mgcv::gam(y[ind2] ~ A[ind2], family = family)
#'exp(coef(tmp2)[2]) #OR comparing treatment A=2 vs. A=1, for subgroup recommended with A=2
#'
#'ind1 <- pred.simml(res)$trt.rule ==1  # subgroup recommended with A=1 under SIMML ITR
#'tmp1 <- mgcv::gam(y[ind1] ~ A[ind1], family = family)
#'exp(coef(tmp1)[2]) #OR comparing treatment A=2 vs. A=1, for subgroup recommended with A=2
#'}
simml <- function(y,   # a n x 1 vector of observed responses; if family = "ordinal", the observed categories should be numerically coded: 1, 2, 3, ... up to the number of categories.
                  A,   # a n x 1 vector of observed discrete treatments; either numeric or factor.
                  X,   # a n x p matrix of covariates.
                  Xm = NULL,  # a n-by-q design matrix associated with the X main effect.
                  aug= NULL,
                  family = "gaussian",  # specifies the distribution of y; e.g., "gaussian", "binomial", "poisson"; can be any family supported by \code{mgcv::gam}; can also be "ordinal" for ordinal catagorical response.
                  R = NULL,    # the number of catergories in the ordinal response y; only needed for the case family = "ordinal".
                  bs = "cr",   # this specifies basis for the treatment-specific link functions; the default is "ps" (p-splines); any basis supported by \code{mgcv::gam} can be used, e.g., "cr" (cubic regression splines).
                  k = 8,       # number of basis function for the treatment-specific link functions
                  sp = NULL,   # smoothing paramters for the treatment-specific link functions; if NULL, then estimated from the data.
                  linear.link = FALSE,
                  method = "GCV.Cp", # the smoothing parameter estimation method; can be "REML" and "ML"; see mgcv::gam for detail.
                  gamma = 1,   # increase this beyond 1 to produce smoother models; gamma multiplies the effective degrees of freedom in the GCV or UBRE/AIC; see mgcv::gam for details.
                  rho= 0,
                  beta.ini=NULL,  # an initial value for beta.coef; a p-by-1 vector; the defult is NULL.
                  ind.to.be.positive = NULL, # for identifiability of the solution beta.coef, we restrict the jth component of beta.coef to be positive; by default j=1.
                  scale.si.01 = FALSE,
                  max.iter = 20,  # an integer specifying the maximum number of iterations for \code{beta.coef} update.
                  eps.iter = 0.01,   # a value specifying the convergence criterion of algorithm.
                  trace.iter = TRUE,
                  lambda = 0,     # a regularziation parameter associated with a penalized least squares estimation of beta.coef.
                  pen.order = 0,  # 0 indicates the ridge penalty; 1 indicates the 1st difference penalty; 2 indicates the 2nd difference penalty, used in a penalized least squares estimation of beta.coef.
                  scale.X = TRUE,
                  center.X = TRUE,
                  ortho.constr = TRUE,
                  si.main.effect = FALSE,
                  random.effect = FALSE,
                  z = NULL, plots = FALSE,
                  bootstrap = FALSE, nboot = 200, boot.conf = 0.95, seed = 1357)   # these are for bootstrap confidnce intervals.
{
  
  simml.obj <- fit.simml(y=y,A=A,X=X, Xm=Xm, aug=aug, rho=rho,
                         family=family, R=R, bs=bs, k=k, sp=sp,
                         linear.link=linear.link, method=method, max.iter=max.iter,
                         eps.iter=eps.iter,lambda=lambda, ind.to.be.positive = ind.to.be.positive,
                         scale.si.01=scale.si.01, center.X= center.X, scale.X= scale.X,
                         ortho.constr = ortho.constr, si.main.effect=si.main.effect,
                         random.effect=random.effect, z=z, plots = plots,
                         beta.ini=beta.ini, trace.iter=trace.iter, pen.order=pen.order)
  
  boot.mat = boot.ci <- NULL
  if(bootstrap){
    set.seed(seed)
    indices <- 1:simml.obj$n
    if(is.null(Xm)) Xm <- rep(0,simml.obj$n)
    Xm <- as.matrix(Xm)
    if(is.null(z))  z <-  rep(0,simml.obj$n)
    if(is.null(aug))  aug <-  rep(0,simml.obj$n)
    if(simml.obj$scale.si.01){
      boot.mat <- matrix(0, nboot, simml.obj$p+1)
    }else{
      boot.mat <- matrix(0, nboot, simml.obj$p)
    }
    
    for(i in 1:nboot){
      boot.indices <- sample(indices, simml.obj$n, replace = TRUE)
      boot.beta  <- fit.simml(y=y[boot.indices], A = A[boot.indices], X = X[boot.indices,],
                              Xm = Xm[boot.indices,], aug = aug[boot.indices], rho =rho,
                              family=family, R=R, bs =bs, k = k, sp= sp, linear.link=linear.link, method= method,
                              beta.ini = beta.ini, random.effect= random.effect, z=z[boot.indices], plots = plots,
                              ind.to.be.positive=ind.to.be.positive, scale.si.01=scale.si.01,
                              pen.order = pen.order, lambda = lambda, max.iter = max.iter, gamma=gamma, trace.iter=trace.iter,
                              center.X= center.X, scale.X= scale.X,  ortho.constr = ortho.constr,
                              si.main.effect= si.main.effect)$beta.coef
      
      if(simml.obj$beta.coef %*% boot.beta > simml.obj$beta.coef %*% (-boot.beta)){
        boot.mat[i,] <-  boot.beta
      }else{
        boot.mat[i,] <- -boot.beta
      }
      if(trace.iter) print(i)
    }
    
    var.t0 <- apply(boot.mat, 2, var)
    boot.ci <- cbind(simml.obj$beta.coef-qnorm((1+boot.conf)/2)*sqrt(var.t0),
                     simml.obj$beta.coef+qnorm((1+boot.conf)/2)*sqrt(var.t0))
    
    boot.ci <- cbind(simml.obj$beta.coef, boot.ci, (boot.ci[,1] > 0 | boot.ci[,2] < 0) )
    colnames(boot.ci) <- c("coef", "LB", "UB", " ***")
    if(simml.obj$scale.si.01){
      rownames(boot.ci) <- c("Int.", colnames(X))
    }else{
      rownames(boot.ci) <- colnames(X)
    }
  }
  simml.obj$boot.mat <- boot.mat
  simml.obj$boot.ci <- boot.ci
  
  return(simml.obj)
}




#' Single-index models with multiple-links (workhorse function)
#'
#' \code{fit.simml} is the workhorse function for Single-index models with multiple-links (SIMML).
#' The function estimates a linear combination (a single-index) of covariates X, and models the treatment-specific outcome y, via treatment-specific nonparametrically-defined link functions.
#'
#' SIMML captures the effect of covariates via a single-index and their interaction with the treatment via nonparametric link functions.
#' Interaction effects are determined by distinct shapes of the link functions.
#' The estimated single-index is useful for comparing differential treatment efficacy.
#' The resulting \code{simml} object can be used to estimate an optimal treatment decision rule
#' for a new patient with pretreatment clinical information.
#'
#' @param y  a n-by-1 vector of treatment outcomes; y is a member of the exponential family; any distribution supported by \code{mgcv::gam}; y can also be an ordinal categorial response with \code{R} categories taking a value from 1 to \code{R}.
#' @param A  a n-by-1 vector of treatment variable; each element is assumed to take a value on a continuum.
#' @param X  a n-by-p matrix of baseline covarates.
#' @param Xm  a n-by-q design matrix associated with an X main effect model; the defult is \code{NULL} and it is taken as a vector of zeros
#' @param family  specifies the distribution of y; e.g., "gaussian", "binomial", "poisson"; can be any family supported by \code{mgcv::gam}; can also be "ordinal", for an ordinal categorical response y.
#' @param R   the number of response categories for the case of family = "ordinal".
#' @param bs basis type for the treatment (A) and single-index domains, respectively; the defult is "ps" (p-splines); any basis supported by \code{mgcv::gam} can be used, e.g., "cr" (cubic regression splines); see \code{mgcv::s} for detail.
#' @param k  basis dimension for the treatment (A) and single-index domains, respectively.
#' @param sp  smoothing paramter for the treatment-specific link functions; if \code{NULL}, then estimated from the data.
#' @param linear.link  if \code{TRUE}, the link function is restricted to be linear.
#' @param method  the smoothing parameter estimation method; "GCV.Cp" to use GCV for unknown scale parameter and Mallows' Cp/UBRE/AIC for known scale; any method supported by \code{mgcv::gam} can be used.
#' @param gamma  increase this beyond 1 to produce smoother models. \code{gamma} multiplies the effective degrees of freedom in the GCV or UBRE/AIC (see \code{mgcv::gam} for detail); the default is 1.
#' @param aug a n-by-1 additional augmentation vector associated with the X main effect; the default is \code{NULL} and it is taken as a vector of zeros
#' @param rho    a tuning parameter associated with the additional augmentation vector \code{aug}; the default is 0.
#' @param beta.ini  an initial value for \code{beta.coef}; a p-by-1 vector; the defult is \code{NULL}, in which case a linear model estimate is used.
#' @param ind.to.be.positive  for identifiability of the solution \code{beta.coef}, the user can restrict the jth (e.g., j=1) component of \code{beta.coef} to be positive; by default, we match the "overall" sign of \code{beta.coef} with that of the linear estimate (i.e., the initial estimate), by restricting the inner product between the two to be positive.
#' @param scale.si.01 if \code{TRUE}, re-scale the index coefficients to restrict the index to the interval [0,1]; in such a case, an intercept term is induced.
#' @param max.iter  an integer specifying the maximum number of iterations for \code{beta.coef} update.
#' @param eps.iter a value specifying the convergence criterion of algorithm.
#' @param trace.iter if \code{TRUE}, trace the estimation process and print the differences in \code{beta.coef}.
#' @param pen.order 0 indicates the ridge penalty; 1 indicates the 1st difference penalty; 2 indicates the 2nd difference penalty, used in a penalized least squares (LS) estimation of \code{beta.coef}.
#' @param lambda  a regularization parameter associated with the penalized LS for \code{beta.coef} update.
#' @param center.X   if \code{TRUE}, center X to have zero mean.
#' @param scale.X    if \code{TRUE}, scale X to have unit variance.
#' @param ortho.constr  separates the interaction effects from the main effect (without this, the interaction effect can be confounded by the main effect; the default is \code{TRUE}.
#' @param si.main.effect  if \code{TRUE}, once the convergence in the estimates of \code{beta.coef} is reached, include the main effect associated with the fitted single-index (beta.coef'X) to the final fit; the default is \code{FALSE}.
#' @param random.effect  if \code{TRUE}, as part of the main effects, the user can incorporate z-specific random intercepts.
#' @param z  a factor that specifies the random intercepts when \code{random.effect = TRUE}.
#' @param plots if \code{TRUE}, produce a plot for the estimated effect contrast (for binary treatment cases) (on a linear predictor scale).
#'
#' @return a list of information of the fitted SIMML including
#'  \item{beta.coef}{ the estimated single-index coefficients.} \item{g.fit}{a \code{mgcv:gam} object containing information about the estimated treatment-specific link functions.} \item{beta.ini}{the initial value used in the estimation of \code{beta.coef}} \item{beta.path}{solution path of \code{beta.coef} over the iterations} \item{d.beta}{records the change in \code{beta.coef} over the solution path, \code{beta.path}} \item{scale.X}{sd of pretreatment covariates X} \item{center.X}{mean of pretreatment covariates X} \item{L}{number of different treatment options} \item{p}{number of pretreatment covariates X} \item{n}{number of subjects} \item{boot.ci}{(1-boot.alpha/2) percentile bootstrap CIs (LB, UB) associated with \code{beta.coef}}
#'
#' @author Park, Petkova, Tarpey, Ogden
#' @import mgcv
#' @seealso \code{pred.simml},  \code{simml}
#' @export
#'
fit.simml <- function(y,   # a n x 1 vector of observed responses; if family = "ordinal", the observed categories should be numerically coded: 1, 2, 3, ... up to the number of categories.
                      A,   # a n x 1 vector of observed discrete treatments; either numeric or factor.
                      X,   # a n x p matrix of covariates.
                      Xm = NULL,  # a n-by-q design matrix assocaited with the X main effect.
                      aug = NULL,
                      rho = 0,
                      family = "gaussian",  # specifies the distribution of y; e.g., "gaussian", "binomial", "poisson"; can be any family supported by \code{mgcv::gam}; can also be "ordinal" for ordinal catagorical response.
                      R = NULL,    # the number of catergories in the ordinal response y; only needed for the case family = "ordinal".
                      bs = "ps",   # this specifies basis for the treatment-specific link functions; the default is "ps" (p-splines); any basis supported by \code{mgcv::gam} can be used, e.g., "cr" (cubic regression splines).
                      k = 8,       # number of basis function for the treatment-specific link functions
                      sp = NULL,   # smoothing paramters for the treatment-specific link functions; if NULL, then estimated from the data.
                      linear.link= FALSE,
                      method = "GCV.Cp", # the smoothing parameter estimation method; can be "REML" and "ML"; see mgcv::gam for detail.
                      gamma = 1,   # increase this beyond 1 to produce smoother models; gamma multiplies the effective degrees of freedom in the GCV or UBRE/AIC; see mgcv::gam for details.
                      max.iter = 20,        # an integer specifying the maximum number of iterations for \code{beta.coef} update.
                      eps.iter = 0.01,   # a value specifying the convergence criterion of algorithm.
                      trace.iter = TRUE,
                      ind.to.be.positive = NULL, # for identifiability of the solution beta.coef, we restrict the jth component of beta.coef to be positive; by default j=1.
                      scale.si.01 = FALSE,
                      lambda = 0,     # a regularziation parameter associated with a penalized least squares estimation of beta.coef.
                      pen.order = 0,  # 0 indicates the ridge penalty; 1 indicates the 1st difference penalty; 2 indicates the 2nd difference penalty, used in a penalized least squares estimation of beta.coef.
                      scale.X = TRUE,
                      center.X = TRUE,
                      ortho.constr = TRUE,
                      beta.ini=NULL,  # an initial value for beta.coef; a p-by-1 vector; the defult is NULL.
                      si.main.effect = FALSE,
                      random.effect = FALSE,
                      z= NULL, plots=FALSE)
{
  
  n <- length(y)
  p <- ncol(X)
  A <- as.factor(A)        # A is a categorical variable.
  prob.a <- summary(A)/n   # randomization probability for A.
  L <- length(levels(A))
  
  ## Center and scale X
  Xc <- scale(X, center = center.X, scale = scale.X)
  X.center <- attr(Xc, "scaled:center")
  X.scale <- attr(Xc, "scaled:scale")
  
  ## If not provided by the user, the covariate matrix associated with the X main effect is set to be a zero matrix.
  if(is.null(Xm)){
    if(is.null(aug)){ Xm <- rep(0,n)}
    else{
      aug.scaled <- (aug - min(aug))/(max(aug)-min(aug))
      Xm <- aug.scaled
    }
  }
  Xm <- as.matrix(Xm)
  
  ## If not provided by the user, the efficiency augmentation vector (corresponding to the X main effect) is set to be a zero vector.
  if(is.null(aug)) aug <- rep(0,n)
  
  
  ## special case of an ordinal categorical response
  if(family=="ordinal"){
    if(is.null(R)) R <- length(unique(y))
    family=ocat(R=R)
  }
  
  ## initialize the single-index coefficient and the single-index.
  if(is.null(beta.ini)){
    Ac <- as.numeric(A)-mean(as.numeric(A))
    if(random.effect){
      tmp <- gam(y~ A + Xm + s(z, bs="re")+ Ac:Xc, family = family)$coef
    }else{
      tmp <- gam(y~ A + Xm + Ac:Xc, family = family)$coef
    }
    beta.ini <- tmp[grep("Ac:" ,names(tmp))]
    beta.ini[which(is.na(beta.ini))] <- 0
    names(beta.ini) <- colnames(Xc)
  }
  gamma.magnitude <- sd(Ac)*sqrt(sum(beta.ini^2))
  
  beta.coef <- beta.ini/sqrt(sum(beta.ini^2))  # enforce unit L2 norm
  if(!is.null(ind.to.be.positive)){
    if(beta.coef[ind.to.be.positive] < 0) beta.coef <- -1*beta.coef      # for the (sign) identifiability
  }
  single.index <- as.vector(Xc %*% beta.coef)
  si.min <- min(single.index)
  si.max <- max(single.index)
  si.ran <- si.max - si.min
  
  
  ## initialize the link function
  if(linear.link){  # the heavy penalization essentially restricts the link functions to be linear.
    
    if(scale.si.01){
      beta.coef <- c("Int."= -range(single.index/si.ran)[1], beta.coef/si.ran)
      single.index <- (single.index - si.min)/si.ran   # cbind(1, Xc) %*% beta.coef
    }
    
    if(random.effect){
      g.fit <- gam(y~ A + Xm + s(z, bs="re")+ s(single.index, by=A, bs="cr", k=3, sp=10^6), family=family)
    }else{
      g.fit <- gam(y~ A + Xm + s(single.index, by=A, bs="cr", k=3, sp=10^6), family=family)
    }
  }else{  # note the initialization strategy: first estimate the link without penalization (sp=0); otherwise it is easy to get trapped in a local optimum in which the smooth is linear.
    if(random.effect){
      g.fit <- gam(y~ A + Xm+ s(z, bs="re") + s(single.index, by=A, bs=bs, k=k, sp=sp), family=family, gamma=gamma, method=method)
    }else{
      g.fit <- gam(y~ A + Xm + s(single.index, by=A, bs=bs, k=k, sp=sp), family=family, gamma=gamma, method=method)
    }
  }
  if(ortho.constr){
    tmp <- g.fit$coefficients
    smooths.coef <- tmp[grep("single.index):A" ,names(tmp))]   # basis coefficients associated with the link functions
    B <- matrix(smooths.coef, ncol = L)
    B <- B - apply(B, 1, function(x) weighted.mean(x, w = prob.a))
    g.fit$coefficients[grep("single.index):A" ,names(tmp))] <- as.vector(B)  # impose the orthogonality constraint on the link funcitons
  }
  
  beta.path <- beta.coef
  d.beta = beta.fit <- NULL
  if(!linear.link)  # if we want a flexible link, we iterate between the estimation of g.fit (i.e., link) and beta.coef (i.e., single.index), based on a local linear approxiamtion of the link.
  {
    
    # take the 1st deriavative of the treatment-specific smooths, w.r.t. the single.index.
    g.der <- der.link(g.fit)
    
    # Specify a penalty matrix associated with the penalized least squares for estimating beta.coef.
    D <- diag(p);  if(pen.order != 0)  for(j in 1:pen.order) D <- diff(D);
    Pen <- sqrt(lambda)*D
    
    for (it in 2:max.iter) {
      ## Update beta.coef and intercept through lsfit
      # adjusted responses, adjusted for the nonlinearity associated with the smooth
      y.star <- residuals(g.fit, "working") + g.der*single.index - rho*aug
      X.star <- diag(g.der)%*%Xc
      nix <- rep(0, nrow(D))
      X.p <- rbind(X.star, Pen)
      y.p <- c(y.star, nix)
      beta.fit <- lsfit(X.p, y.p, wt = c(g.fit$weights, (nix+1)))
      
      # for the identifiability
      beta.new <- beta.fit$coef[-1]/sqrt(sum(beta.fit$coef[-1]^2))
      if(beta.ini %*%beta.new  < 0)  beta.new <- -1*beta.new
      beta.path <- rbind(beta.path, beta.new)
      
      ## Check the convergence of beta
      d.beta   <- c(d.beta, sum((beta.new-beta.coef)^2))
      
      if(trace.iter){
        cat("iter:", it, " "); cat(" difference in beta: ", d.beta[(it-1)], "\n")
      }
      if (d.beta[(it-1)] < eps.iter)
        break
      beta.coef <- beta.new
      single.index <- as.vector(Xc %*% beta.coef)
      
      if(random.effect){
        g.fit <- gam(y~ A + Xm+ s(z, bs="re") + s(single.index, by=A, bs=bs, k=k, sp=sp), family=family, gamma=gamma, method=method)
      }else{
        g.fit <- gam(y~ A + Xm + s(single.index, by=A, bs=bs, k=k, sp=sp), family=family, gamma=gamma, method=method)
      }
      
      tmp <- g.fit$coefficients
      if(ortho.constr){
        smooths.coef <- tmp[grep("single.index):A" ,names(tmp))]   # basis coefficients accosited with the link functions
        B <- matrix(smooths.coef, ncol = L)
        B <- B - apply(B, 1, function(x) weighted.mean(x, w = prob.a))
        g.fit$coefficients[grep("single.index):A" ,names(tmp))] <- as.vector(B)  # impose the orthogonality constraint on the link funcitons
      }
      # take the 1st deriavative of the treatment-specific smooths, w.r.t. the single.index.
      g.der <- der.link(g.fit)
    }
    
    si.min <- min(single.index)
    si.max <- max(single.index)
    si.ran <- si.max - si.min
    if(scale.si.01){
      beta.coef    <- c("Int."= -range(single.index/si.ran)[1], beta.coef/si.ran)
      single.index <- (single.index - si.min)/si.ran
      if(random.effect){
        g.fit <- gam(y~ A + Xm+ s(z, bs="re") + s(single.index, by=A, bs=bs, k=k, sp=sp), family=family, gamma=gamma, method=method)
      }else{
        g.fit <- gam(y~ A + Xm + s(single.index, by=A, bs=bs, k=k, sp=sp), family=family, gamma=gamma, method=method)
      }
      tmp <- g.fit$coefficients
      if(ortho.constr){
        smooths.coef <- tmp[grep("single.index):A" ,names(tmp))]   # basis coefficients accosited with the link functions
        B <- matrix(smooths.coef, ncol = L)
        B <- B - apply(B, 1, function(x) stats::weighted.mean(x, w = prob.a))
        g.fit$coefficients[grep("single.index):A" ,names(tmp))] <- as.vector(B)  # impose the orthogonality constraint on the link funcitons
      }
    }
    
  }
  
  if(si.main.effect) g.fit$coefficients <- tmp
  AIC <- g.fit$aic + 2*ncol(X)
  
  if(plots & L==2){
    A.ord <- ordered(A)
    if(scale.si.01){
      si.grid <- seq(0.025, 0.975,  length.out = 100)
    }else{
      si.grid <- seq(si.min, si.max, length.out = 100)
    }
    Xm.mean <- apply(Xm, 2, mean)
    tmp <- matrix(rep(Xm.mean, each=100), 100)
    colnames(tmp)  <- colnames(Xm)
    if(random.effect){
      gam.fit <- gam(y~ A.ord + s(single.index, bs=bs, k=k, sp=sp) +
                       Xm+ s(z, bs="re")+ s(single.index, by=A.ord, bs=bs, k=k, sp=sp), family=family, method=method)
      dat <- list(Xm = tmp, single.index=si.grid, A.ord = rep(1,length(si.grid)),
                  z=rep(gam.fit$model$z[1], length(si.grid)))
    }else{
      gam.fit <- gam(y~ A.ord + s(single.index, bs=bs, k=k, sp=sp) +
                       Xm+ s(single.index, by=A.ord, bs=bs, k=k, sp=sp), family=family, method=method)
      dat <- list(Xm = tmp, single.index=si.grid, A.ord = rep(1,length(si.grid)), length(si.grid))
    }
    lp <- predict(gam.fit, newdata = dat, type = "lpmatrix")
    coefs <- coef(gam.fit)
    want <- grep("single.index):", colnames(lp))
    fits <- as.vector(cbind(1.414, lp[,want]) %*% coefs[c(2,want)])
    se.fit <- predict(gam.fit, newdata = dat, se.fit =TRUE)$se.fit
    ci.info <- matrix(0, length(si.grid), 3)
    ci.info[,1] <- fits
    ci.info[,2] <- fits +  1.96*se.fit
    ci.info[,3] <- fits -  1.96*se.fit
    colnames(ci.info) <- c("value", "upper", "lower")
    matplot(si.grid, ci.info[,c("value", "upper", "lower")],
            type="l", lty=c(1,2,2), col=1,
            ylab="Treatment effect contrast", xlab="Index")
    rug(single.index, ticksize = 0.02, col="blue")
    abline(0,0, col = 2, lty =2, lwd=0.8)
  }
  
  
  
  list(beta.coef = round(beta.coef,4),
       beta.ini = beta.ini,
       gamma.magnitude=gamma.magnitude,
       d.beta=d.beta, beta.path=beta.path,
       g.fit=g.fit, gam.fit=gam.fit,
       beta.fit=beta.fit,
       X.scale=X.scale, X.center = X.center,
       y=y, A=A, X=X, Xm = Xm, single.index=single.index,
       p=p, n=n, bs=bs, k=k, L=L, linear.link=linear.link,
       random.effect=random.effect, AIC=AIC, scale.si.01=scale.si.01)
}





#' A subfunction used in estimation
#'
#' This function computes the 1st derivative of the treatment-specific link function with respect to the single index, using finite difference.
#'
#' @param g.fit  a \code{mgcv::gam} object
#' @param eps a small finite difference used in numerical differentiation.
#' @seealso \code{fit.simml}, \code{simml}
#'
## a utility funciton used in fit.simml(); this computes the first deriviative of link functions.
der.link <- function(g.fit, eps=10^(-6)){
  m.terms <- attr(terms(g.fit), "term.labels")
  newD  <- model.frame(g.fit)[, m.terms, drop=FALSE]
  newDF <- data.frame(newD)  # needs to be a data frame for predict
  X0 <- predict.gam(g.fit, newDF, type="lpmatrix")
  newDF[,"single.index"] <- newDF[,"single.index"] + eps
  X1 <- predict.gam(g.fit, newDF, type="lpmatrix")
  Xp <- (X1-X0)/eps
  Xi <- Xp*0
  want <- grep("single.index", colnames(X1))
  Xi[,want] <- Xp[,want]
  g.der  <- as.vector(Xi%*%coef(g.fit))  # the first derivative of the link function
  g.der
}


#' SIMML prediction function
#'
#' This function makes predictions from an estimated SIMML, given a (new) set of pretreatment covariates.
#' The function returns a set of predicted outcomes for each treatment condition and a set of recommended treatment assignments (assuming a larger value of the outcome is better).
#'
#' @param simml.obj  a \code{simml} object
#' @param newX  a (n-by-p) matrix of new values for the covariates X at which predictions are to be made.
#' @param newA  a (n-by-L) matrix of new values for the treatment A at which predictions are to be made.
#' @param newXm a (n-by-q) matrix of new values for the covariates associated with the fitted main effect Xm at which predictions are to be made.
#' @param single.index  a length n vector specifying new values for the single-index at which predictions are to be made; the default is \code{NULL}.
#' @param type the type of prediction required; the default "response" is on the scale of the response variable; the alternative "link" is on the scale of the linear predictors.
#' @param maximize the default is \code{TRUE}, assuming a larger value of the outcome is better; if \code{FALSE}, a smaller value is assumed to be prefered.
#'
#' @return
#' \item{pred.new}{a (n-by-L) matrix of predicted values; each column represents a treatment option.}
#' \item{trt.rule}{a (n-by-1) vector of suggested treatment assignments}
#'
#'
#' @author Park, Petkova, Tarpey, Ogden
#' @seealso \code{simml},\code{fit.simml}
#' @export
#'
pred.simml  <-  function(simml.obj, newX=NULL, newA =NULL, newXm =NULL, single.index=NULL, type = "link", maximize=TRUE)
{
  #if(!inherits(simml.obj, "simml"))   # checks input
  #  stop("obj must be of class `simml'")
  
  #if(ncol(newX) != simml.obj$p)
  #  stop("newX needs to be of p columns ")
  
  #####
  if(is.null(single.index)){
    
    if(is.null(newX)){
      newX  <- simml.obj$X
      if(is.null(newXm)) newXm <- simml.obj$Xm
    }else{
      if(is.null(newXm)){
        if(is.matrix(simml.obj$Xm)){ newXm <- matrix(0, nrow(newX), ncol(simml.obj$Xm)) }
        else{ newXm <- rep(0, nrow(newX)) }
      }
    }
    #if(ncol(newX) != simml.obj$p) stop("newX needs to be of p columns ")
    
    if(is.null(simml.obj$X.scale)){
      if(is.null(simml.obj$X.center)){
        newX.scaled <- scale(newX, center = rep(0,simml.obj$p), scale = rep(1,simml.obj$p))
      }else{
        newX.scaled <- scale(newX, center = simml.obj$X.center, scale = rep(1,simml.obj$p))
      }
    }else{
      newX.scaled <- scale(newX, center = simml.obj$X.center, scale = simml.obj$X.scale)
    }
    if(simml.obj$scale.si.01)  newX.scaled <- cbind(1,newX.scaled)
    single.index  <- newX.scaled %*% simml.obj$beta.coef
    
  }else{
    if(is.null(newXm)) newXm <- matrix(0, length(single.index), ncol(simml.obj$Xm))
  }
  
  
  L <- simml.obj$L  # the number of treatment options (levels)
  A.levels <- unique(simml.obj$A)
  
  # compute treatment-specific predicted outcomes
  if(is.null(newA)){
    pred.new <- matrix(NA, length(single.index), L)
    for(a in 1:L){
      if(simml.obj$random.effect){
        newD <- list(A= rep(A.levels[a], length(single.index)), Xm = as.matrix(newXm),
                     single.index=single.index, z=rep(simml.obj$g.fit$model$z[1], length(single.index)))
      }else{
        newD <- list(A= rep(A.levels[a], length(single.index)), Xm = as.matrix(newXm),
                     single.index=single.index)
      }
      pred.new[ ,a] <- predict.gam(simml.obj$g.fit, newD, type =type)
    }
  }else{
    pred.new <- rep(NA, length(single.index))
    if(simml.obj$random.effect){
      newD <- list(A= newA, Xm = newXm, single.index=single.index, z=rep(simml.obj$g.fit$model$z[1], length(single.index)))
    }else{
      newD <- list(A= newA, Xm = newXm, single.index=single.index)
    }
    pred.new <- predict.gam(simml.obj$g.fit, newD, type=type, newdata.guaranteed=TRUE)
  }
  
  
  # compute optimal treatment assignment
  if(maximize){
    opt.trt.index <- apply(pred.new, 1, which.max)
  }else{
    opt.trt.index <- apply(pred.new, 1, which.min)
  }
  
  trt.rule <- rep(NA, nrow(pred.new))
  for(i in 1:nrow(pred.new)){
    trt.rule[i] <- A.levels[opt.trt.index[i]]
  }
  
  if(L==2)  colnames(pred.new) <- c("Trt1", "Trt2")
  
  return(list(trt.rule = trt.rule, pred.new = pred.new, single.index=single.index))
}



#' A data generation function
#'
#' \code{generate.data} generates an example dataset from a mean model that has a "main" effect component and a treatment-by-covariates interaction effect component (and a random component for noise).
#'
#' @param n  sample size.
#' @param p  dimension of covariates.
#' @param family specifies the distribution of the outcome y;  "gaussian", "binomial", "poisson"; the defult is "gaussian"
#' @param sigma  standard deviation of the random noise term (for gaussian response).
#' @param sigmaX  standard deviation of the covariates.
#' @param correlationX  correlation among the covariates.
#' @param pi.1  probability of being assigned to the treatment 1
#' @param s  controls the nonliarity of the treatment-specific link functions that define the interaction effect component.
#' \describe{
#' \item{\code{s=1}}{linear}
#' \item{\code{s=2}}{nonlinear}
#' }
#' @param delta  controls the intensity of the main effect; can take any intermediate value, e.g., \code{delta= 1.4}.
#' \describe{
#' \item{\code{delta=1}}{moderate main effect}
#' \item{\code{delta=2}}{big main effect}
#' }
#'
#' @param true.beta  a p-by-1 vector of the true single-index coefficients (associated with the interaction effect component); if \code{NULL}, \code{true.beta} is set to be \code{(1, 0.5, 0.25, 0.125, 0,...0)}' (only the first 4 elements are nonzero).
#' @param true.eta   a p-by-1 vector of the true main effect coefficients; if \code{NULL}, \code{true.eta} is set to be \code{(0,..., 0.125, 0.25, 0.25, 1)}' (only the last 4 elements are nonzero).
#'
#'
#' @return
#' \item{y}{a n-by-1 vector of treatment outcomes.}
#' \item{A}{a n-by-1 vector of treatment indicators.}
#' \item{X}{a n-by-p matrix of pretreatment covariates.}
#' \item{SNR}{the "signal" (interaction effect) to "nuisance" (main effect) variance ratio (SNR) in the canonical parameter function.}
#' \item{true.beta}{the true single-index coefficient vector.}
#' \item{true.eta}{the true main effect coefficient vector.}
#' \item{optTr}{a n-by-1 vector of treatments, indicating the optimal treatment selections.}
#' \item{value.opt}{the "value" implied by the optimal treatment decision rule, \code{optTr}.}
#' @export
#'
generate.data <- function(n = 200, # number of observations
                          p = 10,  # number of covariates
                          family = "gaussian",  # the distribution of the outcome y
                          correlationX= 0, # correlation among pretreatment covariates X
                          sigmaX = 1, # pretreatment covariate sd
                          sigma = 0.4, # error sd (for gaussian response)
                          s = 2, # shape of the interaction effect curves (1 linear; 2 nonlinear)
                          delta = 1,  # magnitude of the main effect
                          pi.1 = 0.5,  # probability of being assigned to the treatment 1
                          true.beta= NULL,
                          true.eta= NULL)  # "binomial", type of the outcome variable
{
  
  if(is.null(true.beta)){
    true.beta <- c(c(1, 0.5, 0.25, 0.125), rep(0, p-4))  # only the first 4 components are nonzero.
    true.beta <- true.beta/sqrt(sum(true.beta^2))    # the true single index coefficients
  }
  if(length(true.beta)!= p)   stop("true.beta must be of length p")
  
  if(is.null(true.eta)){
    #eta.hold <- rnorm(4, 0, 1);     # randomly generate the coefficients associated with the main effects
    eta.hold <- c(1,2,3,4)
    eta.hold  <- eta.hold /sqrt(sum(eta.hold^2) )
    true.eta <- c(rep(0, p-4), eta.hold)   # only the last 4 components are nonzero.
  }
  if(length(true.eta)!= p)   stop("true.eta must be of length p")
  
  
  # the link function (the curve that defines the interaction effects);
  # s is the nonlinearity parameter (s=1: linear; s=2: nonlinear)
  g <- function(u, s){
    if(s==1) return(0.3* u)
    if(s==2) return( exp(-(u-0.5)^2) - 0.6)
  }
  
  # delta is the intensity parametr (delta = 1: moderate main effect; delta=2: big main effect)
  m <- function(u, delta= 1)   0.5*delta*sin(u*0.5*pi) #delta*cos(u*0.5*pi)   # this curve defines the main effects
  
  
  # Treatment variable
  A <- drop(rbinom(n, 1, pi.1) + 1)  # generate treatment variables
  
  # Pre-treatment covariates
  Psix <- sigmaX*(diag(1 - correlationX, nrow = p, ncol = p) + matrix(correlationX, nrow = p, ncol = p) )   # X covariance matrix.
  ePsix <- eigen(Psix)
  X <- sapply(1:p, function(x) rnorm(n)) %*% diag(sqrt(ePsix$values)) %*% t(ePsix$vectors)
  
  # X main effect
  main.effect  <-  m(drop(X %*% true.eta), delta)
  
  # A-by-X interaction effect
  TIE <- g(drop(X %*% true.beta), s)
  interaction.effect <- 2*(as.numeric(A) + pi.1 -2) * TIE   #  (-1)^A* TIE
  
  # the hypothetical (potential) outcomes, for each treatment
  if(family == "gaussian"){
    mu.inter1 <-  2*(pi.1 - 1) *TIE  #-TIE        # if A =1
    mu.inter2 <-  2*pi.1*TIE  # TIE        # if A =2
  }
  if(family == "binomial"){
    mu.inter1 <-   1/(1+ exp(-(2*(pi.1 - 1) *TIE)))    # 1/(1+ exp(-(-TIE)))  # if Tr =1
    mu.inter2 <-   1/(1+ exp(- 2*pi.1*TIE ))    # 1/(1+ exp(- TIE ))   # if Tr =2
  }
  if(family == "poisson"){
    mu.inter1 <-  exp(2*(pi.1 - 1)*TIE)  # exp(-TIE)   # if Tr =1
    mu.inter2 <-  exp(2*pi.1*TIE)  # exp( TIE)   # if Tr =2
  }
  
  # the canonical parameter
  theta <- main.effect + interaction.effect
  
  # the "signal" to "noise" ratio
  SNR <- var(interaction.effect)/var(main.effect)
  
  if(family == "gaussian"){
    mu <- theta
    y <-  mu  + sigma * rnorm(n)
  }
  if(family == "binomial"){
    mu <- 1/(1+ exp(-theta))
    y <- rbinom(n, size=1, prob= mu)
  }
  if(family == "poisson"){
    mu <-  exp(theta)
    y <- rpois(n, lambda= mu)
  }
  
  optTr <- as.numeric(mu.inter2 > mu.inter1) + 1  # this takes 1 or 2
  value.opt <- mean(mu[A == optTr ])
  value.opt
  
  return(list(y=y, A =A, X=X, SNR=SNR, true.beta=true.beta, true.eta=true.eta, delta=delta, s=s,
              mu.inter1=mu.inter1, mu.inter2=mu.inter2, optTr=optTr, value.opt=value.opt))
}




#' A function for ordinal categorical response data generation.
#'
#' \code{ordinal.data} generates ordered category response data (with p covariates and a treatment variable).
#'
#' @param n  sample size.
#' @param p  dimension of covariates.
#' @param R  number of response levels in y
#' @param delta  magnitude of "main" effect (i.e., "nuisance" effect) of the covariates; a large delta means a larger "nuisance" variance.
#' @param s  type of the treatment-by-covariates interation effect ("linear" or "nonlinear")
#' @param sigma  noise sd in the latent variable representation
#'
#' @return
#' \item{y}{a n-by-1 vector of treatment outcomes.}
#' \item{A}{a n-by-1 vector of treatment indicators.}
#' \item{X}{a n-by-p matrix of pretreatment covariates.}
#' \item{SNR}{the "signal" (interaction effect) to "nuisance" (main effect) variance ratio (SNR) in the canonical parameter function.}
#' \item{true.beta}{the true single-index coefficient vector.}
#' \item{delta}{magnitude of "main" effect.}
#' \item{s}{type of the treatment-by-covariates interation effect.}
#' @export
ordinal.data <- function(n=400,  # number of subjects
                         p=10,   # number of covariates
                         R = 11,  # number of response levels in y
                         delta = 1,  # magnitude of "main" effect (i.e., "nuisance" effect) of the covariates; a large delta means a larger "nuisance" variance.
                         s = "nonlinear", # type of the treatment-by-covariates interation effect
                         sigma=0)
{
  
  # generate p pretreatment scalar-valued covariates X
  X  = matrix(runif(n*p, -pi/2,pi/2), n, p)  # just use uniform distribution to generate X.
  A  = rbinom(n, 1, 0.5) + 1  # treatment variable taking a value in {1,2} with equal prob.
  
  # X main effect on y
  main.effect = rep(0, n)
  for(j in 1:p){
    main.effect <- main.effect + plogis(2*X[,j])-0.5
  }
  # X-by-A interaction effect on y
  true.beta <- c(c(1,-1, 0.5, -0.5), rep(0, p-4))  # the "true" single-index coefficient assocaited with X-by-A interaction effect.
  if(s=="nonlinear"){
    contrast <- 2*exp(-(X%*%true.beta-0.5)^2)-0.5   # this will specify a nonlinear interaction term
  }else{
    contrast <-  X %*% true.beta   # this will specify linear interaction term.
  }
  interaction.effect <- (-1)^A *contrast   # X-ty-A interaction effect term.
  
  f <- delta*main.effect + interaction.effect   + rnorm(n, 0, sd=sigma)  # latent response.
  f <- f - mean(f)  # center the latent response.
  
  var.main =  var(delta*main.effect)
  var.interaction = var(interaction.effect)
  SNR = var.interaction/ var.main  # "interaction effect"-to-"main effect" variance ratio
  
  noise <- rlogis(n, 0, scale = 1)  # standard logistic noise
  y.star <- f + noise # latent response
  theta  <-c(-Inf,seq(-3,3, length.out=R-1), Inf)  # the cut-points
  y <- f
  for(j in 1:R){
    y[y.star>theta[j] & y.star<=theta[j+1]] <- j
  }
  
  list(y = y, A = A, X = X, SNR = SNR, true.beta = true.beta,
       delta = delta, s = s)
}



######################################################################
## END OF THE FILE
######################################################################