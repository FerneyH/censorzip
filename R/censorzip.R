#' Censored Zero-inflated Poisson Regression Model
#' @param formula  symbolic description of the model, see \code{details}.
#' @param data arguments controlling formula.
#' @param status censor indicator variable, the status indicator, normally 0 = censored.
#' Other choices are TRUE/FALSE (FALSE = censored). 
#' @param dist character specification of count model family (a log link is always used).
#' @param control a list of control arguments specified via censorZI.control.
#' @param model,y,x  logicals. If TRUE the corresponding components of the 
#' fit (model frame, response, model matrix) are returned.
#' @param ... arguments passed to zeroinfl.control in the default setup.
#' @import stats
#' @importFrom Rdpack reprompt
#' @importFrom maxLik maxLik
#' @importFrom pscl zeroinfl
#' @return An object of class "censorzip", i.e., a list with components including
#'  
#' 
#' \code{coef} a vector of fitted means,
#' 
#' \code{vcov} covariance matrix of all coefficients in the model (derived from the Hessian),
#' 
#' \code{loglik} log-likelihood of the fitted model,
#' 
#' \code{dist} character string describing the count distribution used,
#' 
#' \code{call} the original function call,
#' 
#' \code{formula} the original formula,
#' 
#' \code{model} the full model frame (if model = TRUE),
#' 
#' \code{y} the response count vector (if y = TRUE),
#' 
#' \code{x} a list with elements "count" and "zero" containing the model matrices from the respective models (if x = TRUE),
#' 
#' @author Ferney Henao-Ceballos
#' 
#' @description
#'    We utilize the simplification of GLMs to derive a new algorithm for maximum
#'    likelihood estimation by the Newton-Rhapson method in the Censored 
#'    Zero-inflated Poisson Regression Model.
#'    
#' @details
#'    Censored Zero-inflated Poisson model is two-component mixture model combining 
#'    a point mass at zero with a Censored Poisson distribution. Thus, there are
#'    two sources of zeros: zeros may come from both the point mass and from the 
#'    count component. For modeling the unobserved state (zero vs. count), a binary
#'    model is used that captures the probability of zero inflation. In the simplest 
#'    case only with an intercept but potentially containing regressors. 
#'    
#'    The formula can be used to specify both components of the model: If a formula 
#'    of type y ~ x1 + x2 is supplied, then the same regressors are employed in both
#'    components. This is equivalent to y ~ x1 + x2 | x1 + x2. Of course, a different
#'    set of regressors could be specified for the count and zero-inflation component,
#'    e.g., y ~ x1 + x2 | z1 + z2 + z3 giving the count data model y ~ x1 + x2
#'    conditional on (|) the zero-inflation model y ~ z1 + z2 + z3. A simple inflation 
#'    model where all zero counts have the same probability of belonging to the zero 
#'    component can by specified by the formula y ~ x1 + x2 | 1.
#'    
#'
#' @references 
#' \insertRef{Rpack:bibtex}{Rdpack}
#' 
#' \insertRef{maxlik}{censorzip}
#' 
#' \insertRef{pscl}{censorzip}
#' 
#' \insertRef{nguyen2021asymptotic}{censorzip}
#' 
#' \insertRef{chen2020censored}{censorzip}
#' 
#' @export

censorzip<-function(formula, data, status, dist = "poisson",
                    control = censorzip_control(...),
                    model = TRUE, y = TRUE, x = FALSE, ...){
  # start value
  start<-as.vector(coef(zeroinfl(formula,data)))
  
  ## call and formula
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  
  ## extended formula processing
  if(length(formula[[3]]) > 1 && identical(formula[[3]][[1]], as.name("|")))
  {
    ff <- formula
    formula[[3]][1] <- call("+")
    mf$formula <- formula
    ffc <- . ~ .
    ffz <- ~ .
    ffc[[2]] <- ff[[2]]
    ffc[[3]] <- ff[[3]][[2]]
    ffz[[3]] <- ff[[3]][[3]]
    ffz[[2]] <- NULL
  } else {
    ffz <- ffc <- ff <- formula
    ffz[[2]] <- NULL
  }
  
  ## call model.frame()
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  ## extract terms, model matrices, response
  mt <- attr(mf, "terms")
  mtX <- terms(ffc, data = data)
  X <- model.matrix(mtX, mf)
  mtZ <- terms(ffz, data = data)
  mtZ <- terms(update(mtZ, ~ .), data = data)
  Z <- model.matrix(mtZ, mf)
  Y <- model.response(mf, "numeric")
  
  
  
  ## sanity checks
  if(length(Y) < 1) stop("empty model")
  if(all(Y > 0)) stop("invalid dependent variable, minimum count is not zero")  
  if(!isTRUE(all.equal(as.vector(Y), as.integer(round(Y + 0.001)))))
    stop("invalid dependent variable, non-integer values")
  Y <- as.integer(round(Y + 0.001))
  if(any(Y < 0)) stop("invalid dependent variable, negative counts")
  
  
  ## convenience variables
  P<-ncol(X);s<-ncol(Z);n<-nrow(X)
  I<-as.numeric(Y==0) # Indicator y_i=0
  di<-status # censored data, status == 0 (censored)
  C<-Y*(1-status) # Ci
  
  
  ##Gradient function  
  gradPois<-function(beta){
    mu<-exp(X%*%beta[1:P]) # mui
    w<- exp(Z%*%beta[(P+1):(P+s)]) #  wi
    b<-1/(w+1)
    t<-1/(w+dpois(0,mu))
    pi<-dpois(C,mu)/(1-ppois(C-1,mu)) # phi 
    
    U1<- t(Z)%*%(di*w*t*I - w*b)
    U2<- t(X)%*%(-di*dpois(0,mu)*mu*t*I+di*(Y-mu)*(1-I) + (1-di)*C*pi)
    as.vector(rbind(U2,U1))} # Derivative vector
  
  
  ## Hessian
  HessPois<-function(beta){
    mu<-exp(X%*%beta[1:P]) # mui
    w<- exp(Z%*%beta[(P+1):(P+s)]) # wi
    b<-1/(w+1)
    t<-1/(w+dpois(0,mu))
    pi<-dpois(C,mu)/(1-ppois(C-1,mu)) # phi
    R<-diag(as.vector(di*dpois(0,mu)*w*t^2*I-w*b^2))
    
    J11<-t(Z)%*%R%*%Z # Second derivative for gamma
    S<-diag(as.vector(-di*(dpois(0,mu)*mu*(w+dpois(0,mu)-w*mu)*t^2*I+mu*(1-I))+
                        (1-di)*C*((C-mu)*pi-C*pi^2)))
    J22<-t(X)%*%S%*%X # Second derivative for beta
    K<-diag(as.vector(di*dpois(0,mu)*mu*w*t^2*I))
    J12<-t(Z)%*%K%*%X # Cross derivative
    as.matrix(rbind(cbind(J22,t(J12)),cbind(J12,J11)))} # Hessian
  
  
  # log-likelihood Function
  likelihood_zipoiss<-function(beta){
    sum(di*(log(exp(Z%*%beta[(P+1):(P+s)])+dpois(0,exp(X%*%beta[1:P])))*I+log(dpois(Y,exp(X%*%beta[1:P])))*(1-I))
        +(1-di)*log(1-ppois(C-1,exp(X%*%beta[1:P])))-log(1+exp(Z%*%beta[(P+1):(P+s)])))}
  
  
  
  
  loglikfun <- switch(dist, "poisson" = likelihood_zipoiss)
  gradfun <- switch(dist, "poisson" = gradPois)
  hessianfun<-switch(dist, "poisson" = HessPois)
  
  
  censorzip_control <- function(tol=-1, ...) {
    rval <- list(tol = tol)
    rval <- c(rval, list(...))
    if(is.null(rval$reltol)) rval$reltol <- .Machine$double.eps
    if(is.null(rval$gradtol)) rval$gradtol <- .Machine$double.eps
    rval
  }
  
  
  fit <- maxLik(logLik = loglikfun, grad = gradfun,
                start = start, hess = hessianfun, control=censorzip_control())
  
 
   ## coefficients and covariances
  coefc <- fit$estimate[1:P]
  names(coefc) <- colnames(X)
  coefz <- fit$estimate[(P+1):(P+s)]
  names(coefz) <- colnames(Z)
  
  vc<-vcov(fit)
  colnames(vc) <- rownames(vc) <- c(paste("count", colnames(X), sep = "_"),
                                    paste("zero",  colnames(Z), sep = "_"))
  
  
  rval<-list(coefficients =list(count = coefc, zero = coefz),
         niter = fit$iterations,
         message = fit$message,
         converged = fit$code,
         gradient = fit$gradient,
         hessian = fit$hessian,
         loglik= logLik(fit),
         AIC= AIC(fit),
         vcov = vc,
         dist = dist,
         formula = ff,
         call = cl)
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- list(count = X, zero = Z)
  
  class(rval) <- "censorzip"
  return(rval)
}

#' @method print censorzip
#' @export


print.censorzip <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  
    cat(paste("Count model coefficients (", x$dist, " with log link):\n", sep = ""))
    print.default(format(x$coefficients$count, digits = digits), print.gap = 2, quote = FALSE)
    
    cat(paste("\nZero-inflation model coefficients (binomial with logit link):\n", sep = ""))
    print.default(format(x$coefficients$zero, digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
  
  invisible(x)
}


#' @method coef censorzip
#' @export


coef.censorzip <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  rval <- object$coefficients
  rval <- switch(model,
                 "full" = structure(c(rval$count, rval$zero),
                                    .Names = c(paste("count", names(rval$count), sep = "_"),
                                               paste("zero", names(rval$zero), sep = "_"))),
                 "count" = rval$count,
                 "zero" = rval$zero)
  rval
}

#' @method vcov censorzip
#' @export

vcov.censorzip <- function(object, ...) {
return(object$vcov)
}


#' @method AIC censorzip
#' @export


AIC.censorzip <- function(object, ...) {
  structure(object$AIC, class = "AIC")
}


#' @method logLik censorzip
#' @export


logLik.censorzip <- function(object, ...) {
structure(object$logLik, class = "maxLik")
}

#' @method summary censorzip
#' @export

summary.censorzip <- function(object,...)
{
  ## compute z statistics
  kc <- length(object$coefficients$count)
  kz <- length(object$coefficients$zero)
  se <- sqrt(diag(object$vcov))
  coef <- c(object$coefficients$count, object$coefficients$zero)  
  zstat <- coef/se
  pval <- 2*pnorm(-abs(zstat))
  coef <- cbind(coef, se, zstat, pval)
  colnames(coef) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  object$coefficients$count <- coef[1:kc,,drop = FALSE]
  object$coefficients$zero <- coef[(kc+1):(kc+kz),,drop = FALSE]
  
  ## delete some slots
  object$gradient <- object$model <- object$y <- object$x <- object$hessian  <- NULL
  
  ## return
  class(object) <- "summary.censorzip"
  object
}

#' @method print summary.censorzip
#' @export


print.summary.censorzip<- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  
    cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
    
    cat(paste("\nCount model coefficients (", x$dist, " with log link):\n", sep = ""))
    printCoefmat(x$coefficients$count, digits = digits, signif.legend = FALSE)
    
    cat(paste("\nZero-inflation model coefficients (binomial with logit link):\n", sep = ""))
    printCoefmat(x$coefficients$zero, digits = digits, signif.legend = FALSE)
    
    if(getOption("show.signif.stars") & any(rbind(x$coefficients$count, x$coefficients$zero)[,4] < 0.1, na.rm=TRUE))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
    
    cat(paste("Number of iterations in Newton-Raphson optimization:", x$niter, "\n"))
    cat(paste("Log-likelihood:", round(x$loglik[1],3), "\n"))
}



