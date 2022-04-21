#' @title Fay-Herriot Model with Measurement Error
#' @description This function gives the EBLUP estimator based on Fay-Herriot model with measurement error.
#' @param formula an object of class \code{\link[stats]{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The variables included \code{formula} must have a length equal to the number of domains \code{m}. This formula can provide auxiliary variable either measured with error or without error or combination between them. If the auxiliary variable are combination between \code{noerror} and \code{witherror} variable, input all \code{witherror} variable first then \code{noerror} variable.
#' @param vardir vector containing the \code{m} sampling variances of direct estimators for each domain. The values must be sorted as the \code{Y}.
#' @param var.x vector containing mean squared error of \code{X} . The values must be sorted as the \code{X}. if you use optional \code{data}, input this parameter use \code{c("")}, example: \code{var.x = c("c1") or var.x = c("c1","c2")}.
#' @param MAXITER maximum number of iterations allowed. Default value is \code{1000} iterations.
#' @param PRECISION convergence tolerance limit. Default value is \code{0.0001}.
#' @param type.x type of auxiliary variable used in the model. Either source measured with \code{noerror}, \code{witherror} and \code{mix}. Default value is \code{witherror}.
#' @param data optional data frame containing the variables named in formula, vardir, and var.x.
#' @details A formula has an implied intercept term. To remove this use either y ~ x - 1 or y ~ 0 + x. See \code{\link[stats]{formula}}  for more details of allowed formulae.
#'
#' @return The function returns a list with the following objects:
#' \describe{
#'    \item{\code{eblup}}{vector with the values of the estimators for the domains.}
#'    \item{\code{fit}}{a list containing the following objects:}
#'    \itemize{
#'     \item \code{method} : type of fitting method.
#'     \item \code{convergence} : a logical value of convergence when calculating estimated beta and estimated random effects.
#'     \item \code{iterations} : number of iterations when calculating estimated beta and estimated random effects.
#'     \item \code{estcoef} : a data frame with the estimated model coefficient (\code{beta}) in the first column, their standard error (\code{std.error}) in the second column, the t-statistics (\code{t.statistics}) in the third column, and the p-values of the significance of each coefficient (\code{pvalue}) in the last column.
#'     \item \code{refvar} : a value of estimated random effects.
#'     \item \code{gamma} : vector with values of the estimated gamma for each domains.
#'     }
#'  }
#' @seealso \code{\link{mse_FHme}}
#' @examples
#' data(dataME)
#' data(datamix)
#' sae.me <- FHme(formula = y ~ x.hat, vardir = vardir, var.x = c("var.x"), data = dataME)
#' sae.mix <- FHme(formula = y ~ x.hat1 + x.hat2 + x3 + x4,
#'             vardir = vardir, var.x = c("var.x1", "var.x2"), type.x = "mix", data = datamix)
#'
#' @export FHme
FHme <- function(formula, vardir, var.x, type.x = "witherror", MAXITER = 1000, PRECISION = 0.0001, data) {
  namevar <- deparse(substitute(vardir))
  if (type.x != "witherror" & type.x != "noerror" & type.x != "mix")
    stop(" type.x=\"", type.x, "\" must be \"witherror\", \"noerror\" or \"mix\".")

  if(!missing(data)){
    formuladata <- model.frame(formula, na.action = na.omit, data)
    X_cap <- model.matrix(formula, data)
    c_dim <- dim(X_cap)[2]
    if (type.x == "witherror") {
      c <- data[, var.x]
    } else if(type.x == "noerror"){
      c <- matrix(0, nrow = dim(X_cap)[1], ncol = c_dim - 1)
    } else{
      c_left <- data[, var.x]
      c_left.tmp <- data.frame(c_left)
      c_right <- matrix(0, nrow = dim(X_cap)[1], ncol = (c_dim - 1) - dim(c_left.tmp)[2])
      c <- cbind(c_left, c_right)
    }
    psi <- data[, namevar]
  } else{
    formuladata <- model.frame(formula, na.action = na.omit)
    X_cap <- model.matrix(formula)
    psi <- vardir
    if(type.x == "witherror"){
      c <- as.matrix(var.x)
    } else if (type.x == "noerror") {
      c <- matrix(0, nrow = dim(X_cap)[1], ncol = dim(X_cap)[2] - 1)
    } else {
      c_left <- as.matrix(var.x)
      c_right <- matrix(0, nrow = dim(X_cap)[1], ncol = (dim(X_cap)[2] - 1) - dim(c_left)[2])
      c <- cbind(c_left, c_right)
    }

  }
  y <- formuladata[, 1]
  c <- cbind(0, c)
  m <- length(y)
  p <- dim(X_cap)[2]

  if (length(na.action(y)) > 0){
    stop("your dependent variable (y) contains NA values.")
  }
  if (length(na.action(X_cap))){
    stop("your independent variable (X_cap) contains NA values.")
  }
  if (length(na.action(psi))){
    stop("your vardir contains NA values.")
  }
  if (length(na.action(c))){
    stop("your mse.x contains NA values.")
  }
  w = rep(1,length(y))
  m <- length(y)
  p <- dim(X_cap)[2]
  sigma2cap <- 0
  betacap <- 0
  R_sigma <- PRECISION
  R_beta <- as.matrix(rep(PRECISION,p))
  max_iter <- MAXITER
  k <- 0
  diff_beta <- as.matrix(rep(1,p))
  diff_sigma <- 1
  convergence <- TRUE

  #iteration for estimated beta and sigma
  while(((any(diff_beta > R_beta)) | (diff_sigma > R_sigma)) & (k < max_iter)){
      betacap_tmp <- betacap
      sigma2cap_tmp <- sigma2cap
      if(all(diff_beta < R_beta)){ #when beta still didnt converge, but sigma was done
        mp <- 1/(m-p)
        point <- sum(sapply(1:m, function(i){
          X_cap_i <- as.matrix(X_cap[i,])
          point <- (y[i] - t(X_cap_i)%*%betacap)^2 - psi[i] - t(betacap)%*%diag(c[i,],nrow = p)%*%betacap
          return(point)}))
        sigma2cap <- mp*point

        #Ybarra & Lohr 2008: if estimated random effect < 0 set to zero
        if(sigma2cap < 0) {
          sigma2cap <- 0}
        w <- sapply(1:m, function(i){
          wi <- 1/(sigma2cap + psi[i] + t(betacap)%*%diag(c[i,],nrow = p)%*%betacap)
          return(wi)
        })
        diff_sigma <- sigma2cap - sigma2cap_tmp

      } else if(diff_sigma < R_sigma){ #when sigma still didnt converge, but beta was done
        wC <- Reduce('+', lapply(1:m, function(i){
          C <- diag(c[i,],nrow = p)
          wC_i <- w[i]*C
          return(wC_i)}))
        Xtwi <- t(w * X_cap)
        chloe <- Xtwi %*% X_cap - wC
        QR.fac <- qr(chloe)

        if (QR.fac$rank == p) {
          Q_matrix <- solve(chloe)
          betacap <- Q_matrix %*% Xtwi %*% y

        }else{
          G <- Xtwi %*% X_cap
          Gsqrt <- qr.R(qr(sqrt(w) * X_cap))
          if(sum(eigen(G)$value > 0) == p){
            invG <- backsolve(Gsqrt, diag(p))
          } else {
            invG <- ginv(Gsqrt)
          }
          eigenGwCG <- eigen(invG%*%wC%*%invG)
          pOrth <- eigenGwCG$vector
          lDiag <- diag(eigenGwCG$values, nrow = p)
          D <- diag(sapply(1:p, function(j)
          {
            Djj <- 0
            if (1-lDiag[j,j] > 1/m)
              Djj <- 1/(1-lDiag[j,j])
            return(Djj)
          }), nrow = p)
          Q_matrix <- invG %*% crossprod(sqrt(D), pOrth) %*% invG
          betacap <- Q_matrix %*% Xtwi %*% y
        }

        w <- sapply(1:m, function(i){
          wi <- 1/(sigma2cap + psi[i] + t(betacap)%*%diag(c[i,],nrow = p)%*%betacap)
          return(wi)
        })
        diff_beta <- betacap - betacap_tmp

      } else { #when beta and sigma didnt converge
        wC <- Reduce('+', lapply(1:m, function(i){
          C <- diag(c[i,],nrow = p)
          wC_i <- w[i]*C
          return(wC_i)}))
        Xtwi <- t(w * X_cap)
        chloe <- Xtwi %*% X_cap - wC
        QR.fac <- qr(chloe)

        if (QR.fac$rank == p) {
          Q_matrix <- solve(chloe)
          betacap <- Q_matrix %*% Xtwi %*% y

        }else{
          G <- Xtwi %*% X_cap
          Gsqrt <- qr.R(qr(sqrt(w) * X_cap))
          if(sum(eigen(G)$value > 0) == p){
            invG <- backsolve(Gsqrt, diag(p))
          } else {
            invG <- ginv(Gsqrt)
          }
          eigenGwCG <- eigen(invG%*%wC%*%invG)
          pOrth <- eigenGwCG$vector
          lDiag <- diag(eigenGwCG$values, nrow = p)
          D <- diag(sapply(1:p, function(j)
          {
            Djj <- 0
            if (1-lDiag[j,j] > 1/m)
              Djj <- 1/(1-lDiag[j,j])
            return(Djj)
          }), nrow = p)
          Q_matrix <- invG %*% crossprod(sqrt(D), pOrth) %*% invG
          betacap <- Q_matrix %*% Xtwi %*% y
        }
        mp <- 1/(m-p)
        point <- sum(sapply(1:m, function(i){
          X_cap_i <- as.matrix(X_cap[i,])
          point <- (y[i] - t(X_cap_i)%*%betacap)^2 - psi[i] - t(betacap)%*%diag(c[i,],nrow = p)%*%betacap
          return(point)}))
        sigma2cap <- mp*point

        #Ybarra & Lohr 2008: if estimated random effect < 0 set to zero
        if(sigma2cap < 0) {
          sigma2cap <- 0}
        w <- sapply(1:m, function(i){
          wi <- 1/(sigma2cap + psi[i] + t(betacap)%*%diag(c[i,],nrow = p)%*%betacap)
          return(wi)
        })
        diff_beta <- betacap - betacap_tmp
        diff_sigma <- sigma2cap - sigma2cap_tmp
      }
      k <- k+1
    }
  if (k >= max_iter & ((any(diff_beta > R_beta)) | (diff_sigma > R_sigma))) {
    convergence <- FALSE}

  gammacap <- sapply(1:m, function(i){
    mse_ri <- sigma2cap + t(betacap)%*%diag(c[i,],nrow = p)%*%betacap
    return(mse_ri/(mse_ri + psi[i]))})

  X.beta <- X_cap %*% betacap
  yme <- gammacap*y + (1-gammacap)*X.beta

  se.b <- sqrt(diag(Q_matrix))
  t.val <- betacap/se.b
  pv <- 2 * pnorm(abs(t.val), lower.tail = FALSE)
  coef <- data.frame(betacap, se.b, t.val, pv)
  colnames(coef) <- c("beta", "std.error", "t.statistics", "p.value")
  result <- list(eblup = NA, fit = list(method = NA, convergence = NA,
                                        iterations = NA,
                                        estcoef =NA,
                                        refvar = NA,
                                        gamma = NA))
  result$eblup <- yme
  result$fit$method <- 'REML'
  result$fit$convergence <- convergence
  if (result$fit$convergence == FALSE) {
    warning("The fitting method does not converge.\n")
  }
  result$fit$iterations <- k
  result$fit$estcoef <- coef
  result$fit$refvar <- sigma2cap
  result$fit$gamma <- gammacap
  return(result)
}







