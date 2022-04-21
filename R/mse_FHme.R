#' @title Mean Squared Error Estimator of the EBLUP under a Fay-Herriot Model with Measurement Error
#' @description This function gives the mean squared error estimator of the EBLUP based on Fay-Herriot model with measurement error using jackknife method.
#' @param formula an object of class \code{\link[stats]{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The variables included  \code{formula} must have a length equal to the number of domains \code{m}. This formula can provide auxiliary variable either measured with error or without error or combination between them. If the auxiliary variable are combination between \code{noerror} and \code{witherror} variable, input all \code{witherror} variable first then \code{noerror} variable.
#' @param vardir vector containing the \code{m} sampling variances of direct estimators for each domain. The values must be sorted as the \code{Y}.
#' @param var.x vector containing mean squared error of \code{X} . The values must be sorted as the \code{X}. if you use optional \code{data}, input this parameter use \code{c("")}, example: \code{var.x = c("c1") or var.x = c("c1","c2")}.
#' @param type.x type of auxiliary variable used in the model. Either source measured with \code{noerror}, \code{witherror} and \code{mix}. Default value is \code{witherror}.
#' @param MAXITER maximum number of iterations allowed. Default value is \code{1000} iterations.
#' @param PRECISION convergence tolerance limit. Default value is \code{0.0001}.
#' @param data optional data frame containing the variables named in formula, vardir, and var.x.
#' @details A formula has an implied intercept term. To remove this use either y ~ x - 1 or y ~ 0 + x. See \code{\link[stats]{formula}}  for more details of allowed formulae.
#'
#' @return The function returns a list with the following objects:
#' \describe{
#'    \item{\code{mse}}{vector with the values of the mean squared errors of the EBLUPs for each domain.}
#'  }
#'
#' @examples
#' data(dataME)
#' data(datamix)
#' \donttest{
#' mse.sae.me <- mse_FHme(formula = y ~ x.hat, vardir = vardir, var.x = c("var.x"), data = dataME)
#' mse.sae.mix <- mse_FHme(formula = y ~ x.hat1 + x.hat2 + x3 + x4,
#'                 vardir = vardir, var.x = c("var.x1", "var.x2"), type.x = "mix", data = datamix)
#' }
#' @export mse_FHme
mse_FHme <- function(formula, vardir, var.x, type.x = "witherror", MAXITER = 1000, PRECISION = 0.0001, data) {
  namevar <- deparse(substitute(vardir))
  #name_c <- deparse(substitute(c))
  if (type.x != "witherror" & type.x != "noerror" & type.x != "mix")
    stop(" type.x=\"", type.x, "\" must be \"witherror\", \"noerror\" or \"mix\".")
  if(!missing(data)){
    formuladata <- model.frame(formula, na.action = na.omit, data)
    y <- formuladata[, 1]
    X_cap <- model.matrix(formula, data)
    c_dim <- dim(X_cap)[2]
    psi <- data[, namevar]
    if (type.x == "witherror") {
      c <- data[, var.x]
      est <- FHme(formula = y ~ X_cap - 1, vardir = psi, var.x = c, type.x = type.x,
                  MAXITER = MAXITER, PRECISION = PRECISION)
    } else if(type.x == "noerror"){
      c <- matrix(0, nrow = dim(X_cap)[1], ncol = c_dim - 1)
      est <- FHme(formula = y ~ X_cap - 1, vardir = psi, var.x = c, type.x = type.x,
                  MAXITER = MAXITER, PRECISION = PRECISION)
    } else{
      c_left <- data[, var.x]
      est <- FHme(formula = y ~ X_cap - 1, vardir = psi, var.x = c_left, type.x = type.x,
                  MAXITER = MAXITER, PRECISION = PRECISION)
      c_left.tmp <- data.frame(c_left)
      c_right <- matrix(0, nrow = dim(X_cap)[1], ncol = (c_dim - 1) - dim(c_left.tmp)[2])
      c <- cbind(c_left, c_right)
    }

  } else{
    formuladata <- model.frame(formula, na.action = na.omit)
    y <- formuladata[, 1]
    X_cap <- model.matrix(formula)
    psi <- vardir
    est <- FHme(formula = y ~ X_cap - 1, vardir = psi, var.x = var.x, type.x = type.x,
                MAXITER = MAXITER, PRECISION = PRECISION)
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

  c <- cbind(0, c)
  m <- length(y)
  p <- dim(X_cap)[2]

  #jackknife function, resampling by excluding area j
  jackknife <- function(y,X_cap,psi,c,j) {
    w = rep(1,length(y))
    m <- length(y)
    p <- dim(X_cap)[2]
    diff_beta <- as.matrix(rep(1,p))
    diff_sigma <- 1
    R_sigma <- 0.001
    R_beta <- as.matrix(rep(0.001,p))
    max_iter <- 100
    betacap <- 0
    sigma2cap <- 0
    k <- 0

    while((any(diff_beta > R_beta) | (diff_sigma > R_sigma)) & (k < max_iter)){
      betacap_tmp <- betacap
      sigma2cap_tmp <- sigma2cap
      yj <- y[-j]
      Xj <- X_cap[-j,]
      Cj <- c[-j,]
      psij <- psi[-j]
      mj <- length(yj)
      pj <- dim(Xj)[2]

      if(all(diff_beta < R_beta)){
        mpj <- 1/(mj-pj)
        point <- sum(sapply(1:mj, function(i){
          X_cap_i <- as.matrix(Xj[i,])
          point <- (yj[i] - t(X_cap_i)%*%betacap)^2 - psij[i] - t(betacap)%*%diag(Cj[i,],nrow = pj)%*%betacap
          return(point)}))
        sigma2cap <- mpj*point

        #Ybarra & Lohr 2008: if estimated random effect < 0 set to zero
        if(sigma2cap < 0) {
          sigma2cap <- 0
        }

        w <- sapply(1:m, function(i){
          wi <- 1/(sigma2cap + psi[i] + t(betacap)%*%diag(c[i,],nrow = p)%*%betacap)
          return(wi)
        })
        diff_sigma <- sigma2cap - sigma2cap_tmp

      } else if(diff_sigma < R_sigma){
        wj <- w[-j]
        wC <- Reduce('+', lapply(1:mj, function(i){
          C <- diag(Cj[i,],nrow = pj)
          wC_i <- w[i]*C
          return(wC_i)}))
        Xtwi <- t(wj * Xj)
        chloe <- Xtwi %*% Xj - wC
        QR.fac <- qr(chloe)

        if (QR.fac$rank == pj) {
          Q_matrix <- solve(chloe)
          betacap <- Q_matrix %*% Xtwi %*% yj

        }else{
          G <- Xtwi %*% Xj
          Gsqrt <- qr.R(qr(sqrt(wj) * Xj))
          if(sum(eigen(G)$value > 0) == pj){
            invG <- backsolve(Gsqrt, diag(pj))
          } else {
            invG <- ginv(Gsqrt)
          }
          eigenGwCG <- eigen(invG%*%wC%*%invG)
          pOrth <- eigenGwCG$vector
          lDiag <- diag(eigenGwCG$values, nrow = pj)
          D <- diag(sapply(1:pj, function(j)
          {
            Djj <- 0
            if (1-lDiag[j,j] > 1/m)
              Djj <- 1/(1-lDiag[j,j])
            return(Djj)
          }), nrow = p)
          Q_matrix <- invG %*% crossprod(sqrt(D), pOrth) %*% invG
          betacap <- Q_matrix %*% Xtwi %*% yj
        }

        w <- sapply(1:m, function(i){
          wi <- 1/(sigma2cap + psi[i] + t(betacap)%*%diag(c[i,],nrow = p)%*%betacap)
          return(wi)
        })
        diff_beta <- betacap - betacap_tmp

      } else {
        wj <- w[-j]
        wC <- Reduce('+', lapply(1:mj, function(i){
          C <- diag(Cj[i,],nrow = pj)
          wC_i <- w[i]*C
          return(wC_i)}))
        Xtwi <- t(wj * Xj)
        chloe <- Xtwi %*% Xj - wC
        QR.fac <- qr(chloe)

        if (QR.fac$rank == pj) {
          Q_matrix <- solve(chloe)
          betacap <- Q_matrix %*% Xtwi %*% yj
        }else{
          G <- Xtwi %*% Xj
          Gsqrt <- qr.R(qr(sqrt(wj) * Xj))
          if(sum(eigen(G)$value > 0) == pj){
            invG <- backsolve(Gsqrt, diag(pj))
          } else {
            invG <- ginv(Gsqrt)
          }
          eigenGwCG <- eigen(invG%*%wC%*%invG)
          pOrth <- eigenGwCG$vector
          lDiag <- diag(eigenGwCG$values, nrow = pj)
          D <- diag(sapply(1:pj, function(j)
          {
            Djj <- 0
            if (1-lDiag[j,j] > 1/m)
              Djj <- 1/(1-lDiag[j,j])
            return(Djj)
          }), nrow = p)
          Q_matrix <- invG %*% crossprod(sqrt(D), pOrth) %*% invG
          betacap <- Q_matrix %*% Xtwi %*% yj
        }
        mpj <- 1/(mj-pj)
        point <- sum(sapply(1:mj, function(i){
          X_cap_i <- as.matrix(Xj[i,])
          point <- (yj[i] - t(X_cap_i)%*%betacap)^2 - psij[i] - t(betacap)%*%diag(Cj[i,],nrow = pj)%*%betacap
          return(point)}))
        sigma2cap <- mpj*point

        #Ybarra & Lohr 2008: if estimated random effect < 0 set to zero
        if(sigma2cap < 0) {
          sigma2cap <- 0
        }
        w <- sapply(1:m, function(i){
          wi <- 1/(sigma2cap + psi[i] + t(betacap)%*%diag(c[i,],nrow = p)%*%betacap)
          return(wi)
        })
        diff_beta <- betacap - betacap_tmp
        diff_sigma <- sigma2cap - sigma2cap_tmp
      }
      k <- k+1
    }
    gammacap <- sapply(1:m, function(i){
      mse_ri <- sigma2cap + t(betacap)%*%diag(c[i,],nrow = p)%*%betacap
      return(mse_ri/(mse_ri + psi[i]))})

    X.beta <- X_cap %*% betacap
    yme <- gammacap*y + (1-gammacap)*X.beta

    result <- list('y_me' = yme,
                   'gamma' = gammacap)
    return(result)
  }
  jk <- lapply(1:m, function(j) jackknife(y,X_cap,psi,c,j))

  #Estimate MSE use jacknife method (M1 + M2)
  m1cap <- sapply(1:m, function(i){
    left <- est$fit$gamma[i]*psi[i]
    right <- ((m-1)/m) * sum(sapply(1:m, function(j){
      jkgamma <- jk[[j]]$gamma[i]
      return(left-jkgamma*psi[i])}))
    return(left+right)
  })
  m2cap <- sapply(1:m, function(i){
    m2 <- ((m-1)/m) * sum(sapply(1:m, function(j){
      return((jk[[j]]$y_me[i] - est$eblup[i])^2)}))
    return(m2)
  })
  mse <- m1cap + m2cap
  return(list("mse" = mse))
}


