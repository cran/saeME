#' @title Fay-Herriot Model with Measurement Error of Nonsampled Area
#' @description This function gives the EBLUP estimator of nonsampled area using cluster information.
#' @param formula an object of class \code{\link[stats]{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The variables included  \code{formula} must have a length equal to the number of domains \code{m}. This formula can provide auxiliary variable either measured with error or without error or combination between them. If the auxiliary variable are combination between \code{noerror} and \code{witherror} variable, input all \code{witherror} variable first then \code{noerror} variable.
#' @param vardir vector containing the \code{m} sampling variances of direct estimators for each domain. The values must be sorted as the \code{Y}.
#' @param var.x vector containing mean squared error of \code{X} . The values must be sorted as the \code{X}. if you use optional \code{data}, input this parameter use \code{c("")}, example: \code{var.x = c("c1") or var.x = c("c1","c2")}.
#' @param type.x type of auxiliary variable used in the model. Either source measured with \code{noerror}, \code{witherror} and \code{mix}. Default value is \code{witherror}.
#' @param n.cluster either the number of clusters, say \code{k}, or a set of initial (distinct) cluster centers.
#' @param data data frame containing the variables named in formula, vardir, and var.x.
#' @details A formula has an implied intercept term. To remove this use either y ~ x - 1 or y ~ 0 + x. See \code{\link[stats]{formula}}  for more details of allowed formulae.
#'
#' @return The function returns a list with the following objects:
#' \describe{
#'    \item{sampled_data}{data frame of nonsampled area containing cluster information and mean of random effect for each cluster.}
#'    \item{nonsampled_data}{data frame of sampled area containing cluster information and mean of random effect for each cluster.}
#'    \item{full_data}{data frame of observed area containing cluster information and mean of random effect for each cluster.}
#'    \item{result_sampled}{a list containing result of small area estimation for sampled area, containing following objects: \code{eblup}, \code{fit}, and \code{ref} values of the random effect for each area.}
#'    \item{result_nonsampled}{a list containing result of small area estimation for nonsampled area, containing following objects: \code{eblup} and \code{estcoef}.}
#'    \item{mse_sample}{a list containing mean squared error of sampled area and the values of \code{g1}, referring to \code{g1} in MSE by Prasad-Rao (1990).}
#'    \item{mse_nonsample}{a list containing mean squared error of nonsampled area and the values of \code{g1}, referring to \code{g1} in MSE by Prasad-Rao (1990).}
#'    \item{cluster}{data frame containing cluster information and mean of random effect for each cluster.}
#'  }
#' @examples
#' \donttest{
#' data(nonsample)
#' test <- FHme_nonsamples(formula = y ~ x.hat, var.x = c("var.x"),
#'                         vardir = vardir, n.cluster = 3, data = nonsample)
#' }
#'
#' @export FHme_nonsamples
FHme_nonsamples <- function(formula, var.x, vardir, type.x = "witherror", n.cluster, data){
  formuladata <- model.frame(formula, na.action = NULL, data)
  y <- names(formuladata)[1]
  auxiliary <- formuladata[,-1]
  formula.a <- formula
  var.xa <- var.x
  vardir <- deparse(substitute(vardir))
  type.xa <- type.x

  get_cluster <- as.matrix(kmeans(auxiliary, n.cluster)$cluster)
  data_cluster <- data.frame(data, get_cluster)
  colnames(data_cluster)[ncol(data_cluster)] <- "cluster"

  full <- na.omit(data_cluster)
  used_full <- na.omit(data)
  zero <- data_cluster[which(is.na(data_cluster[,y])),]
  formula_full <- model.frame(formula, data = full)
  y_full <- names(formula_full)[1]
  aux_full <- model.matrix(formula, data = full)


  zero_tmp <- zero
  zero_tmp[,1] <- rnorm(nrow(zero_tmp), mean = 2, sd = 0.1)
  aux_zero <- model.matrix(formula, data = zero_tmp)

  estimated_full <- FHme_edit(formula = formula.a, var.x = var.xa, vardir = vardir,
                              type.x = type.xa, data = used_full)
  ref <- estimated_full$ref
  gamma <- estimated_full$fit$gamma
  g1 <- sapply(1:nrow(full), function(i){
    return(gamma[i]*full[,vardir][i])
  })
  cluster_ref <- data.frame(full, ref, g1)
  byCluster <- aggregate(cluster_ref[,c("ref", "g1")], by = list(cluster_ref$cluster), mean)
  colnames(byCluster) <- c("cluster", "refmean", "g1mean")
  full_r <- left_join(full, byCluster, by = "cluster")

  beta <- estimated_full$fit$estcoef$beta
  refvar <- estimated_full$fit$refvar

  X.beta_zero <- aux_zero %*% beta
  zero_r <- left_join(zero, byCluster, by = "cluster")

  calculate <- data.frame(zero_r, X.beta_zero)
  colnames(calculate)[ncol(calculate)] <- "X.beta"

  EBLUP_zero <- sapply(1:nrow(calculate), function(i){
    cal <- calculate$X.beta[i] + calculate$refmean[i]
    return(cal)
  })
  EBLUP_zero <- data.frame(EBLUP_zero)

  g1_zero <- calculate$g1mean
  g1_zero <- data.frame(g1_zero)
  jk <- lapply(1:nrow(full), function(j) jackknife(full[,y_full],aux_full,full[,vardir],full[,var.xa],j, type.x = type.xa))
  g1_jack <- lapply(1:nrow(full_r), function(i){
    g1 <- sapply(1:nrow(full_r), function(j){
      return(full_r[,vardir][i] * jk[[i]]$result$gamma[j])
    })
    return(data.frame(g1))
  })
  for (i in 1:nrow(full)) {
    jk[[i]]$result$cluster <- full_r$cluster
    jk[[i]]$result$g1 <- g1_jack[[i]]
  }
  Cluster_jack <- lapply(1:nrow(full), function(i){
    tmp <- aggregate(jk[[i]]$result$ref, by = list(jk[[i]]$result$cluster), FUN = mean)
    colnames(tmp) <- c("cluster","refmean")
    tmp2 <- aggregate(jk[[i]]$result$g1, by = list(jk[[i]]$result$cluster), FUN = mean)
    colnames(tmp2) <- c("cluster","g1mean")
    tmp1 <- aux_zero %*% jk[[i]]$beta$betacap
    zero_jack <- left_join(left_join(zero, tmp, by = "cluster"),tmp2, by = "cluster")
    zero_jack$Xbeta <- tmp1
    EBLUP_zero <- sapply(1:nrow(zero_jack), function(j){
      return(zero_jack$Xbeta[j] + zero_jack$refmean[j])
    })
    EBLUP_zero <- data.frame(EBLUP_zero)

    result <- cbind(zero_jack, EBLUP_zero)
    return(result)
  })
  m <- nrow(zero)
  m1cap <- sapply(1:m, function(i){
    left <- g1_zero[i,]
    right <- (m-1)/m * sum(sapply(1:m, function(j){
    jkgamma <- Cluster_jack[[j]]$g1mean[i]
    return(left - jkgamma)
    }))
    return(left + right)
  })

  m2cap <- sapply(1:m, function(i){
    m2 <- ((m-1)/m) * sum(sapply(1:m, function(j)
    {
      return((Cluster_jack[[j]]$EBLUP_zero[i] - EBLUP_zero[i,])^2)
    }))
    return(m2)
  })

  mse_sample <- mse_FHme_e(formula = formula.a, vardir = vardir, var.x = var.xa, type.x = type.xa,
                         data = full)

  msenonsample <- m1cap + m2cap
  mse_nonsample <- list("g1" = g1_zero,
                        "mse" = msenonsample)

  result_nonsampled <- list("eblup" = NA, "estcoef" = NA)
  result_nonsampled$eblup <- EBLUP_zero
  result_nonsampled$estcoef <- estimated_full$fit$estcoef
  zero_r <- zero_r[,-ncol(zero_r)]
  full_r <- full_r[,-ncol(full_r)]
  full_data <- left_join(data_cluster, byCluster, by = "cluster")

  result <- list("nonsampled_data" = zero_r,
                 "sampled_data" = full_r,
                 "full_data" = full_data[,-ncol(full_data)],
                 "result_sampled" = estimated_full,
                 "result_nonsampled" = result_nonsampled,
                 "mse_sample" = mse_sample,
                 "mse_nonsample" = mse_nonsample,
                 "cluster" = byCluster)

  return(result)
}
jackknife <- function(y,X_cap,psi,c,j, type.x = "witherror", w = rep(1,length(y))) {
  m <- length(y)
  p <- dim(X_cap)[2]
  if(type.x == "witherror"){
    c <- as.matrix(c)
  } else if(type.x == "noerror"){
    c <- matrix(0, nrow = dim(X_cap)[1], ncol = p - 1)
  } else {
    c_left <- as.matrix(c)
    c_right <- matrix(0, nrow = dim(X_cap)[1], ncol = (p - 1) - dim(c_left)[2])
    c <- cbind(c_left, c_right)
  }
  c <- cbind(0, c)

  diff_beta <- as.matrix(rep(1,p))
  diff_sigma <- 1
  R_sigma <- 0.001
  R_beta <- as.matrix(rep(0.001,p))
  max_iter <- 100
  betacap_b <- 0
  sigma2cap_b <- 0
  k <- 0
  betacap <- function(y,X_cap,c,w) {
    m <- length(y)
    p <- dim(X_cap)[2]

    wX_capy <- Reduce('+', lapply(1:m, function(i)
    {
      X_cap_i <- as.matrix(X_cap[i,])
      wX_capy <- w[i]*X_cap_i*y[i]
      return(wX_capy)
    }))

    wX_capX_cap <- Reduce('+', lapply(1:m, function(i)
    {
      X_cap_i <- as.matrix(X_cap[i,])
      wX_capX_cap <- w[i]*(X_cap_i%*%t(X_cap_i))
      return(wX_capX_cap)
    }))

    wC <- Reduce('+', lapply(1:m, function(i)
    {
      C <- diag(c[i,],nrow = p)
      wC_i <- w[i]*C
      return(wC_i)
    }))

    betacap <- 0


    chloe <- wX_capX_cap - wC
    if (det(chloe) != 0) {
      Q_matrix <- solve(chloe)
      betacap_a <- Q_matrix %*% wX_capy
      betacap <- list("betacap" = betacap_a,
                      "Q_matrix" = Q_matrix)

    }else{
      G <- wX_capX_cap
      Gsqrt <- sqrtm(G)

      if(sum(eigen(G)$value > 0) == p)
      {
        invG <- solve(Gsqrt)

      }else{
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
      Q1_matrix <- invG%*%pOrth%*%D%*%t(pOrth)%*%invG
      betacap_a <- Q1_matrix%*%wX_capy
      betacap <- list("betacap" = betacap_a,
                      "Q_matrix" = Q1_matrix)
    }

    return(betacap)
  }
  sigma2cap <- function(y,X_cap,c,betacap_i,psi) {
    m <- length(y)
    betacap <- betacap_i
    p <- dim(X_cap)[2]

    pengali <- 1/(m-p)
    poin <- sum(sapply(1:m, function(i)
    {
      X_cap_i <- as.matrix(X_cap[i,])
      poin <- (y[i] - t(X_cap_i)%*%betacap)^2 - psi[i] - t(betacap)%*%diag(c[i,],nrow = p)%*%betacap
      return(poin)
    }))
    sigma2cap <- pengali*poin

    if(sigma2cap < 0) {
      sigma2cap <- 0
    }

    return(sigma2cap)
  }
  gammacap <- function(y,X_cap,betacap_i,sigma2cap,c,psi) {
    m <- length(y)
    p <- dim(X_cap)[2]
    betacap <- betacap_i$betacap
    gammacap <- sapply(1:m, function(i)
    {
      mse_ri <- sigma2cap + t(betacap)%*%diag(c[i,],nrow = p)%*%betacap
      return(mse_ri/(mse_ri + psi[i]))
    })
    return(gammacap)
  }
  y_me <- function(y,X_cap,betacap_i,gammacap) {
    m <- length(y)
    betacap <- betacap_i$betacap
    resid <- y - X_cap %*% betacap
    ref <- gammacap * resid
    yme <- sapply(1:m, function(i){
      yme <- t(as.matrix(X_cap[i,]))%*%betacap + ref[i]
      return(yme)
    })
    sum_up <- list("yme" = yme,
                   "ref" = ref)
    return(sum_up)
  }

  while((any(diff_beta > R_beta) | (diff_sigma > R_sigma)) & (k < max_iter)){
    betacap_a <- betacap_b
    sigma2cap_a <- sigma2cap_b
    if(diff_beta[2] < R_beta[2]){
      sigma2cap_b <- sigma2cap(y[-j],X_cap[-j,],c[-j,],betacap_b,psi[-j])
      w <- sapply(1:m, function(i){
        wi <- 1/(sigma2cap_b + psi[i] + t(betacap_b)%*%diag(c[i,],nrow = p)%*%betacap_b)
        return(wi)
      })
      diff_sigma <- sigma2cap_b - sigma2cap_a

    } else if(diff_sigma < R_sigma){
      betacap_b <- betacap(y[-j],X_cap[-j,],c[-j,],w[-j])$betacap
      w <- sapply(1:m, function(i){
        wi <- 1/(sigma2cap_b + psi[i] + t(betacap_b)%*%diag(c[i,],nrow = p)%*%betacap_b)
        return(wi)
      })
      diff_beta <- betacap_b - betacap_a

    } else {
      betacap_b <- betacap(y[-j],X_cap[-j,],c[-j,],w[-j])$betacap
      sigma2cap_b <- sigma2cap(y[-j],X_cap[-j,],c[-j,],betacap_b,psi[-j])
      w <- sapply(1:m, function(i){
        wi <- 1/(sigma2cap_b + psi[i] + t(betacap_b)%*%diag(c[i,],nrow = p)%*%betacap_b)
        return(wi)
      })
      diff_beta <- betacap_b - betacap_a
      diff_sigma <- sigma2cap_b - sigma2cap_a
    }
    k <- k+1
  }
  betacap_b <- list("betacap" = betacap_b)
  gcap <- gammacap(y,X_cap,betacap_b,sigma2cap_b,c,psi)
  yme <- y_me(y,X_cap,betacap_b,gcap)$yme
  ref <- y_me(y,X_cap,betacap_b,gcap)$ref
  proper <- data.frame(yme, gcap, ref)
  colnames(proper) <- c("EBLUP", "gamma", "ref")
  result <- list("result" = proper,
                 "sigma_2"= sigma2cap_b,
                 'beta' = betacap_b)
  return(result)
}
mse_FHme_e <- function(formula, vardir, var.x, type.x = "witherror", data) {
  #namevar <- deparse(substitute(vardir))
  #name_c <- deparse(substitute(c))
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
    psi <- data[, vardir]
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
  betacap <- function(y,X_cap,c,w) {
    m <- length(y)
    p <- dim(X_cap)[2]

    wX_capy <- Reduce('+', lapply(1:m, function(i)
    {
      X_cap_i <- as.matrix(X_cap[i,])
      wX_capy <- w[i]*X_cap_i*y[i]
      return(wX_capy)
    }))

    wX_capX_cap <- Reduce('+', lapply(1:m, function(i)
    {
      X_cap_i <- as.matrix(X_cap[i,])
      wX_capX_cap <- w[i]*(X_cap_i%*%t(X_cap_i))
      return(wX_capX_cap)
    }))

    wC <- Reduce('+', lapply(1:m, function(i)
    {
      C <- diag(c[i,],nrow = p)
      wC_i <- w[i]*C
      return(wC_i)
    }))

    betacap <- 0


    chloe <- wX_capX_cap - wC
    if (det(chloe) != 0) {
      Q_matrix <- solve(chloe)
      betacap_a <- Q_matrix %*% wX_capy
      betacap <- list("betacap" = betacap_a,
                      "Q_matrix" = Q_matrix)

    }else{
      G <- wX_capX_cap
      Gsqrt <- sqrtm(G)

      if(sum(eigen(G)$value > 0) == p)
      {
        invG <- solve(Gsqrt)

      }else{
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
      Q1_matrix <- invG%*%pOrth%*%D%*%t(pOrth)%*%invG
      betacap_a <- Q1_matrix%*%wX_capy
      betacap <- list("betacap" = betacap_a,
                      "Q_matrix" = Q1_matrix)
    }

    return(betacap)
  }
  sigma2cap <- function(y,X_cap,c,betacap_i,psi) {
    m <- length(y)
    betacap <- betacap_i
    p <- dim(X_cap)[2]

    pengali <- 1/(m-p)
    poin <- sum(sapply(1:m, function(i)
    {
      X_cap_i <- as.matrix(X_cap[i,])
      poin <- (y[i] - t(X_cap_i)%*%betacap)^2 - psi[i] - t(betacap)%*%diag(c[i,],nrow = p)%*%betacap
      return(poin)
    }))
    sigma2cap <- pengali*poin

    if(sigma2cap < 0) {
      sigma2cap <- 0
    }

    return(sigma2cap)
  }
  gammacap <- function(y,X_cap,betacap_i,sigma2cap,c,psi) {
    m <- length(y)
    p <- dim(X_cap)[2]
    betacap <- betacap_i$betacap
    gammacap <- sapply(1:m, function(i)
    {
      mse_ri <- sigma2cap + t(betacap)%*%diag(c[i,],nrow = p)%*%betacap
      return(mse_ri/(mse_ri + psi[i]))
    })
    return(gammacap)
  }
  beta_sigma_conv <- function(y,X_cap,psi,c, w = rep(1,length(y))) {
    m <- length(y)
    p <- dim(X_cap)[2]
    sigma2cap_b <- 0
    betacap_b <- 0
    R_sigma <- 0.001
    R_beta <- as.matrix(rep(0.001,p))
    max_iter <- 1000
    k <- 0
    diff_beta <- as.matrix(rep(1,p))
    diff_sigma <- 1
    convergence <- TRUE

    while(((any(diff_beta > R_beta)) | (diff_sigma > R_sigma)) & (k < max_iter)){
      betacap_a <- betacap_b
      sigma2cap_a <- sigma2cap_b
      if(all(diff_beta < R_beta)){
        sigma2cap_b <- sigma2cap(y,X_cap,c,betacap_b,psi)
        w <- sapply(1:m, function(i){
          wi <- 1/(sigma2cap_b + psi[i] + t(betacap_b)%*%diag(c[i,],nrow = p)%*%betacap_b)
          return(wi)
        })
        diff_sigma <- sigma2cap_b - sigma2cap_a

      } else if(diff_sigma < R_sigma){
        betacap_b <- betacap(y,X_cap,c,w)$betacap
        Q_matrix <- betacap(y, X_cap, c, w)$Q_matrix
        w <- sapply(1:m, function(i){
          wi <- 1/(sigma2cap_b + psi[i] + t(betacap_b)%*%diag(c[i,],nrow = p)%*%betacap_b)
          return(wi)
        })
        diff_beta <- betacap_b - betacap_a

      } else {
        betacap_b <- betacap(y,X_cap,c,w)$betacap
        Q_matrix <- betacap(y, X_cap, c, w)$Q_matrix
        sigma2cap_b <- sigma2cap(y,X_cap,c,betacap_b,psi)
        w <- sapply(1:m, function(i){
          wi <- 1/(sigma2cap_b + psi[i] + t(betacap_b)%*%diag(c[i,],nrow = p)%*%betacap_b)
          return(wi)
        })
        diff_beta <- betacap_b - betacap_a
        diff_sigma <- sigma2cap_b - sigma2cap_a
      }
      k <- k+1

    }
    if (k >= max_iter & ((any(diff_beta > R_beta)) | (diff_sigma > R_sigma))) {
      convergence <- FALSE
    }
    beta_sigma <- list('betacap'= betacap_b,
                       'sigma2cap'= sigma2cap_b,
                       "Q_matrix" = Q_matrix,
                       "convergence" = convergence,
                       "iterations" = k)
    return(beta_sigma)
  }

  beta_sigma <- beta_sigma_conv(y,X_cap,psi,c)
  betacap_b <- beta_sigma
  sigma2cap_b <- beta_sigma$sigma2cap
  gcap <- gammacap(y,X_cap,betacap_b,sigma2cap_b,c,psi)

  y_me <- function(y,X_cap,betacap_i,gammacap) {
    m <- length(y)
    betacap <- betacap_i$betacap
    yme <- sapply(1:m, function(i){
      yme <- gammacap[i]*y[i] + (1-gammacap[i])*t(as.matrix(X_cap[i,]))%*%betacap
      return(yme)
    })
    return(yme)
  }

  yme <- y_me(y,X_cap,betacap_b,gcap)
  yME <- list("y_me" = yme,
              "gamma" = gcap)

  jackknife <- function(y,X_cap,psi,c,j, w = rep(1,length(y))) {
    m <- length(y)
    p <- dim(X_cap)[2]
    diff_beta <- as.matrix(rep(1,p))
    diff_sigma <- 1
    R_sigma <- 0.001
    R_beta <- as.matrix(rep(0.001,p))
    max_iter <- 100
    betacap_b <- 0
    sigma2cap_b <- 0
    k <- 0

    while((any(diff_beta > R_beta) | (diff_sigma > R_sigma)) & (k < max_iter)){
      betacap_a <- betacap_b
      sigma2cap_a <- sigma2cap_b
      if(diff_beta[2] < R_beta[2]){
        sigma2cap_b <- sigma2cap(y[-j],X_cap[-j,],c[-j,],betacap_b,psi[-j])
        w <- sapply(1:m, function(i){
          wi <- 1/(sigma2cap_b + psi[i] + t(betacap_b)%*%diag(c[i,],nrow = p)%*%betacap_b)
          return(wi)
        })
        diff_sigma <- sigma2cap_b - sigma2cap_a

      } else if(diff_sigma < R_sigma){
        betacap_b <- betacap(y[-j],X_cap[-j,],c[-j,],w[-j])$betacap
        w <- sapply(1:m, function(i){
          wi <- 1/(sigma2cap_b + psi[i] + t(betacap_b)%*%diag(c[i,],nrow = p)%*%betacap_b)
          return(wi)
        })
        diff_beta <- betacap_b - betacap_a

      } else {
        betacap_b <- betacap(y[-j],X_cap[-j,],c[-j,],w[-j])$betacap
        sigma2cap_b <- sigma2cap(y[-j],X_cap[-j,],c[-j,],betacap_b,psi[-j])
        w <- sapply(1:m, function(i){
          wi <- 1/(sigma2cap_b + psi[i] + t(betacap_b)%*%diag(c[i,],nrow = p)%*%betacap_b)
          return(wi)
        })
        diff_beta <- betacap_b - betacap_a
        diff_sigma <- sigma2cap_b - sigma2cap_a
      }
      k <- k+1
    }
    betacap_b <- list("betacap" = betacap_b)
    gcap <- gammacap(y,X_cap,betacap_b,sigma2cap_b,c,psi)
    yme <- y_me(y,X_cap,betacap_b,gcap)
    result <- list('y_me' = yme,
                   'gamma' = gcap,
                   'sigma^2' = sigma2cap_b,
                   'beta' = betacap_b)
    return(result)
  }

  jk <- lapply(1:m, function(j) jackknife(y,X_cap,psi,c,j))

  m1cap <- sapply(1:m, function(i)
  {
    left <- yME$gamma[i]*psi[i]

    right <- ((m-1)/m) * sum(sapply(1:m, function(j){
      jkgamma <- jk[[j]]$gamma[i]
      return(left-jkgamma*psi[i])
    }))
    return(left+right)
  })

  m2cap <- sapply(1:m, function(i){
    m2 <- ((m-1)/m) * sum(sapply(1:m, function(j)
    {
      return((jk[[j]]$y_me[i] - yME$y_me[i])^2)
    }))
    return(m2)
  })

  mse <- m1cap + m2cap
  g1 <- sapply(1:m, function(i){
    return(yME$gamma[i]*psi[i])
  })
  return(list("g1" = g1,
              "mse" = mse))
}
FHme_edit <- function(formula, vardir, var.x, type.x = "witherror", data) {
  #namevar <- deparse(substitute(vardir))
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
    psi <- data[, vardir]
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

  betacap <- function(y,X_cap,c,w) {
    m <- length(y)
    p <- dim(X_cap)[2]

    wX_capy <- Reduce('+', lapply(1:m, function(i)
    {
      X_cap_i <- as.matrix(X_cap[i,])
      wX_capy <- w[i]*X_cap_i*y[i]
      return(wX_capy)
    }))

    wX_capX_cap <- Reduce('+', lapply(1:m, function(i)
    {
      X_cap_i <- as.matrix(X_cap[i,])
      wX_capX_cap <- w[i]*(X_cap_i%*%t(X_cap_i))
      return(wX_capX_cap)
    }))

    wC <- Reduce('+', lapply(1:m, function(i)
    {
      C <- diag(c[i,],nrow = p)
      wC_i <- w[i]*C
      return(wC_i)
    }))

    betacap <- 0

    chloe <- wX_capX_cap - wC
    if (det(chloe) != 0) {
      Q_matrix <- solve(chloe)
      betacap_a <- Q_matrix %*% wX_capy
      betacap <- list("betacap" = betacap_a,
                      "Q_matrix" = Q_matrix)

    }else{
      G <- wX_capX_cap
      Gsqrt <- sqrtm(G)

      if(sum(eigen(G)$value > 0) == p)
      {
        invG <- solve(Gsqrt)

      }else{
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
      Q1_matrix <- invG%*%pOrth%*%D%*%t(pOrth)%*%invG
      betacap_a <- Q1_matrix%*%wX_capy
      betacap <- list("betacap" = betacap_a,
                      "Q_matrix" = Q1_matrix)
    }

    return(betacap)
  }
  sigma2cap <- function(y,X_cap,c,betacap_i,psi) {
    m <- length(y)
    betacap <- betacap_i
    p <- dim(X_cap)[2]

    pengali <- 1/(m-p)
    poin <- sum(sapply(1:m, function(i)
    {
      X_cap_i <- as.matrix(X_cap[i,])
      poin <- (y[i] - t(X_cap_i)%*%betacap)^2 - psi[i] - t(betacap)%*%diag(c[i,],nrow = p)%*%betacap
      return(poin)
    }))
    sigma2cap <- pengali*poin
    #Ybarra & Lohr 2008: if estimated random effect < 0 set to zero
    if(sigma2cap < 0) {
      sigma2cap <- 0
    }

    return(sigma2cap)
  }

  beta_sigma_conv <- function(y,X_cap,psi,c, w = rep(1,length(y))) {
    m <- length(y)
    p <- dim(X_cap)[2]
    sigma2cap_b <- 0
    betacap_b <- 0
    R_sigma <- 0.001
    R_beta <- as.matrix(rep(0.001,p))
    max_iter <- 1000
    k <- 0
    diff_beta <- as.matrix(rep(1,p))
    diff_sigma <- 1
    convergence <- TRUE

    while(((any(diff_beta > R_beta)) | (diff_sigma > R_sigma)) & (k < max_iter)){
      betacap_a <- betacap_b
      sigma2cap_a <- sigma2cap_b
      if(all(diff_beta < R_beta)){
        sigma2cap_b <- sigma2cap(y,X_cap,c,betacap_b,psi)
        w <- sapply(1:m, function(i){
          wi <- 1/(sigma2cap_b + psi[i] + t(betacap_b)%*%diag(c[i,],nrow = p)%*%betacap_b)
          return(wi)
        })
        diff_sigma <- sigma2cap_b - sigma2cap_a

      } else if(diff_sigma < R_sigma){
        betacap_b <- betacap(y,X_cap,c,w)$betacap
        Q_matrix <- betacap(y, X_cap, c, w)$Q_matrix
        w <- sapply(1:m, function(i){
          wi <- 1/(sigma2cap_b + psi[i] + t(betacap_b)%*%diag(c[i,],nrow = p)%*%betacap_b)
          return(wi)
        })
        diff_beta <- betacap_b - betacap_a

      } else {
        betacap_b <- betacap(y,X_cap,c,w)$betacap
        Q_matrix <- betacap(y, X_cap, c, w)$Q_matrix
        sigma2cap_b <- sigma2cap(y,X_cap,c,betacap_b,psi)
        w <- sapply(1:m, function(i){
          wi <- 1/(sigma2cap_b + psi[i] + t(betacap_b)%*%diag(c[i,],nrow = p)%*%betacap_b)
          return(wi)
        })
        diff_beta <- betacap_b - betacap_a
        diff_sigma <- sigma2cap_b - sigma2cap_a
      }
      k <- k+1

    }
    if (k >= max_iter & ((any(diff_beta > R_beta)) | (diff_sigma > R_sigma))) {
      convergence <- FALSE
    }
    beta_sigma <- list('betacap'= betacap_b,
                       'sigma2cap'= sigma2cap_b,
                       "Q_matrix" = Q_matrix,
                       "convergence" = convergence,
                       "iterations" = k)
    return(beta_sigma)
  }

  gammacap <- function(y,X_cap,betacap_i,sigma2cap,c,psi) {
    m <- length(y)
    p <- dim(X_cap)[2]
    betacap <- betacap_i$betacap
    gammacap <- sapply(1:m, function(i)
    {
      mse_ri <- sigma2cap + t(betacap)%*%diag(c[i,],nrow = p)%*%betacap
      return(mse_ri/(mse_ri + psi[i]))
    })
    return(gammacap)
  }
  y_me <- function(y,X_cap,betacap_i,gammacap) {
    m <- length(y)
    betacap <- betacap_i$betacap
    resid <- y - X_cap %*% betacap
    ref <- gammacap * resid
    yme <- sapply(1:m, function(i){
      yme <- t(as.matrix(X_cap[i,]))%*%betacap + ref[i]
      return(yme)
    })
    sum_up <- list("yme" = yme,
                   "ref" = ref)
    return(sum_up)
  }

  beta_sigma <- beta_sigma_conv(y,X_cap,psi,c)
  se.b <- sqrt(diag(beta_sigma$Q_matrix))
  betacap_b <- beta_sigma
  t.val <- betacap_b$betacap/se.b
  pv <- 2 * pnorm(abs(t.val), lower.tail = FALSE)
  coef <- data.frame(betacap_b$betacap, se.b, t.val, pv)
  colnames(coef) <- c("beta", "std.error", "t.statistics", "p.value")
  sigma2cap_b <- beta_sigma$sigma2cap

  gcap <- gammacap(y,X_cap,betacap_b,sigma2cap_b,c,psi)
  yme <- y_me(y,X_cap,betacap_b,gcap)$yme
  ref <- y_me(y,X_cap,betacap_b,gcap)$ref
  result <- list(eblup = NA, fit = list(method = NA, convergence = NA,
                                        iterations = NA,
                                        estcoef =NA,
                                        refvar = NA,
                                        gamma = NA))
  result$eblup <- data.frame(yme)
  result$ref <- ref
  result$fit$method <- 'REML'
  result$fit$convergence <- beta_sigma$convergence
  if (result$fit$convergence == FALSE) {
    warning("The fitting method does not converge.\n")
  }
  result$fit$iterations <- beta_sigma$iterations
  result$fit$estcoef <- coef
  result$fit$refvar <- sigma2cap_b
  result$fit$gamma <- gcap

  return(result)
}
