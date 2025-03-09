#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### SRM EEG analysis
#-----------------------------------------------------------------------------#
### Data taken from: https://openneuro.org/datasets/ds003775/versions/1.2.1
### - Specifically we use the derivatives/cleaned_epochs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#Important note: install Difdtl using devtools::install_github("SusanYuan/Difdtl"). 
#  For Mac, requires R <= 4.1.2. The C compiler changed in later versions of R and code does not compile

## If run in Rstudio
library(Difdtl)
library(R.matlab)
library(igraph)
library(ggplot2)
library(glassoFast)
library(astsa)
library(data.table)
library(JGL)
res_folder <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/results")
setwd(paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/data/derivatives/cleaned_epochs"))

## If run in CMD line - run in srm_eeg64 folder
args <- commandArgs(trailingOnly = F)
stopifnot(length(args)==1) #Only 1 arg and it should be subj
subj <- args[length(args)]
subj <- sub("-","",subj)
## Load libraries here
res_folder <- paste0(getwd(), "/results")
setwd("./data/derivatives/cleaned_epochs")

#~~~~~~~~~~~~~~~~~~~~#
#     Helper fns     #
#~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### Code for JGL but allow to take list of S (covariances) instead of list of Y
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### - Editing to do FGL penalty only and takes S = list of K cov matrices (or in our case K expanded spectral densities)
JGL_cov <- function (S, penalty = "fused", lambda1, lambda2, n = NA, rho = 1, weights = "equal", 
                     penalize.diagonal = FALSE, maxiter = 500, tol = 1e-05, warm = NULL, 
                     return.whole.theta = FALSE, screening = "fast", truncate = 1e-05) 
{
  p = dim(S[[1]])[2]
  K = length(S)
  
  # for (k in 1:K) {
  #   n[k] = dim(Y[[k]])[1]
  # }
  if (length(dimnames(S[[1]])[[2]]) == 0) {
    for (k in 1:K) {
      dimnames(S[[k]])[[2]] = paste("V", 1:p, sep = "")
    }
  }
  # for (k in 1:K) {
  #   for (j in 1:p) {
  #     Y[[k]][, j] = Y[[k]][, j] - mean(Y[[k]][, j])
  #   }
  # }
  if (length(weights) == 1) {
    if (weights == "equal") {
      weights = rep(1, K)
    }
  }
  if ((length(weights) == 1) & (!is.na(n))) {
    if (weights == "sample.size") {
      weights = n/sum(n)
    }
  }
  # if (screening == "memory.efficient") {
  #   if (penalty == "fused") {
  #     connected = screen.fgl(Y, lambda1, lambda2, weights)
  #   }
  #   if (penalty == "group") {
  #     connected = screen.ggl(Y, lambda1, lambda2, weights)
  #   }
  # }
  if (screening == "fast") {
    connected = rep(TRUE, p)
  }
  # if (!((screening == "memory.efficient") | (screening == "fast"))) {
  #   stop("screening must equal \"fast\" or \"memory.efficient\".")
  # }
  # S = vector("list", length = K)
  # for (k in 1:K) {
  #   ntemp = dim(Y[[k]])[1]
  #   S[[k]] = cov(Y[[k]][, connected]) * (ntemp - 1)/ntemp
  # }
  lam1 = lambda1
  lam2 = lambda2
  if (length(lam1) > 1) {
    lam1 = lam1[connected, connected]
  }
  if (length(lam2) > 1) {
    lam2 = lam2[connected, connected]
  }
  if (penalty == "fused") {
    if (K == 2) {
      crit1 = list()
      for (k in 1:K) {
        crit1[[k]] = abs(S[[k]]) * weights[k] > lam1 + 
          lam2
      }
      S.sum = matrix(0, sum(connected), sum(connected))
      for (k in 1:K) {
        S.sum = S.sum + weights[k] * S[[k]]
      }
      S.sum = abs(S.sum)
      crit2 = S.sum > 2 * lam1
    }
    if (K > 2) {
      crit1 = list()
      for (k in 1:K) {
        crit1[[k]] = abs(S[[k]]) * weights[k] > lam1
      }
      crit2 = matrix(0, sum(connected), sum(connected))
    }
    critboth = crit2
    for (k in 1:K) {
      critboth = critboth + crit1[[k]]
    }
    
    critboth = (critboth != 0)
    diag(critboth) = 1
  }
  if (penalty == "group") {
    tempsum = matrix(0, sum(connected), sum(connected))
    for (k in 1:K) {
      tempsum = tempsum + (pmax(weights[k] * abs(S[[k]]) - 
                                  lam1, 0))^2
    }
    critboth = tempsum > lam2^2
    diag(critboth) = 1
  }
  g1 <- graph_from_adjacency_matrix(critboth)
  cout = components(g1)
  blocklist = list()
  unconnected = c()
  if (min(cout$membership) == 0) {
    cout$membership = cout$membership + 1
  }
  for (i in 1:(cout$no)) {
    if (sum(cout$membership == i) == 1) {
      unconnected <- c(unconnected, which(cout$membership == 
                                            i))
    }
    if (sum(cout$membership == i) > 1) {
      blocklist[[length(blocklist) + 1]] <- which(cout$membership == 
                                                    i)
    }
  }
  connected[unconnected] = FALSE
  connected.index = rep(0, p)
  connected.index[connected] = 1:sum(connected)
  unconnected = !connected
  theta = list()
  for (k in 1:K) {
    theta[[k]] = matrix(0, sum(connected), sum(connected))
    if (sum(connected) > 0) {
      dimnames(theta[[k]])[[1]] = dimnames(theta[[k]])[[2]] = dimnames(S[[k]])[[2]][connected]
    }
  }
  Su = list()
  for (k in 1:K) {
    Su[[k]] = diag(S[[k]])[unconnected]
  }
  if (length(lambda1) == 1) {
    lam1.unconnected = lambda1
  }
  if (length(lambda1) > 1) {
    lam1.unconnected = diag(lambda1)[unconnected]
  }
  if (length(lambda2) == 1) {
    lam2.unconnected = lambda2
  }
  if (length(lambda2) > 1) {
    lam2.unconnected = diag(lambda2)[unconnected]
  }
  if (!penalize.diagonal) {
    lam1.unconnected = lam1.unconnected * 0
    if (penalty == "group") {
      lam2.unconnected = lam2.unconnected * 0
    }
  }
  if (sum(unconnected) > 0) {
    theta.unconnected = admm.iters.unconnected.cov(Su, p = sum(unconnected), lambda1 = lam1.unconnected, 
                                                   lambda2 = lam2.unconnected, penalty = penalty, rho = rho, 
                                                   weights = weights, maxiter = maxiter, tol = tol)$Z
    for (k in 1:K) {
      names(theta.unconnected[[k]]) = dimnames(S[[k]])[[2]][!connected]
    }
  }
  if (sum(unconnected) == 0) {
    theta.unconnected = NULL
  }
  if (length(blocklist) > 0) {
    for (i in 1:length(blocklist)) {
      bl <- blocklist[[i]]
      Sbl = list()
      for (k in 1:K) {
        Sbl[[k]] = S[[k]][bl, bl]
      }
      if (length(lambda1) == 1) {
        lam1.bl = lambda1
      }
      if (length(lambda1) > 1) {
        lam1.bl = lambda1[bl, bl]
      }
      if (length(lambda2) == 1) {
        lam2.bl = lambda2
      }
      if (length(lambda2) > 1) {
        lam2.bl = lambda2[bl, bl]
      }
      
      lam1.bl = JGL:::penalty.as.matrix(lam1.bl, dim(Sbl[[1]])[2], 
                                        penalize.diagonal = penalize.diagonal)
      if (penalty == "fused") {
        lam2.bl = JGL:::penalty.as.matrix(lam2.bl, dim(Sbl[[1]])[2], 
                                          penalize.diagonal = TRUE)
      }
      if (penalty == "group") {
        lam2.bl = JGL:::penalty.as.matrix(lam2.bl, dim(Sbl[[1]])[2], 
                                          penalize.diagonal = penalize.diagonal)
      }
      if (length(warm) == 0) {
        warm.bl = NULL
      }
      if (length(warm) > 0) {
        warm.bl = list()
        for (k in 1:K) {
          warm.bl[[k]] = warm[[k]][bl, bl]
        }
      }
      
      Thetabl = admm.iters.cov(Sbl, lam1.bl, lam2.bl, penalty = penalty, 
                               rho = rho, weights = weights, penalize.diagonal = TRUE, 
                               maxiter = maxiter, tol = tol, warm = warm.bl)
      for (k in 1:K) {
        theta[[k]][connected.index[bl], connected.index[bl]] = Thetabl$Z[[k]]
      }
    }
  }
  if (dim(theta[[1]])[1] > 0) {
    for (k in 1:K) {
      rounddown = abs(theta[[k]]) < truncate
      diag(rounddown) = FALSE
      theta[[k]] = theta[[k]] * (1 - rounddown)
    }
  }
  if (!return.whole.theta) {
    out = list(theta = theta, theta.unconnected = theta.unconnected, 
               connected = connected)
  }
  if (return.whole.theta) {
    whole.theta = list()
    for (k in 1:K) {
      whole.theta[[k]] = matrix(0, p, p)
      diag(whole.theta[[k]])[unconnected] = theta.unconnected[[k]]
      whole.theta[[k]][connected, connected] = theta[[k]]
      dimnames(whole.theta[[k]])[[1]] = dimnames(whole.theta[[k]])[[2]] = dimnames(S[[k]])[[2]]
    }
    out = list(theta = whole.theta, connected = connected)
  }
  class(out) = "jgl"
  return(out)
}


### Need to edit this function to take list of S, not Y shouldn't be too bad to do this
admm.iters.cov <- function (S, lambda1, lambda2, penalty = "fused", rho = 1, rho.increment = 1, 
                            weights, penalize.diagonal, maxiter = 1000, tol = 1e-05, 
                            warm = NULL) 
{
  
  K = length(S)
  p = dim(S[[1]])[2]
  n = weights
  # ns = c()
  # for (k in 1:K) {
  #   ns[k] = dim(Y[[k]])[1]
  # }
  # S = list()
  # for (k in 1:K) {
  #   S[[k]] = cov(Y[[k]]) * (ns[k] - 1)/ns[k]
  # }
  theta = list()
  for (k in 1:K) {
    theta[[k]] = diag(1/diag(S[[k]]))
  }
  Z = list()
  for (k in 1:K) {
    Z[[k]] = matrix(0, p, p)
  }
  W = list()
  for (k in 1:K) {
    W[[k]] = matrix(0, p, p)
  }
  
  lam1 = JGL:::penalty.as.matrix(lambda1, p, penalize.diagonal = penalize.diagonal)
  if (penalty == "fused") {
    lam2 = JGL:::penalty.as.matrix(lambda2, p, penalize.diagonal = TRUE)
  }
  if (penalty == "group") {
    lam2 = JGL:::penalty.as.matrix(lambda2, p, penalize.diagonal = penalize.diagonal)
  }
  
  iter = 0
  diff_value = 10
  while ((iter == 0) || (iter < maxiter && diff_value > tol)) {
    
    if (FALSE) {
      print(paste("iter=", iter))
      if (penalty == "fused") {
        print(paste("crit=", crit(theta, S, n = rep(1, 
                                                    K), lam1, lam2, penalize.diagonal = penalize.diagonal)))
        print(paste("crit=", crit(Z, S, n = rep(1, K), 
                                  lam1, lam2, penalize.diagonal = penalize.diagonal)))
      }
      if (penalty == "group") {
        print(paste("crit=", gcrit(theta, S, n = rep(1, 
                                                     K), lam1, lam2, penalize.diagonal = penalize.diagonal)))
      }
    }
    theta.prev = theta
    for (k in 1:K) {
      edecomp_mat <- S[[k]] - rho * Z[[k]]/n[k] + rho * W[[k]]/n[k]
      max_t_diff <- max(abs(edecomp_mat - t(edecomp_mat)))
      if(max_t_diff > 1e-12){stop("Maximum diff between edecomp_mat and t(edecomp_mat) is: ", max_t_diff)}
      ## Enforce symmetry so we don't get complex eigenvalues due to numerical issues
      edecomp_mat <- (edecomp_mat + t(edecomp_mat))/2
      edecomp = eigen(edecomp_mat)
      D = edecomp$values
      V = edecomp$vectors
      D2 = n[k]/(2 * rho) * (-D + sqrt(D^2 + 4 * rho/n[k]))
      theta[[k]] = V %*% diag(D2) %*% t(V)
    }
    A = list()
    for (k in 1:K) {
      A[[k]] = theta[[k]] + W[[k]]
    }
    if (penalty == "fused") {
      if (K == 2) {
        Z = JGL:::flsa2(A, rho, lam1, lam2, penalize.diagonal = TRUE)
      }
      if (K > 2) {
        Z = JGL:::flsa.general(A, rho, lam1, lam2, penalize.diagonal = TRUE)
      }
    }
    if (penalty == "group") {
      Z = JGL:::dsgl(A, rho, lam1, lam2, penalize.diagonal = TRUE)
    }
    for (k in 1:K) {
      W[[k]] = W[[k]] + (theta[[k]] - Z[[k]])
    }
    iter = iter + 1
    diff_value = 0
    for (k in 1:K) {
      diff_value = diff_value + sum(abs(theta[[k]] - theta.prev[[k]]))/sum(abs(theta.prev[[k]]))
    }
    rho = rho * rho.increment
  }
  diff = 0
  for (k in 1:K) {
    diff = diff + sum(abs(theta[[k]] - Z[[k]]))
  }
  out = list(theta = theta, Z = Z, diff = diff, iters = iter)
  return(out)
}


### Need to edit this function to take list of S, not Y shouldn't be too bad to do this
### S should be list where each element is vector of variances of the unconnected nodes
### p = # of unconnected nodes
admm.iters.unconnected.cov <- function (S, p, lambda1, lambda2, penalty = "fused", rho = 1, rho.increment = 1, 
                                        weights, maxiter = 1000, tol = 1e-05, warm = NULL) 
{
  K = length(S)
  # for (k in 1:K) {
  #   S[[k]] = as.matrix(S[[k]])
  # }
  # p = dim(S[[1]])[2]
  n = weights
  # S = list()
  # for (k in 1:K) {
  #   ntemp = dim(Y[[k]])[1]
  #   S[[k]] = rep(0, p)
  #   for (j in 1:p) {
  #     S[[k]][j] = var(Y[[k]][, j]) * (ntemp - 1)/ntemp
  #   }
  # }
  theta = list()
  for (k in 1:K) {
    theta[[k]] = 1/S[[k]]
  }
  Z = list()
  for (k in 1:K) {
    Z[[k]] = rep(0, p)
  }
  W = list()
  for (k in 1:K) {
    W[[k]] = rep(0, p)
  }
  lam1 = lambda1
  lam2 = lambda2
  iter = 0
  diff_value = 10
  while ((iter == 0) || (iter < maxiter && diff_value > tol)) {
    theta.prev = theta
    for (k in 1:K) {
      B = n[k] * S[[k]] - rho * (Z[[k]] - W[[k]])
      theta[[k]] = 1/(2 * rho) * (-B + sqrt(B^2 + 4 * rho * 
                                              n[k]))
    }
    A = list()
    for (k in 1:K) {
      A[[k]] = theta[[k]] + W[[k]]
    }
    if (penalty == "fused") {
      if (K == 2) {
        Z = JGL:::flsa2(A, rho, lam1, lam2, penalize.diagonal = TRUE)
      }
      if (K > 2) {
        Z = JGL:::flsa.general(A, rho, lam1, lam2, penalize.diagonal = TRUE)
      }
    }
    if (penalty == "group") {
      Z = JGL:::dsgl(A, rho, lam1, lam2, penalize.diagonal = TRUE)
    }
    for (k in 1:K) {
      W[[k]] = W[[k]] + (theta[[k]] - Z[[k]])
    }
    iter = iter + 1
    diff_value = 0
    for (k in 1:K) {
      diff_value = diff_value + sum(abs(theta[[k]] - theta.prev[[k]]))/sum(abs(theta.prev[[k]]))
    }
    rho = rho * rho.increment
  }
  diff = 0
  for (k in 1:K) {
    diff = diff + sum(abs(theta[[k]] - Z[[k]]))
  }
  out = list(theta = theta, Z = Z, diff = diff, iters = iter)
  return(out)
}

## AIC for FGL
aic_FGL <- function(S_l, theta_l, n_l){
  k <- length(S_l)
  aics <- sapply(1:k, function(cond){
    S <- S_l[[cond]]
    theta <- theta_l[[cond]]
    n <- n_l[cond]
    E <- (abs(theta)>1e-10)
    #E <- E[upper.tri(E)] ## The paper says to use the number of non-zero entries so do that
    edges <- sum(E)
    n*sum(diag(S %*% theta)) - n*log(det(theta)) + 2*edges
  })
  return(sum(aics))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### End JGL code
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

offdiags <- function(X){
  return(X[row(X) != col(X)])
}


pickfxx <- function(mvspec_res, freqi){
  specfreqs <- mvspec_res$freq
  use_idx <- if(any(abs(freqi - specfreqs)<1e-10)){which(abs(freqi - specfreqs)<1e-10)} else{ max(which(freqi > specfreqs))}
  return(list(fxx = mvspec_res$fxx[,,use_idx], freq = specfreqs[use_idx], orig.n = mvspec_res$orig.n))
}

calc_eBIC <- function(sigX, sigY, delta, nx, ny, isComplex = TRUE, edgeMag = 1e-8){
  gamma <- 0.5
  p_ebic <- nrow(delta)
  
  if(isComplex){
    p <- nrow(delta)/2
    reDelta <- delta[1:p,1:p]
    imDelta <- delta[1:p,(p+1):(2*p)]
    n_edges <- sum(abs(reDelta[upper.tri(reDelta, diag = TRUE)]) > edgeMag) + sum(abs(imDelta[upper.tri(imDelta, diag = TRUE)]) > edgeMag)
  } else{
    n_edges <- sum(abs(delta[upper.tri(delta, diag = TRUE)]) > edgeMag)
  }
  
  ee <- (1/2)*((sigX %*% delta %*% sigY) + (sigY %*% delta %*% sigX)) - (sigX - sigY)
  ee_max <- max(abs(ee))
  min_ebic_max_correct <- min(nx,ny)*ee_max + log(min(nx,ny))*n_edges + 4*n_edges*gamma*log(p_ebic)
  
  return(min_ebic_max_correct)
}


eBIC_glasso <- function(theta, S, n, gamma){
  p <- nrow(S)
  E <- (abs(theta)>1e-10)
  E <- E[upper.tri(E)]
  cardE <- sum(E)
  log_det <- determinant(theta, logarithm = TRUE)$modulus
  attributes(log_det) <- NULL
  lik <- (n/2)*(log_det - sum(diag(theta %*% S)))
  return((-2*lik + cardE*log(n) + 4*cardE*gamma*log(p))) 
}


# Hard vs soft thresholding
threshold <- function(X, thr, hard){
  if(hard){
    X[abs(X) < thr] <- 0
  } else{
    X <- ifelse(abs(X) < thr, 0, ifelse(X > 0, X - thr, X + thr))
  }
  return(X)
}

runDtl <- function(sigX, sigY, nx, ny, nlambdas = 20){
  p <- dim(sigX)[1]
  sigXcor <- diag(diag(sigX)^(-1/2)) %*% sigX %*% diag(diag(sigX)^(-1/2))
  sigYcor <- diag(diag(sigY)^(-1/2)) %*% sigY %*% diag(diag(sigY)^(-1/2))
  
  realSigX <- rbind(cbind(Re(sigXcor), -Im(sigXcor)), cbind(Im(sigXcor), Re(sigXcor)))
  realSigY <- rbind(cbind(Re(sigYcor), -Im(sigYcor)), cbind(Im(sigYcor), Re(sigYcor)))
  
  realSigX <- (realSigX + t(realSigX))/2
  realSigY <- (realSigY + t(realSigY))/2
  
  max_lam <- max(abs(realSigX - realSigY))*2
  lam_seq <- exp(seq(log(max_lam*0.001), log(max_lam), length.out = nlambdas))
  
  dtl_ests <- sapply(1:nlambdas, function(i){
    print(paste0("---------- DTL Lambda: ", i))
    lam <- lam_seq[i]
    est <- L1_dts(SigmaX = round(realSigX,8), SigmaY = round(realSigY,8), rho = 1, lambda = lam) #Remeber L1_dts returns est of SigmaY^(-1) - SigmaX^(-1). See ?L1_dts. The example makes this clear
    bic <- calc_eBIC(realSigX, realSigY, est, nx = nx, ny = ny)
    return(list(est = est, bic = bic, lambda = lam, lam_num = i))
  }, simplify = F, USE.NAMES = T)
  
  bics <- sapply(dtl_ests, "[[", "bic")
  stopifnot(sum(bics == min(bics))==1)
  lam_opt_num <- which.min(bics)
  opt_est <- dtl_ests[[lam_opt_num]]$est
  
  
  ## Naive GLASSO estimator
  naive_lam_maxX <- max(abs(realSigX))
  naive_lam_maxY <- max(abs(realSigY))

  naive_lam_seqX <- exp(seq(log(naive_lam_maxX*0.001), log(naive_lam_maxX), length.out = nlambdas))
  naive_lam_seqY <- exp(seq(log(naive_lam_maxY*0.001), log(naive_lam_maxY), length.out = nlambdas))

  naiveXs <- lapply(naive_lam_seqX, function(lamX){
    thetaX <- glassoFast::glassoFast(realSigX, rho = lamX)
    ebic <- eBIC_glasso(thetaX$wi, realSigX, nx, gamma = 0.5) #gamma = 0.5 based on eBIC paper - this seems to give decent performance in paper
    return(list(thetaX = thetaX, ebic = ebic))
  })
  naiveYs <- lapply(naive_lam_seqY, function(lamY){
    thetaY <- glassoFast::glassoFast(realSigY, rho = lamY)
    ebic <- eBIC_glasso(thetaY$wi, realSigY, ny, gamma = 0.5) #gamma = 0.5 based on eBIC paper - this seems to give decent performance in paper
    return(list(thetaY = thetaY, ebic = ebic))
  })

  eBICXnum = which.min(sapply(naiveXs, "[[", "ebic"))
  eBICYnum = which.min(sapply(naiveYs, "[[", "ebic"))

  naiveX = naiveXs[[eBICXnum]]$thetaX$wi
  naiveY = naiveYs[[eBICYnum]]$thetaY$wi
  naiveGLASSODiff = naiveY - naiveX
  naiveGLASSODiff = naiveGLASSODiff[1:p,1:p] + complex(imaginary = naiveGLASSODiff[(p+1):(2*p), 1:p])
  
  ## FGL estimator
  nlambda1 <- 2
  max_lambda1 <- max(abs(offdiags(realSigX) + offdiags(realSigY)))/2 ## point where should just be disconnected (atleast for this condition, for actual disconnected we'd need to satisfy other conditions too)
  lambda1s <- seq(0.01, 0.1, length.out = nlambda1)*max_lambda1
  
  nlambda2 <- nlambdas/nlambda1
  max_lambda2 <- max(max(abs(offdiags(realSigX))) - min(lambda1s), max(abs(offdiags(realSigY))) - min(lambda1s)) ## point where just disconnected (for these)
  lambda2s <-  exp(seq(log(max_lambda2*0.0001), log(max_lambda2), length.out = nlambda2))
  
  lambda12 <- expand.grid(lambda1 = lambda1s, lambda2 = lambda2s)
  
  ## Fails for n = 1000, sim3, seed = 303 at obs 49 & lam_num 2
  ### This takes an annoyingly long time
  fgl_ests <- lapply(1:nrow(lambda12), function(lam_num) {
    lam1 <- lambda12[lam_num,"lambda1"]
    lam2 <- lambda12[lam_num,"lambda2"]
    fgl_res <- JGL_cov(list(realSigX, realSigY), penalty = "fused", lambda1 = lam1, lambda2 = lam2, return.whole.theta = TRUE)
    return(list(fgl = fgl_res, lam_num = lam_num, lambda1 = lam1, lambda2 = lam2))
  })
  
  aic_fgls <- sapply(fgl_ests, function(res) aic_FGL(list(realSigX, realSigY), res$fgl[["theta"]], c(nx,ny)))
  fgl_aic_num <- max(which(aic_fgls == min(aic_fgls)))
  FGLDiff <- fgl_ests[[fgl_aic_num]]$fgl$theta[[2]] - fgl_ests[[fgl_aic_num]]$fgl$theta[[1]]
  FGLDiff <- FGLDiff[1:p,1:p] + complex(imaginary = FGLDiff[(p+1):(2*p), 1:p])
  
  ## If it's invertible
  if(!inherits(try(solve(realSigY) - solve(realSigX), silent = TRUE), "try-error")){ 
    
    inv_diff <- solve(realSigY) - solve(realSigX)
    inv_diff_max <- max(abs(inv_diff))
    inv_diff_min <- min(abs(inv_diff))
    thr_diff_lam_seq <- exp(seq(log(inv_diff_min), log(inv_diff_max), length.out = nlambdas))
    
    ## Naive Hardthrehsold
    naiveHard <- lapply(thr_diff_lam_seq, function(thr){
      diff_thr <- threshold(inv_diff, thr, hard = TRUE)
      bic <- calc_eBIC(realSigX, realSigY, diff_thr, nx = nx, ny = ny)
      return(list(bic = bic, est = diff_thr))
    })
    hard_opt_num <- which.min(sapply(naiveHard, "[[", "bic"))
    naiveHardDiff <- naiveHard[[hard_opt_num]]$est
    naiveHardDiff = naiveHardDiff[1:p,1:p] + complex(imaginary = naiveHardDiff[(p+1):(2*p), 1:p])
    
    ## Naive Softthreshold
    naiveSoft <- lapply(thr_diff_lam_seq, function(thr){
      diff_thr <- threshold(inv_diff, thr, hard = FALSE)
      bic <- calc_eBIC(realSigX, realSigY, diff_thr, nx = nx, ny = ny)
      return(list(bic = bic, est = diff_thr))
    })
    soft_opt_num <- which.min(sapply(naiveSoft, "[[", "bic"))
    naiveSoftDiff <- naiveSoft[[soft_opt_num]]$est
    naiveSoftDiff = naiveSoftDiff[1:p,1:p] + complex(imaginary = naiveSoftDiff[(p+1):(2*p), 1:p])
  } else{
    ## Not invertible
    hard_opt_num <- NULL
    naiveHardDiff <- NULL
    soft_opt_num <- NULL
    naiveSoftDiff <- NULL
  }
  
  return(list(dtl_est = opt_est[1:p,1:p] + complex(imaginary = opt_est[(p+1):(2*p), 1:p]), 
              lambdas = lam_seq, lam_range = c(min(lam_seq), max(lam_seq)), 
              lam_opt = lam_seq[lam_opt_num], lam_opt_num = lam_opt_num,
              eBICXnum = eBICXnum, eBICYnum = eBICYnum,
              naiveX = naiveX, naiveY = naiveY, naiveGLASSODiff = naiveGLASSODiff,
              FGLDiff = FGLDiff, fgl_res = fgl_ests[fgl_aic_num], fgl_aic_num = fgl_aic_num,
              realSigY = realSigY, realSigX = realSigX,
              hard_opt_num = hard_opt_num, naiveHardDiff = naiveHardDiff, 
              soft_opt_num = soft_opt_num, naiveSoftDiff = naiveSoftDiff))
}

# Read in EEG MATLAB data 
## If two sessions read in both
## If one session split in half
readEEG <- function(subj_dir, ds_rate, analyze_within_session, use_seconds){
  
  ## Get bad channels
  bad_ch_files <- list.files(subj_dir, "*.tsv",recursive = TRUE, full.names = TRUE)
  bad_chs <- unique(do.call(c, lapply(bad_ch_files, function(f) which(data.table::fread(f)$status == "bad"))))
  
  ## Sampling frequency - read in and adjust based on desired downsample rate
  srate <- unique(data.table::fread(bad_ch_files[1])$sampling_frequency)
  stopifnot( (length(srate) == 1) & (srate == 1024) )
  srate <- srate/ds_rate
  
  data_files <- list.files(subj_dir, "*.set",recursive = TRUE, full.names = TRUE)
  ## Read in t1 data, rearrange, drop bad channels
  t1_file <- data_files[grepl("ses-t1", data_files)]
  t1dat <- readMat(t1_file)$data
  t1dat <- t(matrix(t1dat, nrow = dim(t1dat)[1], ncol = prod(dim(t1dat)[2:3])))
  if(length(bad_chs)!=0){
    t1dat <- t1dat[,-bad_chs] 
  }
  t1dat <- t1dat[seq(1,nrow(t1dat), ds_rate), ]
  if(analyze_within_session){
    ## Use first use_seconds, take a use_seconds break, then get next use_seconds
    Xdat <- t1dat[1:(use_seconds*srate),]
    Ydat <- t1dat[(use_seconds*srate*2 + 1):(use_seconds*srate*3),]
  } else{
    ## Read in t2 data, rearrange, drop bad channels
    t2_file <- data_files[grepl("ses-t2", data_files)]
    t2dat <- readMat(t2_file)$data
    t2dat <- t(matrix(t2dat, nrow = dim(t2dat)[1], ncol = prod(dim(t2dat)[2:3])))
    t2dat <- t2dat[,-bad_chs]
    t2dat <- t2dat[seq(1,nrow(t2dat), ds_rate), ]
    
    Xdat <- t1dat
    Ydat <- t2dat
  }
  
  
  Xdat <- ts(scale(Xdat), frequency = srate)
  Ydat <- ts(scale(Ydat), frequency = srate)
  return(list(Xdat = Xdat, Ydat = Ydat))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#    Finding subjects corresponding to each analysis   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

subjs_w_s1 <- unique(gsub(".*(sub-\\d{3}).*", "\\1", list.dirs()[grepl("ses-t1",list.dirs())]))
subjs_w_s2 <- unique(gsub(".*(sub-\\d{3}).*", "\\1", list.dirs()[grepl("ses-t2",list.dirs())]))

subj_across_session <- intersect(subjs_w_s1, subjs_w_s2)
subj_within_session <- setdiff(subjs_w_s1, subjs_w_s2)
all_subjs <- sort(c(subj_across_session, subj_within_session))

ds_rate = 2
smooth_window_fn <- function(n){ceiling(n^(2/3))}
freqs <-  c(4,5,6,7,8, # Theta
            12,16,20,24,28, # Beta
            30,40,50,60,70, # Gamma
            80,95,110,125,140,150) # High-gamma #Freq in Hz
use_seconds = 60
options(warn = 1)

#~~~~~~~~~~~~~~~~~~~~~~~#
#    Analyze subjects   #
#~~~~~~~~~~~~~~~~~~~~~~~#
# Skip subj 100 b/c only has 43*4 = 172s but we need 180s min
# Skip subj 104 b/c t2 has no data
skip_subjs <- c("sub-100", "sub-104")
if(!(subj %in% skip_subjs)){
  print(paste0("---  Start analyzing ", subj))
  
  subj_timing <- system.time({
    within_session <- (subj %in% subj_within_session)
    print(paste0("------- Use within session analysis: ", within_session))
    XYdat <- readEEG(paste0(getwd(),"/",subj), ds_rate = ds_rate, analyze_within_session = within_session, use_seconds = use_seconds)
    Xdat <- XYdat$Xdat
    Ydat <- XYdat$Ydat
    
    specX <- astsa::mvspec(Xdat, kernel = kernel("daniell", smooth_window_fn(nrow(Xdat))), taper = 0, pad = 0, demean = FALSE, detrend = FALSE, plot = FALSE)
    specY <- astsa::mvspec(Ydat, kernel = kernel("daniell", smooth_window_fn(nrow(Ydat))), taper = 0, pad = 0, demean = FALSE, detrend = FALSE, plot = FALSE)
    
    diffnetests <- lapply(1:length(freqs), function(i){
      print(paste0("------- Frequency: ", i))
      freqi <- freqs[i]
      X = pickfxx(specX, freqi)
      Y = pickfxx(specY, freqi)
      
      X_Y_diff <- runDtl(X$fxx, Y$fxx, X$orig.n, Y$orig.n, nlambdas = 20)
      return(list(X_Y_diff = X_Y_diff))
    })
    diffnetests <- setNames(diffnetests, freqs)
    
    saveRDS(diffnetests, paste0(res_folder, "/", subj))
  })
  print(paste0("---  Finished analyzing ", subj, ". Time taken: ", round(subj_timing["elapsed"],1), "s"))
}

