#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### Calculate true spectral density from known VAR(1) process
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
true_spec_dens <- function(lam, A, sig_eps){
  scriptA_inv = solve( diag(1, dim(A)[1]) - A*exp(-(1i)*lam) )
  return(1/(2*pi) * (scriptA_inv) %*% sig_eps %*% Conj(t(scriptA_inv)))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### eBIC for GLASSO
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

eBIC_glasso <- function(theta, S, n, gamma){
  p <- nrow(S)
  E <- (abs(theta)>1e-10)
  E <- E[upper.tri(E)]
  cardE <- sum(E)
  lik <- (n/2)*(log(det(theta)) - sum(diag(theta %*% S)))
  return((-2*lik + cardE*log(n) + 4*cardE*gamma*log(p))) 
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### AIC for FGL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## S_l = list of S's
## theta_l = list of thetas
## n = vector of ns for each corresponding theta
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### Simulate VAR(1) processes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

genAMats <- function(p, minA, maxA, propNonZero_full, propNonZero_add, #prop of A1 that is nonzero
                     addMatDim){
  base_mat = matrix(0,nrow = p, ncol = p)
  samp_idxs <- sample(1:length(base_mat), floor(propNonZero_full*length(base_mat)))
  base_mat[samp_idxs] <- ifelse(runif(length(samp_idxs)) <= 0.5, runif(length(samp_idxs), minA, maxA), runif(length(samp_idxs), -maxA, -minA))
  
  addMat <- matrix(0, nrow = addMatDim, ncol = addMatDim)
  samp_idxs2 <- sample(1:length(addMat), floor(propNonZero_add*length(addMat)))
  addMat[samp_idxs2] <- ifelse(runif(length(samp_idxs2)) <= 0.5, runif(length(samp_idxs2), minA, maxA), runif(length(samp_idxs2), -maxA, -minA))
  
  A1 <- as.matrix(Matrix::bdiag(base_mat, addMat))
  A2 <- as.matrix(Matrix::bdiag(base_mat, addMat*(-1)))
  
  max_evals <- max(abs(eigen(A1)$values), abs(eigen(A2)$values))
  if(max_evals > 1){
    print(paste0("Max eval is: ", max_evals))
    A1 <- A1/(max_evals*1.2)
    A2 <- A2/(max_evals*1.2)
    stopifnot(max(abs(eigen(A1)$values), abs(eigen(A2)$values)) < 1)
  }
  
  return(list(A1=A1, A2=A2))
}

genAMatsSunEtAl <- function(p){
  stopifnot( (p %% 3) == 0)
  mat1 <- matrix(c(0.5,0,0,0.9,0.5,0,0,0.9,0.5), nrow = 3)
  mat2 <- -mat1
  A1 <- as.matrix(Matrix::bdiag(rep(list(mat1), times = p/3)))
  A2 <- as.matrix(Matrix::bdiag(rep(list(mat1, mat2), times = c(p/3-1, 1))))
  return(list(A1=A1, A2=A2, p =p))
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Simulate real VAR(1) from a transition matrix   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
sim_VAR1 <- function(p, n, burn_in,
                     mu_eps, sig_eps, #Mean/Cov matrix of error terms
                     errors_type, A, freqs #freqs are fourier frequency numbers divided by n
                     ){
  stopifnot(max(abs(eigen(A)$values))<1)

  #Generate data
  X_init <- mu_eps #Initialize with mean of process. This will almost always be 0 because always assume its standardized
  if(errors_type == "normal"){eps = MASS::mvrnorm(n + burn_in, mu = mu_eps, Sigma = sig_eps)} 
  else{stop("Must specify normal errors")}
  
  X_sim <- matrix(NA, nrow = n + burn_in, ncol = p)
  X_sim[1,] <- X_init
  for(i in 2:(n+burn_in)){
    X_sim[i,] <- A %*% X_sim[i-1,] + eps[i,]
  }
  X_obs <- X_sim[(burn_in + 1):(n+burn_in),]
  
  stopifnot(n%%2 ==0) #n must be even just for ease of use
  #freqs = c(1:(n/2), n)/n
  true_specs <- lapply(freqs, function(freq){
    lam = freq*2*pi
    true_spec <- true_spec_dens(lam=lam, A=A, sig_eps = sig_eps)
    true_spec_exp <- rbind(cbind(Re(true_spec), -Im(true_spec)), cbind(Im(true_spec), Re(true_spec)))
    return(list(spec = true_spec, spec_exp = true_spec_exp, lam = lam, freq = freq))
  })
  
  #Store parameters just so we have for later
  parms = list(p = p, n = n, burn_in = burn_in, mu_eps = mu_eps, sig_eps = sig_eps, errors_tpe = errors_type)
  
  return(list(A = A, X = X_obs, true_specs = true_specs, freqs = freqs, parms = parms))
}


mvspec_ff <- function(x, ff_nums, #ff_nums must be between 0 and floor(N/2)
                      spans = NULL, kernel = NULL, taper = 0, pad = 0,
                      fast = TRUE, demean = FALSE, detrend = FALSE, 
                      log = "n", type = NULL, na.action = na.fail){
  na.fail = stats::na.fail
  as.ts = stats::as.ts
  frequency = stats::frequency
  is.tskernel = stats::is.tskernel
  spec.taper = stats::spec.taper
  nextn = stats::nextn
  mvfft = stats::mvfft
  kernapply = stats::kernapply
  df.kernel = stats::df.kernel
  series <- deparse(substitute(x))
  x <- na.action(as.ts(x))
  xfreq <- frequency(x)
  x <- as.matrix(x)
  N <- N0 <- nrow(x)
  ff_nums <- sort(ff_nums)
  
  nser <- ncol(x)
  if (!is.null(spans)) 
    kernel <- {
      if (is.tskernel(spans)) 
        spans
      else kernel("modified.daniell", spans%/%2)
    }
  if (!is.null(kernel) && !is.tskernel(kernel)) 
    stop("must specify 'spans' or a valid kernel")
  if (detrend) {
    t <- 1:N - (N + 1)/2
    sumt2 <- N * (N^2 - 1)/12
    for (i in 1:ncol(x)) x[, i] <- x[, i] - mean(x[, i]) - 
      sum(x[, i] * t) * t/sumt2
  }else if (demean) {
    x <- sweep(x, 2, colMeans(x))
  }
  x <- spec.taper(x, taper)
  u2 <- (1 - (5/8) * taper * 2)
  u4 <- (1 - (93/128) * taper * 2)
  if (pad > 0) {
    x <- rbind(x, matrix(0, nrow = N * pad, ncol = ncol(x)))
    N <- nrow(x)
  }
  NewN <- if (fast) nextn(N) else N
  x <- rbind(x, matrix(0, nrow = (NewN - N), ncol = ncol(x)))
  N <- nrow(x)
  Nspec <- floor(N/2)
  stopifnot(all(ff_nums %in% 0:Nspec))
  freq <- seq(from = 0, by = xfreq/N, length = Nspec+1)[ff_nums+1]
  
  xfft <- mvfft(x)
  
  ff_nums2 <- ff_nums + 1 # add 1 for compatibility with astsa::mvspec
  kernlen <- kernel$m
  idxs <- lapply(ff_nums2, function(ff_num){
    idxs <- (ff_num - kernlen):(ff_num + kernlen)
  })
  idxs2 <- sort(unique(do.call(c, idxs)))
  idxs_all <- ifelse(idxs2 <= 0, N+idxs2, ifelse(idxs2 > N,idxs2 - N, idxs2)) #Turns out ifelse is faster than an sapply
  
    pgrams_to_smooth <- array(NA, dim = c(length(idxs_all), ncol(x), ncol(x)))
    pgrams_final <- array(NA, dim = c(length(ff_nums2), ncol(x), ncol(x)))
    xfft_subs <- xfft[idxs_all,]
    for (i in 1:ncol(x)) {
      for (j in 1:ncol(x)) {
        pgrams_to_smooth[, i, j] <- xfft_subs[, i] * Conj(xfft_subs[, j])/(N0 * xfreq)
        
        if(1 %in% idxs_all){
          pgrams_to_smooth[idxs_all == 1, i, j] <- 0.5 * (xfft[2, i] * Conj(xfft[2, j])/(N0 * xfreq) + xfft[N, i] * Conj(xfft[N, j])/(N0 * xfreq))
        }
      }
    }
    if (!is.null(kernel)) {
      want_idxs <- which(idxs_all %in% ff_nums2) - kernlen # Select only pgrams of interest
      for (i in 1:ncol(x)) for (j in 1:ncol(x)) pgrams_final[,i,j] <- kernapply(pgrams_to_smooth[,i,j], kernel, circular = FALSE)[want_idxs]
      Lh = 1/sum(kernel[-kernel$m:kernel$m]^2)
    }
    else {
      Lh <- 1
    }
  
  df <- 2 * Lh
  df <- df * (N0/N)
  bandwidth <- Lh * xfreq/N
  fxx = base::aperm(pgrams_final, c(2, 3, 1))/(2*pi) #Rescale by 2pi
  fxx_to_smooth = base::aperm(pgrams_to_smooth, c(2, 3, 1))/(2*pi)
  spg.out <- list(freq = freq, ff_nums = ff_nums,
                  kernel = kernel, df = df, bandwidth = bandwidth, 
                  fxx = fxx, fxx_to_smooth = fxx_to_smooth,
                  Lh = Lh, n.used = N, orig.n = N0, 
                  series = series, snames = colnames(x), method = ifelse(!is.null(kernel), 
                                                                         "Smoothed Periodogram", "Raw Periodogram"), taper = taper, 
                  pad = pad, detrend = detrend, demean = demean)
  class(spg.out) <- "spec"
  return(spg.out)
}

threshold <- function(X, thr, hard){
  if(hard){
    X[abs(X) < thr] <- 0
  } else{
    ifelse(abs(X) < thr, 0, ifelse(X > 0, X - thr, X + thr))
  }
  return(X)
}


calcBIC <- function(sigX, sigY, delta, nx, ny, norm, edgeMag = 1e-8){ #BIC for sigY^{-1} - sigX^{-1}
  mat <- (1/2)*((sigX %*% delta %*% sigY) + (sigY %*% delta %*% sigX)) - (sigX - sigY)
  mat_norm <- if(norm == "frob") {sum(mat**2)^(1/2)} 
  else if (norm == "max") {max(abs(mat))} 
  else {stop("Matrix norm: ", norm, " not implemented.")}
  bic <- (nx + ny)*mat_norm + log(nx + ny)*sum(abs(delta)>edgeMag)
  return(bic)
}

offdiags <- function(X){
  return(X[row(X) != col(X)])
}



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
  g1 <- graph.adjacency(critboth)
  cout = clusters(g1)
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

