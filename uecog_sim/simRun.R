## If run in CMD line - run in oslocation folder
args <- commandArgs(TRUE)
print(args)

for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}

library(Difdtl)
library(data.table)

res_folder <- paste0(getwd(), "/results")

smooth_window_fn <- function(n){ceiling(n^(2/3))}
freqs <-  c(4,5,6,7,8, # Theta
            12,16,20,24,28, # Beta
            30,40,50,60,70, # Gamma
            80,95,110,125,140,150) # High-gamma #Freq in Hz

## Inputs
blockX_nms <- c("RecBlock1", "RecBlock2", "RecBlock3", "RecBlock4", "RecBlock5",
                "RecBlock1", "RecBlock2", "RecBlock3", "RecBlock4", "RecBlock5")

blockY_nms <- c("RecBlock2", "RecBlock3", "RecBlock4", "RecBlock5", "RecBlock6",
                "CondBlock1", "CondBlock2", "CondBlock3", "CondBlock4", "CondBlock5")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  Edited astsa::mvspec function to take specific freqs of interest   #
#  should make it much faster for large datasets                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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
    n_add_zeros <- nextn(dim(pgrams_to_smooth)[1]) - dim(pgrams_to_smooth)[1] # Number of zeros to add. Do this zero-padding b/c will make MUCH faster (up to 20x)
    for (i in 1:ncol(x)){
      for (j in 1:ncol(x)){ 
        #print(paste0("i = ", i, " j = ", j))
        i_j_to_smooth <- c(pgrams_to_smooth[,i,j], rep(0, n_add_zeros))
        pgrams_final[,i,j] <- kernapply(i_j_to_smooth, kernel, circular = FALSE)[want_idxs]
      }
    }
    Lh = 1/sum(kernel[-kernel$m:kernel$m]^2)
  }
  else {
    Lh <- 1
  }
  
  df <- 2 * Lh
  df <- df * (N0/N)
  bandwidth <- Lh * xfreq/N

  spg.out <- list(freq = freq, ff_nums = ff_nums,
                  kernel = kernel, df = df, bandwidth = bandwidth, 
                  fxx = base::aperm(pgrams_final, c(2, 3, 1))/(2*pi), ## Rescale by 2pi - note astsa::mvspec does not scale by 2pi
                  Lh = Lh, n.used = N, orig.n = N0, 
                  series = series, snames = colnames(x), method = ifelse(!is.null(kernel), "Smoothed Periodogram", "Raw Periodogram"), taper = taper, 
                  pad = pad, detrend = detrend, demean = demean)
  class(spg.out) <- "spec"
  return(spg.out)
}


# find ff closest to desired frequency (freqi) using logic from pickfxx
findff <- function(x, freqs){
  stopifnot(all(!is.na(x)))
  fast <- TRUE
  N <- nrow(x)
  NewN <- if (fast) nextn(N) else N
  N <- NewN
  Nspec <- floor(N/2)
  xfreq <- frequency(x)
  ff_freqs <- seq(from = 0, by = xfreq/N, length = Nspec+1)
  ff_freq_use <- sapply(freqs, function(freqi) {if(any(abs(freqi - ff_freqs)<1e-10)){which(abs(freqi - ff_freqs)<1e-10)} else{ max(which(freqi > ff_freqs))}})
  return(ff_freq_use-1) #subtract 1 b/c we add 1 in mvspec_ff: freq <- seq(from = 0, by = xfreq/N, length = Nspec+1)[ff_nums+1]s
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


runDtl <- function(sigX, sigY, nx, ny, bic_norm, nlambdas = 20){
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
  
  return(list(dtl_est = opt_est[1:p,1:p] + complex(imaginary = opt_est[(p+1):(2*p), 1:p]), 
              lambdas = lam_seq, lam_range = c(min(lam_seq), max(lam_seq)), 
              lam_opt = lam_seq[lam_opt_num], lam_opt_num = lam_opt_num,
              realSigY = realSigY, realSigX = realSigX
              ))
}


readuEcoG <- function(session_key, blockX_nm, blockY_nm){
  
  ## Read in data
  blockX_data <- readRDS(paste0("./data/", session_key, paste0("/", blockX_nm, ".rds")))
  blockY_data <- readRDS(paste0("./data/", session_key, paste0("/", blockY_nm, ".rds")))

  params <- readRDS(paste0("./data/", session_key, "/params.rds"))
  freq <- params$fs / params$ds
  
  # only use 2nd half of Y block if Cond block for computational reasons
  if(grepl("Cond",blockY_nm)){
    blockY_data <- blockY_data[floor(nrow(blockY_data)/2):nrow(blockY_data), ]
  }
  
  Xdat <- ts(scale(blockX_data), frequency = freq)
  Ydat <- ts(scale(blockY_data), frequency = freq)
  
  return(list(Xdat = Xdat, Ydat = Ydat))
}


diffnetests_all <- Map(function(blockX_nm, blockY_nm){
  print(paste0("---  Start analyzing ", blockX_nm, " -> ", blockY_nm))
  
  timing <- system.time({
    XYdat <- readuEcoG(session_key = session_key, blockX_nm = blockX_nm, blockY_nm = blockY_nm)
    
    ff_freqsX <- findff(XYdat$Xdat, freqs)
    specX_time <- system.time({specX <- mvspec_ff(XYdat$Xdat, ff_freqsX, kernel = kernel("daniell", smooth_window_fn(nrow(XYdat$Xdat))), taper = 0, pad = 0, demean = FALSE, detrend = FALSE)})
    specX_freqs <- lapply(freqs, function(freqi) pickfxx(specX, freqi))
    
    ff_freqsY <- findff(XYdat$Ydat, freqs)
    specY_time <- system.time({specY <- mvspec_ff(XYdat$Ydat, ff_freqsY, kernel = kernel("daniell", smooth_window_fn(nrow(XYdat$Ydat))), taper = 0, pad = 0, demean = FALSE, detrend = FALSE)})
    specY_freqs <- lapply(freqs, function(freqi) pickfxx(specY, freqi))
    
    nx <- nrow(XYdat$Xdat)
    ny <- nrow(XYdat$Ydat)
    
    diffnetests <- lapply(1:length(freqs), function(i){
      print(paste0("------- Frequency: ", i))
      freqi <- freqs[i]
      X = specX_freqs[[i]]
      Y = specY_freqs[[i]]
      
      X_Y_diff <- runDtl(X$fxx, Y$fxx, X$orig.n, Y$orig.n, nlambdas = 20, bic_norm = bic_norm)
      colnames(X_Y_diff$dtl_est) <- rownames(X_Y_diff$dtl_est) <- colnames(XYdat$Xdat)
      return(list(X_Y_diff = X_Y_diff))
    })
    
    diffnetests <- setNames(diffnetests, freqs)
  })
  
  print(paste0("---  Finished analyzing ", blockX_nm, " -> ", blockY_nm, ". Time taken: ", round(timing["elapsed"],1), "s"))
  return(diffnetests)
}, blockX_nms, blockY_nms)

res <- setNames(diffnetests_all, paste0(blockX_nms, "_", blockY_nms))

saveRDS(res, paste0(res_folder, "/", session_key))
