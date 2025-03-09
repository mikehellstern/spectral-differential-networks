#------------------------------------------------------------------------------------------#
# As along as sim_functions.R and a "results" folder are available in this directory
# and the packages are installed, these experiments should run without issue
#------------------------------------------------------------------------------------------#

library(Difdtl)  # devtools::install_github("SusanYuan/Difdtl")
library(data.table)
library(igraph)
library(parallel)
ncores <- ifelse(detectCores() >= 8, 6, 1)
setwd("./simulations")
source("./sim_functions.R")

eval_sims <- function(input_list, n, idxs, use_true_diags = FALSE, seed = NULL){
  if(!is.null(seed)){set.seed(seed)}
  p <- input_list[["p"]]
  burn_in <- input_list[["burn_in"]]
  mu_eps <- input_list[["mu_eps"]]
  sig_eps <- input_list[["sig_eps"]]
  errors_type <- input_list[["errors_type"]]
  A1 <- input_list[["A1"]]
  A2 <- input_list[["A2"]]
  smooth_window_fn <- function(n) {ceiling(n^(2/3))}
  smooth_window_fn_nm <- "n_2/3"
  fnums <- -n/2 - 1 + idxs
  
  DataSim_timing <- system.time({
    cond1_data <- sim_VAR1(p = p, n = n, burn_in = burn_in, mu_eps = mu_eps, 
                           sig_eps = sig_eps, errors_type = errors_type, A = A1, freqs = fnums/n)
    cond2_data <- sim_VAR1(p = p, n = n, burn_in = burn_in, mu_eps = mu_eps, 
                           sig_eps = sig_eps, errors_type = errors_type, A = A2, freqs = fnums/n)
  })
  
  ### Estimate periodograms for cond 1 and 2
  #### Use sapply() for USE.NAMES arg. names are same as those in "smooth_window_fns"
  SmoothPeri_timing <- system.time({
    smooth_peri_cond1 <- mvspec_ff(x = cond1_data$X, ff_nums = fnums, kernel = kernel("daniell", smooth_window_fn(n)), taper = 0, pad = 0, demean = FALSE, detrend = FALSE)
    smooth_peri_cond2 <- mvspec_ff(x = cond2_data$X, ff_nums = fnums, kernel = kernel("daniell", smooth_window_fn(n)), taper = 0, pad = 0, demean = FALSE, detrend = FALSE)  
  })
  
  
  smooth_nums <- smooth_window_fn(n)
  
  AllRes_timing <- system.time({
    #### Loop over IDXs and do L1 D-trace loss
    results <- mclapply(1:length(fnums), function(idx){
      print(paste0("Idx #: ", idx))
      ## Create correct trueX/trueY
      trueX_exp <- cond1_data$true_specs[[idx]]$spec_exp
      trueY_exp <- cond2_data$true_specs[[idx]]$spec_exp
      
      cond1Spec <- smooth_peri_cond1$fxx[,,idx]
      sigCond1Hat <- rbind(cbind(Re(cond1Spec), -Im(cond1Spec)), cbind(Im(cond1Spec), Re(cond1Spec)))
      
      cond2Spec <- smooth_peri_cond2$fxx[,,idx]
      sigCond2Hat <- rbind(cbind(Re(cond2Spec), -Im(cond2Spec)), cbind(Im(cond2Spec), Re(cond2Spec)))
      
      ## Enforce symmetry for numerical reasons
      sigCond1Hat <- (sigCond1Hat + t(sigCond1Hat))/2
      sigCond2Hat <- (sigCond2Hat + t(sigCond2Hat))/2
      
      ### For ref with p = 30, n = 2000 we get around 0.01 MAD and 0.3 relative MAD
      mean_abs_devX = mean(abs(trueX_exp - sigCond1Hat)); mean_absX = mean(abs(trueX_exp))
      mean_abs_devY = mean(abs(trueY_exp - sigCond2Hat)); mean_absY = mean(abs(trueY_exp))
      mseX = mean((trueX_exp - sigCond1Hat)^2); mean_sqX = mean(trueX_exp^2)
      mseY = mean((trueY_exp - sigCond2Hat)^2); mean_sqY = mean(trueY_exp^2)
      
      true_diff <- (solve(trueY_exp) - solve(trueX_exp))
      
      
      ### D-trace tuning with BIC
      n_lambdas <- 20
      max_lam <- max(abs(sigCond1Hat - sigCond2Hat))*2
      lam_seq <- exp(seq(log(max_lam*0.001), log(max_lam), length.out = n_lambdas))
      
      Dtrace_est_timing <- system.time({
        dtl_ests <- sapply(1:n_lambdas, function(i){
          lam <- lam_seq[i]
          est <- L1_dts(SigmaX = round(sigCond1Hat,8), SigmaY = round(sigCond2Hat,8), rho = 1, lambda = lam) #Remeber L1_dts returns est of SigmaY^(-1) - SigmaX^(-1). See ?L1_dts. The example makes this clear
          return(list(est = est, lambda = lam, smooth_window_fn_nm = smooth_window_fn_nm, lam_num = i))
        }, simplify = F, USE.NAMES = T)
      })
      
      ### Naive difference method
      naive_lam_maxX <- max(abs(sigCond1Hat))
      naive_lam_maxY <- max(abs(sigCond2Hat))
      
      naive_lam_seqX <- exp(seq(log(naive_lam_maxX*0.0001), log(naive_lam_maxX), length.out = n_lambdas))
      naive_lam_seqY <- exp(seq(log(naive_lam_maxY*0.0001), log(naive_lam_maxY), length.out = n_lambdas))
      
      naiveXs <- lapply(naive_lam_seqX, function(lamX){
        thetaX <- glassoFast::glassoFast(sigCond1Hat, rho = lamX)
        ebic <- eBIC_glasso(thetaX$wi, sigCond1Hat, n, gamma = 0.5) #gamma = 0.5 based on eBIC paper - this seems to give decent performance in paper
        return(list(thetaX = thetaX, ebic = ebic))
      })
      naiveYs <- lapply(naive_lam_seqY, function(lamY){
        thetaY <- glassoFast::glassoFast(sigCond2Hat, rho = lamY)
        ebic <- eBIC_glasso(thetaY$wi, sigCond2Hat, n, gamma = 0.5) #gamma = 0.5 based on eBIC paper - this seems to give decent performance in paper
        return(list(thetaY = thetaY, ebic = ebic))
      })
      
      eBICXnum = which.min(sapply(naiveXs, "[[", "ebic"))
      eBICYnum = which.min(sapply(naiveYs, "[[", "ebic"))
      
      naiveX = naiveXs[[eBICXnum]]$thetaX$wi
      naiveY = naiveYs[[eBICYnum]]$thetaY$wi
      glasso_est = naiveY - naiveX
      
      ### Fused Graphical Lasso estimates
      Fgl_est_timing <- system.time({
        nlambda1 <- 2
        max_lambda1 <- max(abs(offdiags(sigCond1Hat) + offdiags(sigCond2Hat)))/2 ## point where should just be disconnected (atleast for this condition, for actual disconnected we'd need to satisfy other conditions too)
        lambda1s <- seq(0.01, 0.1, length.out = nlambda1)*max_lambda1
        
        nlambda2 <- n_lambdas/nlambda1
        max_lambda2 <- max(max(abs(offdiags(sigCond1Hat))) - min(lambda1s), max(abs(offdiags(sigCond2Hat))) - min(lambda1s)) ## point where just disconnected (for these)
        lambda2s <-  exp(seq(log(max_lambda2*0.0001), log(max_lambda2), length.out = nlambda2))
        
        lambda12 <- expand.grid(lambda1 = lambda1s, lambda2 = lambda2s)
        
        ## Fails for n = 1000, sim3, seed = 303 at obs 49 & lam_num 2
        ### This takes an annoyingly long time
        fgl_ests <- lapply(1:nrow(lambda12), function(lam_num) {
          lam1 <- lambda12[lam_num,"lambda1"]
          lam2 <- lambda12[lam_num,"lambda2"]
          fgl_res <- JGL_cov(list(sigCond1Hat, sigCond2Hat), penalty = "fused", lambda1 = lam1, lambda2 = lam2, return.whole.theta = TRUE)
          return(list(fgl = fgl_res, lam_num = lam_num, lambda1 = lam1, lambda2 = lam2))
        })
      })
      
      if(!inherits(try( solve(sigCond2Hat) - solve(sigCond1Hat), silent = TRUE), "try-error")){ #It's invertible
        inv_diff <- solve(sigCond2Hat) - solve(sigCond1Hat)
      } else{
        inv_diff <- NULL
      }
      
      extras_dt = data.table(n = n, p = p, smooth = list(smooth_nums), smooth_window_fns = list(smooth_window_fn),
                             est_usertime = Dtrace_est_timing[1], est_systemtime = Dtrace_est_timing[2], est_elapsedtime = Dtrace_est_timing[3], 
                             fnum = fnums[idx], idx = idx, n_lambdas = n_lambdas)
      
      
      return(list(extras_dt = extras_dt,
                  trueX_exp = trueX_exp,
                  trueY_exp = trueY_exp,
                  sigCond1Hat = sigCond1Hat,
                  sigCond2Hat = sigCond2Hat,
                  dtl_ests = dtl_ests,
                  glasso_est = glasso_est,
                  fgl_ests = fgl_ests,
                  inv_diff = inv_diff,
                  true_diff = true_diff
      ))
    }, mc.cores = ncores)
  })
  
  return(list(results = results, AllRes_timing = AllRes_timing, SmoothPeri_timing = SmoothPeri_timing, DataSim_timing = DataSim_timing))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Standard setup - using A matrices from Sun et al's paper    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
set.seed(1)
p=54;sim1 <- c(genAMatsSunEtAl(p),list(burn_in = 500, mu_eps = rep(0,p), sig_eps = diag(1,p), errors_type = "normal"))
## IDXs can range from 1 (corresponds to freq -n/2) to n + 1 (corresponds to freq = n/2). Since pgrams are conj symmetric and it takes a while to run Dtrace_tune_lambda
## do 100 evenly spaced IDXs from n/2 + 1 to n + 1 (freq=0 to freq=n/2)
sim1_100 <- eval_sims(sim1, 100, seq.int(from = 100/2 + 1, to = 100, length.out = 50), use_true_diags = FALSE, seed = 100) #Only do 50 for n = 100
sim1_200 <- eval_sims(sim1, 200, seq.int(from = 200/2 + 1, to = 200, length.out = 100), use_true_diags = FALSE, seed = 101)
sim1_500 <- eval_sims(sim1, 500, round(seq.int(from = 500/2 + 1, to = 500, length.out = 100),0), use_true_diags = FALSE, seed = 102 )
sim1_1000 <- eval_sims(sim1, 1000, round(seq.int(from = 1000/2 + 1, to = 1000, length.out = 100),0), use_true_diags = FALSE, seed = 103 )
sim1_2000 <- eval_sims(sim1, 2000, round(seq.int(from = 2000/2 + 1, to = 2000, length.out = 100),0), use_true_diags = FALSE, seed = 104 )

saveRDS(sim1_100, "./results/sim1_100")
saveRDS(sim1_200, "./results/sim1_200")
saveRDS(sim1_500, "./results/sim1_500")
saveRDS(sim1_1000, "./results/sim1_1000")
saveRDS(sim1_2000, "./results/sim1_2000")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Setup V2 - random sparse A matrices    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
set.seed(2)
p=54;sim2 <- c(genAMats(p-3, minA = 0.2, maxA = 0.5, propNonZero_full = 0.4, propNonZero_add = 0.6, addMatDim = 3), 
               list(burn_in = 500, mu_eps = rep(0,p), sig_eps = diag(1,p), errors_type = "normal", p =p))
sim2_100 <- eval_sims(sim2, 100, seq.int(from = 100/2 + 1, to = 100, length.out = 50), use_true_diags = FALSE) #Only do 50 for n = 100
sim2_200 <- eval_sims(sim2, 200, seq.int(from = 200/2 + 1, to = 200, length.out = 100), use_true_diags = FALSE)
sim2_500 <- eval_sims(sim2, 500, round(seq.int(from = 500/2 + 1, to = 500, length.out = 100),0), use_true_diags = FALSE )
sim2_1000 <- eval_sims(sim2, 1000, round(seq.int(from = 1000/2 + 1, to = 1000, length.out = 100),0), use_true_diags = FALSE )
sim2_2000 <- eval_sims(sim2, 2000, round(seq.int(from = 2000/2 + 1, to = 2000, length.out = 100),0), use_true_diags = FALSE )

saveRDS(sim2_100, "./results/sim2_100")
saveRDS(sim2_200, "./results/sim2_200")
saveRDS(sim2_500, "./results/sim2_500")
saveRDS(sim2_1000, "./results/sim2_1000")
saveRDS(sim2_2000, "./results/sim2_2000")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Setup V3 - random even more sparse A matrices    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
set.seed(3)
p=54;sim3 <- c(genAMats(p-3, minA = 0.4, maxA = 0.8, propNonZero_full = 0.05, propNonZero_add = 0.6, addMatDim = 3),
               list(burn_in = 500, mu_eps = rep(0,p), sig_eps = diag(1,p), errors_type = "normal", p =p))

sim3_100 <- eval_sims(sim3, 100, seq.int(from = 100/2 + 1, to = 100, length.out = 50), use_true_diags = FALSE) #Only do 50 for n = 100
sim3_200 <- eval_sims(sim3, 200, seq.int(from = 200/2 + 1, to = 200, length.out = 100), use_true_diags = FALSE)
sim3_500 <- eval_sims(sim3, 500, round(seq.int(from = 500/2 + 1, to = 500, length.out = 100),0), use_true_diags = FALSE )
sim3_1000 <- eval_sims(sim3, 1000, round(seq.int(from = 1000/2 + 1, to = 1000, length.out = 100),0), use_true_diags = FALSE )
sim3_2000 <- eval_sims(sim3, 2000, round(seq.int(from = 2000/2 + 1, to = 2000, length.out = 100),0), use_true_diags = FALSE )

saveRDS(sim3_100, "./results/sim3_100")
saveRDS(sim3_200, "./results/sim3_200")
saveRDS(sim3_500, "./results/sim3_500")
saveRDS(sim3_1000, "./results/sim3_1000")
saveRDS(sim3_2000, "./results/sim3_2000")
