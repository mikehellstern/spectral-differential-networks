setwd(paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/data"))


# Simulate Session Data
sim_uecog_data <- function(session_nm, seed, burn_in = 1000){
  
  set.seed(seed)
  params <- readRDS(paste0("./", session_nm, "/params.rds"))
  p <- dim(params$As[[1]])[1]
  As <- params$As
  freq <- params$fs / params$ds
  
  #Generate data
  mu_eps <- rep(0,p)
  
  block_data <- vector("list", length(As))
  
  # Sim data for each block & save
  for (block in 1:length(As)){
    
    blockn <- params$ns[block]
    blockA <- As[[block]]
    
    # initialize with mean of process or last obs from previous block
    if(block == 1){
      X_init <- mu_eps
      simn <- blockn + burn_in
      throwoutn <- burn_in + 1
    } else{
      X_init <- tail(block_data[[block-1]], 1)
      simn <- blockn
      throwoutn <- 1
    }
    
    X_sim <- rbind(X_init, matrix(NA, nrow = simn, ncol = p))
    eps = MASS::mvrnorm(simn+1, mu = mu_eps, Sigma = params$errcov)
    
    for(i in 2:nrow(X_sim)){
      X_sim[i,] <- blockA %*% X_sim[i-1,] + eps[i,]
    }
    
    X_obs <- X_sim[(throwoutn + 1):(nrow(X_sim)),]
    stopifnot(nrow(X_obs)==blockn)
    stopifnot(all(!is.na(X_obs)))
    block_data[[block]] <- X_obs
    
    block_nm <- names(As)[block]
    saveRDS(X_obs, paste0("./", session_nm, "/", block_nm, ".rds"))
  }
  
  saveRDS(sessionInfo(), paste0("./", session_nm, "/sessionInfo.rds"))
  
  return(block_data)
}

# ------ Simulate data
session_nms <- paste0("session", 1)
for (session_nmi in 1:length(session_nms)){
  tmp <- sim_uecog_data(session_nm = session_nms[session_nmi], seed = session_nmi)
}








