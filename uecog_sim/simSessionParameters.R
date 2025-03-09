setwd(paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/data"))

# Simulate A matrices
sim_As <- function(nchs, stimto, stimfrom, delay, nrecbl = 6, nstimbl = 5){
  
  As <- vector("list", nrecbl + nstimbl)
  baseA <- matrix(rnorm(nchs*nchs), nchs)
  baseA <- make_A_stationary(baseA)
  
  stopifnot(delay %in% c("10ms", "100ms", "control"))
  
  stimeffectMult <- ifelse(delay == "10ms", 2, 0.5)
  
  if(is.null(stimto) & is.null(stimfrom)){
    As <- rep(list(baseA), nrecbl + nstimbl)
  } else{
    As[[1]] <- baseA
    
    # generate effect of stim. Info flows from stimFrom -> stimTo thus we modify stimA[stimto, stimfrom] entry
    stimeffect <- max(abs(baseA))*stimeffectMult
    stimA <- baseA
    stimA[stimto, stimfrom] <- baseA[stimto, stimfrom] + stimeffect
    stimA <- make_A_stationary(stimA)
    
    # effect of stim diminishes into rec. This is how much of that effect remains on average over all sessions
    stimdim <- 1/3
    # this is how much each addtl session increases the stim effect. This effect is cumulative
    deltadim <- 1/7
    
    for (i in 2:(nrecbl + nstimbl)){
      if( (i %% 2) == 0){ # even blocks = stim
        As[[i]] <- stimA  
      } else{ # odd blocks = rec
        # If nrecbl = 6, then resideffect in recblock1, 2, ... is:  (stimdim - 2delta, stimdim - delta, stimdim, stimdim + delta)
        # E.g. for nrecbl = 6, prop of stim effect is: (stimdim + deltadim * (1:5 - nrecbl/2))
        resideffect <- stimeffect * (stimdim + deltadim * ((i %/% 2) - nrecbl/2))
        recAtmp <- baseA
        recAtmp[stimto, stimfrom] <- baseA[stimto, stimfrom] + resideffect
        As[[i]] <- recAtmp
      }
    }
  }
  
  blnames <- paste0(c("RecBlock", "CondBlock"), cumsum(rep_len(c(1,0), nrecbl + nstimbl)))
  
  return(setNames(As, blnames))
}

# Take A and make it so max eigenval is < 1 to ensure stationarity
make_A_stationary <- function(A){
  max_eval <- max(abs(eigen(A)$values))
  if(max_eval > 1){
    A <- A/(1.2*max_eval)
  }
  stopifnot(max(abs(eigen(A)$values)) < 1)
  return(A)
}

# Simulate Session Parameters
sim_one_session_params <- function(name, nchs, rangebad, delay, fs, ds){
  
  nbadchs <- sample(rangebad[1]:rangebad[2], 1)
  goodchs <- sort(sample(1:nchs, nchs - nbadchs))
  
  # simulate As
  if(delay == "control"){
    # simulate As without stim
    As <- sim_As(length(goodchs), stimto = NULL, stimfrom = NULL, delay = "control")
    stimto <- NULL
    stimfrom <- NULL
  } else{
    stimto_stimfrom <- sample(goodchs, 2)
    stimto <- which(goodchs == stimto_stimfrom[1])
    stimfrom <- which(goodchs == stimto_stimfrom[2])
    if(stimto == stimfrom){stop("stimto and stimfrom are same channel")}
    # simulate As with stim
    As <- sim_As(nchs = length(goodchs), stimto = stimto, stimfrom = stimfrom, delay = delay)
    
  }
  
  ns <- rep(c(5*60*fs/ds, 10*60*fs/ds), times = 5)
  ns <- c(ns, 5*60*fs/ds)

  errcov <- diag(1, length(goodchs))
  
  params <- list(As = As, stimto = stimto, stimfrom = stimfrom, delay = delay, ns = ns, errcov = errcov, fs = fs, ds = ds)
  
  # save
  savefolder <- paste0("./", name)
  if(!dir.exists(savefolder)){
    dir.create(savefolder)
  } 
  saveRDS(params, paste0(savefolder, "/params.rds"))
  return(params)
}


# ------ Experimental sessions
experiment_sessions <- 1:32
nchannels = 96
bad_ch_range = c(2,40)
fs <- 1000
ds <- 2
delays <- c(rep("10ms", 23), rep("100ms", 9))

exp_params_l <- vector("list", length(experiment_sessions))

for (sessioni in experiment_sessions){
  seshnm <- paste0("session", sessioni)
  exp_params_l[[sessioni]] <- sim_one_session_params(name = seshnm, nchs = nchannels, rangebad = bad_ch_range, delay = delays[sessioni], fs = fs, ds = ds)
}


# ------ Control sessions
control_sessions <- 33:36

ctrl_params_l <- vector("list", length(experiment_sessions))

for (ctrl_sessioni in 1:length(control_sessions)){
  seshnm <- paste0("session", control_sessions[ctrl_sessioni])
  ctrl_params_l[[ctrl_sessioni]] <- sim_one_session_params(name = seshnm, nchs = nchannels, rangebad = bad_ch_range, delay = "control", fs = fs, ds = ds)
}
# 
# 
# i <- 31
# test_stimto <- exp_params_l[[i]]$stimto
# test_stimfrom <- exp_params_l[[i]]$stimfrom
# test_Astim <- exp_params_l[[i]]$As[[2]]
# test_Arec2 <- exp_params_l[[i]]$As[[1]]
# 
# test_Astim[test_stimto, test_stimfrom]
# test_Arec2[test_stimto, test_stimfrom]
