#------------------------------------------#
# Load packages and set wd
#------------------------------------------#
library(knitr)
library(data.table)
library(kableExtra)
library(ggplot2)
library(tikzDevice)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("./sim_functions.R")

#------------------------------------------#
# Start of helper functions
#------------------------------------------#

calc_eBIC <- function(sigX, sigY, delta, nx, ny, isComplex, edgeMag = 1e-8){
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

computeErrorMetrics <- function(est, true, edgeMag = 1e-8){ ## edgeMag is Magnitude used to define edge. If < this its not considered an edge
  trueEdges <- (abs(true) > edgeMag)
  estEdges <- (abs(est) > edgeMag)
  tp = sum(estEdges & trueEdges); fp = sum(estEdges & !trueEdges); tn = sum(!estEdges & !trueEdges); fn = sum(!estEdges & trueEdges)
  tpr = tp/(tp + fn); ppv = tp/(tp + fp)
  if((tp + fp)==0) ppv = 0
  
  RMSEedge = sqrt(mean((est[trueEdges] - true[trueEdges])^2))
  rRMSEedge = RMSEedge/sqrt(mean((true[trueEdges])^2)) 
  RMSEnonEdge = sqrt(mean((est[!trueEdges] - true[!trueEdges])^2))
  RMSEall = sqrt(mean((est - true)^2))
  rRMSEall = RMSEall/sqrt(mean((true)^2)) 
  acc = (tp + tn)/(prod(dim(true)))
  return(data.table(n_true_edge = sum(trueEdges), n_est_edge = sum(estEdges), tpr = tpr, ppv = ppv,
                    RMSEedge = RMSEedge, rRMSEedge = rRMSEedge, RMSEnonEdge = RMSEnonEdge, RMSEall = RMSEall, rRMSEall = rRMSEall, acc = acc))
}

summResultsOneSim <- function(sim_nm){
  sim_res <- readRDS(paste0("./results/", sim_nm))
  error_res <- lapply(sim_res$results, function(sim_res_i){
    n <- sim_res_i$extras_dt$n
    setnames(sim_res_i$extras_dt, "smooth", "smooth_n")
    true_diff <- sim_res_i[["true_diff"]]
    
    ## Error metrics for DTL
    sigmaX <- sim_res_i[["sigCond1Hat"]]
    sigmaY <- sim_res_i[["sigCond2Hat"]]
    dtls <- sim_res_i[["dtl_ests"]]
    bics <- sapply(dtls, function(dtl_i) calc_eBIC(sigX = sigmaX, sigY = sigmaY, delta = dtl_i[["est"]], nx = n, ny = n, isComplex = TRUE))
    if(any(bics > 1e7)){stop("BIC > 1e7 detected. Please double check results.")}
    lams <- sapply(dtls, function(dtl_i) dtl_i[["lambda"]])
    lams_valid <- lams[which(bics == min(bics))]
    est_diff_dtl <- dtls[[which(lams == max(lams_valid))]]
    error_diff_dtl <- computeErrorMetrics(est_diff_dtl[["est"]], true_diff)
    error_diff_dtl[,c("n", "fnum", "smooth_n", "idx", "p") := sim_res_i$extras_dt[,c("n","fnum","smooth_n","idx","p")]]
    
    
    ## Error metrics for Naive GLASSO
    est_diff_naive <- sim_res_i[["glasso_est"]]
    if(!is.null(est_diff_naive)){
      error_diff_naive <- computeErrorMetrics(est_diff_naive, true_diff)
      error_diff_naive[,c("n", "fnum", "smooth_n", "idx", "p") := sim_res_i$extras_dt[,c("n","fnum","smooth_n","idx","p")]]  
    } else{
      error_diff_naive = NULL
    }
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ## Error metrics for Fused LASSO
    ## - note in FGL case we have two penalties lam1/lam2.
    ## - This is important b/c we can have ties and need to pick "larger" of two lambda numbers
    ##   lam_num = 19 might be lam1=0.0038;lam2 = 0.517 and lam_num = 20 might be lam1=0.038;lam2=0.517
    ##   so as lam_num increases lam1/lam2 is non-decreasing.
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    fgls <- sim_res_i[["fgl_ests"]]
    fgl_aics <- sapply(fgls, function(fgl) { aic_FGL(list(sigmaX, sigmaY), fgl$fgl$theta, c(n,n))})
    fgl_aic_num <- max(which(fgl_aics == min(fgl_aics)))
    fgl_diff <- fgls[[fgl_aic_num]]$fgl$theta[[2]] - fgls[[fgl_aic_num]]$fgl$theta[[1]]
    error_diff_fgl <- computeErrorMetrics(fgl_diff, true_diff)
    error_diff_fgl[,c("n", "fnum", "smooth_n", "idx", "p") := sim_res_i$extras_dt[,c("n","fnum","smooth_n","idx","p")]]
    error_diff_fgl[,c("tuning_type") := "fgl_aic"]
    
    
    inv_diff <- sim_res_i[["inv_diff"]]
    if(!is.null(inv_diff)){
      inv_diff_max <- max(abs(inv_diff))
      inv_diff_min <- min(abs(inv_diff[abs(inv_diff)!=0])) # covers special cases where the min value is actually 0 (we need it to not be 0)
      nlambdas = length(lams)
      thr_diff_lam_seq <- exp(seq(log(inv_diff_min), log(inv_diff_max), length.out = nlambdas))
      
      ## Error metrics for Naive Hard threshold
      naiveHard <- lapply(thr_diff_lam_seq, function(thr){
        diff_thr <- threshold(inv_diff, thr, hard = TRUE)
        bic <- calc_eBIC(sigX = sigmaX, sigY = sigmaY, delta = diff_thr, nx = n, ny = n, isComplex = TRUE)
        return(list(bic = bic, est = diff_thr))
      })
      hard_opt_num <- which.min(sapply(naiveHard, "[[", "bic"))
      naiveHardDiff <- naiveHard[[hard_opt_num]]$est
      
      stopifnot(all(dim(naiveHardDiff) == dim(true_diff)))
      
      error_diff_hard <- computeErrorMetrics(naiveHardDiff, true_diff)
      error_diff_hard[,c("n", "fnum", "smooth_n", "idx", "p") := sim_res_i$extras_dt[,c("n","fnum","smooth_n","idx","p")]]  
      
      
    } else{
      error_diff_hard = NULL
    }
    
    
    dtl_est_timing <- sim_res_i$extras_dt[,c("n", "fnum", "smooth_n", "idx", "p", "est_usertime", "est_systemtime", "est_elapsedtime")]
    
    return(list(error_diff_dtl = error_diff_dtl, 
                error_diff_naive = error_diff_naive, 
                error_diff_hard = error_diff_hard,
                error_diff_fgl = error_diff_fgl,
                dtl_est_timing = dtl_est_timing))
  })
  nms <- c("error_diff_dtl", "error_diff_naive", "error_diff_hard", "error_diff_fgl", "dtl_est_timing")
  error_res_sum <- sapply(nms, function(nm) rbindlist(lapply(error_res, "[[", nm)), USE.NAMES = TRUE, simplify = FALSE)
  return(error_res_sum)
}

summResults <- function(sim_nms){
  sims_res_list <- lapply(sim_nms, summResultsOneSim)
  nms <- c("error_diff_dtl", "error_diff_naive", "error_diff_hard", "error_diff_fgl","dtl_est_timing")
  sim_res_combined <- sapply(nms, function(nm) rbindlist(lapply(sims_res_list, "[[", nm)), USE.NAMES = TRUE, simplify = FALSE)
  return(sim_res_combined)
}

format_col <- function(x, mean_digs = 2, se_digs = 2){
  nx <- sum(!is.na(x))
  mx <- mean(x, na.rm = T)
  varx <- var(x, na.rm = T)*((nx-1)/nx)
  se <- sqrt(varx/nx)
  if(abs(mx) > 10 & abs(mx) < 100){mean_digs <- 1} else if(abs(mx) >= 100){mean_digs <- 0}
  if(se > 10 & se < 100){se_digs <- 1} else if (se>=100){se_digs <- 0}
  mean_rounded = sprintf(paste0("%.", mean_digs, "f"), mx)
  se_rounded = sprintf(paste0("%.", se_digs, "f"), se)
  return(paste0(mean_rounded, " (", se_rounded, ")"))
}

add_miss_rows <- function(x_to_add, x_ref){
  n_all <- x_ref[["n"]]
  n_to_add <- setdiff(n_all, x_to_add[["n"]])
  add_rows <- x_ref[1:length(n_to_add),]
  add_rows[] <- "-"
  add_rows$n <- n_to_add
  x_added <- rbind(add_rows, x_to_add)
  return(x_added)
}

make_tables <- function(sim_res, cols){
  cols_nms <- names(cols)
  diff_dtl_fmt <- sim_res$error_diff_dtl[,setNames(lapply(.SD, format_col), cols_nms), by = .(n), .SDcols = cols]
  diff_naive_glasso_fmt <- sim_res$error_diff_naive[,setNames(lapply(.SD, format_col), cols_nms), by = .(n), .SDcols = cols]
  diff_fgl_fmt <- sim_res$error_diff_fgl[,setNames(lapply(.SD, format_col), cols_nms), by = .(n), .SDcols = cols]
  diff_naive_hard_fmt <- sim_res$error_diff_hard[,setNames(lapply(.SD, format_col), cols_nms), by = .(n), .SDcols = cols]
  diff_naive_hard_fmt <- add_miss_rows(diff_naive_hard_fmt, diff_dtl_fmt)
  
  diff_fmt_top <- cbind(diff_dtl_fmt, 
                        rep("",nrow(diff_dtl_fmt)), diff_naive_glasso_fmt[,-c("n","# True edges"),with=FALSE])
  diff_fmt_bot <- cbind(diff_naive_hard_fmt,
                        rep("",nrow(diff_dtl_fmt)), diff_fgl_fmt[,-c("n","# True edges"),with=FALSE])
  diff_fmt_bot[["# True edges"]] <- " "
  diff_fmt_bot[["n"]] <- " "
  colnames(diff_fmt_bot) <- rep(" ", ncol(diff_fmt_bot))
  
  table_header_top <- c(2, length(cols)-1, 1, length(cols)-1)
  names(table_header_top) <- c(" ", "SDD", " ", "GLASSO difference")
  kbl_res_top <- kable(diff_fmt_top, format = "latex", align = paste0(rep("c",ncol(diff_fmt_top))), booktabs = T) %>% add_header_above(table_header_top)
  kbl_simp_top <- kable(diff_fmt_top, format = "pipe", align = paste0(rep("c",ncol(diff_fmt_top))), booktabs = T) #%>% add_header_above(table_header_top)
  
  table_header_bot <- c(2, length(cols)-1, 1, length(cols)-1)
  names(table_header_bot) <- c(" ", "Hard threhsold difference", " ", "FGL")
  kbl_res_bot <- kable(diff_fmt_bot, format = "latex", align = paste0(rep("c",ncol(diff_fmt_bot))), booktabs = T) %>% add_header_above(table_header_bot)
  kbl_simp_bot <- kable(diff_fmt_bot, format = "pipe", align = paste0(rep("c",ncol(diff_fmt_bot))), booktabs = T) #%>% add_header_above(table_header_bot)
  
  kbl_simp_dtl <- kable(diff_dtl_fmt, format = "pipe", align = paste0(rep("c", ncol(diff_dtl_fmt))))
  kbl_simp_naive_glasso <- kable(diff_naive_glasso_fmt, format = "pipe", align = paste0(rep("c", ncol(diff_naive_glasso_fmt))))
  kbl_simp_naive_hard <- kable(diff_naive_hard_fmt, format = "pipe", align = paste0(rep("c", ncol(diff_naive_hard_fmt))))
  kbl_simp_fgl <- kable(diff_fgl_fmt, format = "pipe", align = paste0(rep("c", ncol(diff_fgl_fmt))))
  
  kbl_simp_latex_dtl <- kable(diff_dtl_fmt, format = "latex", align = paste0(rep("c", ncol(diff_dtl_fmt))))
  kbl_simp_latex_naive_glasso <- kable(diff_naive_glasso_fmt, format = "latex", align = paste0(rep("c", ncol(diff_naive_glasso_fmt))))
  kbl_simp_latex_naive_hard <- kable(diff_naive_hard_fmt, format = "latex", align = paste0(rep("c", ncol(diff_naive_hard_fmt))))
  kbl_simp_latex_fgl <- kable(diff_fgl_fmt, format = "latex", align = paste0(rep("c", ncol(diff_fgl_fmt))))
  
  ## Difference tables
  diff_dtl_glasso <- sim_res$error_diff_dtl[sim_res$error_diff_naive, on = .(n, fnum, p), nomatch = 0][,c(n = list(n), lapply(cols, function(col) get(col) - get(paste0("i.",col))))]
  diff_dtl_glasso_fmt <- diff_dtl_glasso[,lapply(.SD, format_col), by = .(n), .SDcols = names(cols)]
  diff_dtl_glasso_fmt[["# True edges"]] <- NULL
  kbl_diff_dtl_glasso <- kable(diff_dtl_glasso_fmt, format = "latex", align = paste0(rep("c",ncol(diff_dtl_glasso_fmt))), booktabs = T)
  
  diff_dtl_hard <- sim_res$error_diff_dtl[sim_res$error_diff_hard, on = .(n, fnum, p), nomatch = 0][,c(n = list(n), lapply(cols, function(col) get(col) - get(paste0("i.",col))))]
  diff_dtl_hard_fmt <- diff_dtl_hard[,lapply(.SD, format_col), by = .(n), .SDcols = names(cols)]
  diff_dtl_hard_fmt[["# True edges"]] <- NULL
  diff_dtl_hard_fmt <- add_miss_rows(diff_dtl_hard_fmt, diff_dtl_glasso_fmt)
  kbl_diff_dtl_hard <- kable(diff_dtl_hard_fmt, format = "latex", align = paste0(rep("c",ncol(diff_dtl_hard_fmt))), booktabs = T)
  
  
  diff_dtl_fgl <- sim_res$error_diff_dtl[sim_res$error_diff_fgl, on = .(n, fnum, p), nomatch = 0][,c(n = list(n), lapply(cols, function(col) get(col) - get(paste0("i.",col))))]
  diff_dtl_fgl_fmt <- diff_dtl_fgl[,lapply(.SD, format_col), by = .(n), .SDcols = names(cols)]
  diff_dtl_fgl_fmt[["# True edges"]] <- NULL
  kbl_diff_dtl_fgl <- kable(diff_dtl_fgl_fmt, format = "latex", align = paste0(rep("c",ncol(diff_dtl_fgl_fmt))), booktabs = T)
  
  return(list(table_top = kbl_res_top, table_bot = kbl_res_bot, 
              table_simple_top = kbl_simp_top, table_simple_bot = kbl_simp_bot, 
              table_dtl = kbl_simp_dtl,
              table_fgl = kbl_simp_fgl,
              table_naive_glasso = kbl_simp_naive_glasso,
              table_naive_hard = kbl_simp_naive_hard,
              table_latex_dtl = kbl_simp_latex_dtl,
              table_latex_naive_glasso = kbl_simp_latex_naive_glasso,
              table_latex_naive_hard = kbl_simp_latex_naive_hard,
              table_latex_fgl = kbl_simp_latex_fgl,
              table_diff_dtl_glasso = kbl_diff_dtl_glasso,
              table_diff_dtl_hard = kbl_diff_dtl_hard,
              table_diff_dtl_fgl = kbl_diff_dtl_fgl,
              table_long = cbind(diff_fmt_top, diff_fmt_bot)))
}



theme_Publication <- function(base_size=14, base_family="serif", opaque = FALSE) {
  library(grid)
  library(ggthemes)
  theme_use <- (theme_foundation(base_size=base_size, base_family=base_family)
                + theme(plot.title = element_text(face = "bold",
                                                  size = rel(1.2), hjust = 0.5),
                        text = element_text(family = "serif"),
                        panel.background = if(opaque){element_rect(fill='transparent', color = NA)}else{element_rect(colour = NA)},
                        panel.spacing = unit(2, "lines"),
                        plot.background = if(opaque){element_rect(fill='transparent', color = NA)}else{element_rect(colour = NA)},
                        panel.border = element_rect(colour = NA),
                        axis.title = element_text( size = rel(2)),
                        axis.title.y = element_text(angle=90,vjust =2),
                        axis.title.x = element_text(vjust = -0.2),
                        axis.text = element_text(size = rel(1.75)), 
                        axis.line = element_line(colour="black"),
                        axis.ticks = element_line(),
                        panel.grid.major = element_line(colour="#f0f0f0"),
                        panel.grid.minor = element_blank(),
                        legend.key = element_rect(colour = NA),
                        legend.text = element_text(size = rel(1.3)),
                        legend.position = "bottom",
                        legend.direction = "horizontal",
                        legend.box = "vertical",
                        legend.key.size= unit(0.2, "cm"),
                        legend.margin = margin(0,0,0,0,"cm"),
                        legend.title = element_text(face="italic", size = rel(1.5)),
                        plot.margin=unit(c(10,5,5,5),"mm"),
                        #strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                        strip.background=element_blank(),
                        strip.text = element_text(size = rel(1.8))
                ))
  if(opaque){
    theme_use <- theme_use + theme(legend.background = element_rect(fill='transparent', color = NA), legend.box.background = element_rect(fill='transparent', color = NA))
  }
  return(theme_use)
}


## NOTE - This doesn't look good at all if we include Soft Threshold because it messes up scale.
##      - Looks just OK if exclude Soft Threshold, but I still think we should go with tables
make_plots_all <- function(sim_res, cols, sim_nm){

  diff_dtl_mean <- sim_res$error_diff_dtl[,c(list(type = "mean"), lapply(.SD, mean)), by = .(n), .SDcols = cols]
  diff_dtl_se <- sim_res$error_diff_dtl[,c(list(type = "se"), lapply(.SD, function(x) sqrt(var(x)/.N))), by = .(n), .SDcols = cols]
  diff_dtl <- rbind(diff_dtl_mean, diff_dtl_se)
  
  diff_naive_glasso_mean <- sim_res$error_diff_naive[,c(list(type = "mean"), lapply(.SD, mean)), by = .(n), .SDcols = cols]
  diff_naive_glasso_se <- sim_res$error_diff_naive[,c(list(type = "se"), lapply(.SD, function(x) sqrt(var(x)/.N))), by = .(n), .SDcols = cols]
  diff_naive_glasso <- rbind(diff_naive_glasso_mean, diff_naive_glasso_se)
  
  diff_fgl_mean <- sim_res$error_diff_fgl[,c(list(type = "mean"), lapply(.SD, mean)), by = .(n), .SDcols = cols]
  diff_fgl_se <- sim_res$error_diff_fgl[,c(list(type = "se"), lapply(.SD, function(x) sqrt(var(x)/.N))), by = .(n), .SDcols = cols]
  diff_fgl <- rbind(diff_fgl_mean, diff_fgl_se)
  
  diff_naive_hard_mean <- sim_res$error_diff_hard[,c(list(type = "mean"), lapply(.SD, mean)), by = .(n), .SDcols = cols]
  diff_naive_hard_se <- sim_res$error_diff_hard[,c(list(type = "se"), lapply(.SD, function(x) sqrt(var(x)/.N))), by = .(n), .SDcols = cols]
  diff_naive_hard <- rbind(diff_naive_hard_mean, diff_naive_hard_se)
  
  res_all_mean <- rbindlist(list("SDD" = diff_dtl_mean, "Naive" = diff_naive_glasso_mean, "Hard threshold" = diff_naive_hard_mean, "FGL" = diff_fgl_mean), idcol = "Method")
  res_all_se <- rbindlist(list("SDD" = diff_dtl_se, "Naive" = diff_naive_glasso_se, "Hard threshold" = diff_naive_hard_se, "FGL" = diff_fgl_se), idcol = "Method")
  res_all_mean_m <- melt.data.table(res_all_mean, id.vars = c("Method", "n"), measure.vars = cols, variable.factor = FALSE)
  res_all_se_m <- melt.data.table(res_all_se, id.vars = c("Method", "n"), measure.vars = cols, variable.factor = FALSE)
  res_all_mean_m[,c("Method") := ifelse(variable == "n_true_edge", variable, Method)][,c("variable") := ifelse(variable == "n_true_edge", "n_est_edge", variable)]
  res_all_se_m[,c("Method") := ifelse(variable == "n_true_edge", variable, Method)][,c("variable") := ifelse(variable == "n_true_edge", "n_est_edge", variable)]
  res_all_m <- res_all_mean_m[res_all_se_m ,on = c("Method", "n", "variable"), .(Method, n, variable, value, value_se=i.value)]
  res_all_m[, c("variable") := setNames(names(cols), cols)[variable]]
  
  
  col_scale <- c("SDD" = "#785ef0", "Naive" = "#648FFF", "Hard threshold" = "#FE6100", "FGL" = "#DC267F")
  res_plt <- res_all_m
  
  tikz(file = paste0("./figures/", sim_nm, "_all_metrics", ".tex"), width = 3.5, height = 4.5, standAlone = TRUE)
  
  plt <- ggplot(res_plt, aes(x = n, y = value, color = Method)) + geom_point(size = 2, position = position_dodge(width = 30)) + geom_line(linewidth = 1) + 
    geom_errorbar(aes(ymin = value - value_se, ymax = value + value_se), width = NA, linewidth=0.66, position = position_dodge(width = 30)) + 
    xlab("T") + ylab("") + 
    scale_color_manual(values = col_scale) +
    scale_y_continuous(n.breaks = 5) + 
    facet_wrap(~variable, nrow = 5, scales = "free_y")+
    theme_Publication(base_size = 6) +   theme(legend.background = element_rect(fill='transparent', color = NA), 
                                                legend.box.background = element_rect(fill='transparent', color = NA),
                                                panel.background = element_rect(fill='transparent', color = NA),
                                                plot.background = element_rect(fill='transparent', color = NA))
  print(plt)
  dev.off()
}


#------------------------------------------#
# Start of code to summarize results
#------------------------------------------#

sd_cols <- c("n_true_edge", "n_est_edge", "tpr", "ppv", "RMSEedge", "rRMSEedge", "RMSEnonEdge", "RMSEall", "rRMSEall", "auroc", "fnum", "idx")
cols1 <- c("# True edges" = "n_true_edge", "# Est edges" = "n_est_edge", "Precision" = "ppv", "Recall" = "tpr","Accuracy" = "acc", "RRMSE" = "rRMSEall")
## Sim 1
sim1_nms <- paste0("sim1_", c("100","200","500","1000","2000"))
sim1_res <- summResults(sim1_nms)
sim1_tables <- make_tables(sim1_res, cols1)
make_plots_all(sim1_res, c("Accuracy" = "acc", "RRMSE" = "rRMSEall"), sim_nm = "sim1_acc_rrmse")
make_plots_all(sim1_res, c("Precision" = "ppv", "Recall" = "tpr"), sim_nm = "sim1_prec_reca")

## Sim 2
sim2_nms <- paste0("sim2_", c("100","200","500","1000","2000"))
sim2_res <- summResults(sim2_nms)
sim2_tables <- make_tables(sim2_res, cols1)
make_plots_all(sim2_res, c("Accuracy" = "acc", "RRMSE" = "rRMSEall"), sim_nm = "sim2_acc_rrmse")
make_plots_all(sim2_res, c("Precision" = "ppv", "Recall" = "tpr"), sim_nm = "sim2_prec_reca")


## Sim 3
sim3_nms <- paste0("sim3_", c("100","200","500","1000","2000"))
sim3_res <- summResults(sim3_nms)
sim3_tables <- make_tables(sim3_res, cols1)
make_plots_all(sim3_res, c("Accuracy" = "acc", "RRMSE" = "rRMSEall"), sim_nm = "sim3_acc_rrmse" )
make_plots_all(sim3_res, c("Precision" = "ppv", "Recall" = "tpr"), sim_nm = "sim3_prec_reca" )
