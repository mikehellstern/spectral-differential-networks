### SRM EEG analysis
#Important note: install Difdtl using devtools::install_github("SusanYuan/Difdtl"). 
#  For Mac, requires R <= 4.1.2. The C compiler changed in later versions of R and code does not compile

## If run in Rstudio
library(R.matlab)
library(igraph)
library(ggplot2)
library(glassoFast)
library(astsa)
library(data.table)
library(kableExtra)
library(tikzDevice)
res_folder <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/results")
fig_folder <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/figures")
dat_folder <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/data/derivatives/cleaned_epochs") 


#~~~~~~~~~~~~~~~~~~~~#
#     Helper fns     #
#~~~~~~~~~~~~~~~~~~~~#

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

analyzeOneFreq <- function(freqRes, freq, ch_nms, analyzeMethods = c("dtl_est", "naiveGLASSODiff", "naiveHardDiff", "FGLDiff")){
  
  sparsity <- sapply(analyzeMethods, function(method){method_res <- freqRes[[method]]
                                                      if(!is.null(method_res)){return(mean(abs(method_res)==0))}
                                                      else{return(NA)} })
  sparsity_table <- data.table(sparsity = sparsity, method = names(sparsity))
  sparsity_table$frequency = freq
  edges <- lapply(analyzeMethods, function(method){method_res <- freqRes[[method]]
                                                   if(!is.null(method_res) && (!all(abs(method_res)==0))){
                                                     undir_edges <- which(abs(method_res) != 0 & lower.tri(method_res), arr.ind = T) 
                                                     return(data.table(from = ch_nms[undir_edges[,1]], to = ch_nms[undir_edges[,2]], method = method, frequency = freq))
                                                   } else if((!is.null(method_res) && (all(abs(method_res)==0)))){
                                                      return(NULL)
                                                    }
                                                   else{return(data.table(from = NA, to = NA, method = method, frequency = freq))} })
  edges <- rbindlist(edges)
  
  return(list(sparsity = sparsity_table, edges = edges))
}

analyzeOneSub <- function(sub, bands_dt, dat_folder, res_folder){
  
  ## Get channel names - I have to do it this way w/ bad channels b/c I didn't add channel names in the actual analysis
  bad_ch_files <- list.files(paste0(dat_folder, "/", sub), "*.tsv",recursive = TRUE, full.names = TRUE)
  bad_chs <- unique(do.call(c, lapply(bad_ch_files, function(f) which(data.table::fread(f)$status == "bad"))))
  if(length(bad_chs) > 0){
    ch_nms <- data.table::fread(bad_ch_files[1])$name[-bad_chs]  
  } else{
    ch_nms <- data.table::fread(bad_ch_files[1])$name
  }
  
  ## Load in subjects results
  sub_res <- readRDS(paste0(res_folder, "/", sub))
  sub_summ <- lapply(names(sub_res), function(freq) analyzeOneFreq(sub_res[[freq]]$X_Y_diff, as.numeric(freq), ch_nms))

  ## Edges are undirected
  edges <- rbindlist(lapply(sub_summ, "[[", "edges"))
  edges$sub <- sub
  edges[bands_dt, "band" := i.band_nm, on = .(frequency <= band_max, frequency >= band_min)]
  
  sparsity <- rbindlist(lapply(sub_summ, "[[", "sparsity"))
  sparsity$sub <- sub
  sparsity[bands_dt, "band" := i.band_nm, on = .(frequency <= band_max, frequency >= band_min)]
  
  return(list(edges = edges, sparsity = sparsity))
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

order_cols_rows <- function(dt, col_reorder, row_reorder = c("theta"=1, "beta"=2, "gamma"=3, "high_gamma"=4)){
  setcolorder(dt, col_reorder)
  dt2 <- dt[, c("ord") := row_reorder[band]][order(ord)]
  dt2[,ord:=NULL]
  return(dt2)
}

make_sparsity_tables <- function(dt, dig = 2){
  method_recode <- c("dtl_est" = "SDD", "naiveGLASSODiff" = "GLASSO", "naiveHardDiff" = "Hard threshold", "FGLDiff" = "FGL")
  
  ## Sparsity analysis
  sparse_dt <- copy(dt)
  sparse_dt[,method := method_recode[method]]
  sparse_summ <- sparse_dt[,.("sparsity" = format_col(sparsity)), by = .(method, band)]
  sparse_summ_tbl <- data.table::dcast(sparse_summ, band ~ method, value.var = "sparsity")
  sparse_summ_tbl <- order_cols_rows(sparse_summ_tbl, c("band", method_recode))

  
  table_header_spars <- c(1, ncol(sparse_summ_tbl)-1)
  names(table_header_spars) <- c(" ", "Sparsity")
  kbl_latex_spars <- kable(sparse_summ_tbl, format = "latex", align = paste0(rep("c",ncol(sparse_summ_tbl))), booktabs = T) %>% add_header_above(table_header_spars)
  kbl_spars <- kable(sparse_summ_tbl, format = "pipe", align = paste0(c("l", rep("c",ncol(sparse_summ_tbl)-1))), booktabs = T)
  
  ## Invertibility analysis
  invertible_summ <- sparse_dt[,.("% not solvable" = round(sum(is.na(sparsity))/(uniqueN(sub)*uniqueN(frequency)),2)), by = .(method, band)]
  invertible_summ_tbl <- data.table::dcast(invertible_summ, band ~ method, value.var = "% not solvable")
  invertible_summ_tbl <- order_cols_rows(invertible_summ_tbl, c("band", method_recode))
  
  table_header_inv <- c(1, ncol(invertible_summ_tbl)-1)
  names(table_header_inv) <- c(" ", "Sparsity")
  kbl_latex_inv <- kable(invertible_summ_tbl, format = "latex", align = paste0(rep("c",ncol(invertible_summ_tbl))), booktabs = T) %>% add_header_above(table_header_spars)
  
  return(list(table_sparsity = kbl_latex_spars, table_inv = kbl_latex_inv, table_sparsity_html = kbl_spars))
}

make_sparsity_figures <- function(dt, analysis_type){
  method_recode <- c("dtl_est" = "SDD", "naiveGLASSODiff" = "Naive", "naiveHardDiff" = "Hard threshold", "FGLDiff" = "FGL")
  col_scale <- c("SDD" = "#785ef0", "Naive" = "#648FFF", "Hard threshold" = "#FE6100", "FGL" = "#DC267F")
  band_recode <- c("theta" = "Theta", "beta" = "Beta", "gamma" = "Gamma", "high_gamma" = "High-gamma")
  
  
  ## Sparsity analysis
  sparse_dt <- copy(dt)
  sparse_dt[,method := method_recode[method]]
  sparse_summ <- sparse_dt[,.("sparsity" = mean(sparsity, na.rm = T), se = sqrt(var(sparsity, na.rm = T)/sum(!is.na(sparsity))), n_all = .N, n_notna = sum(!is.na(sparsity))), by = .(method, band)]
  sparse_summ$band <- factor(band_recode[sparse_summ$band], levels = c("Theta", "Beta", "Gamma", "High-gamma"))
  
  tikz(file = paste0(fig_folder, "/sparsity_", analysis_type, ".tex"), width = 4, height = 4, standAlone = TRUE)
  plt <- ggplot(sparse_summ, aes(x = band, y = sparsity, color = method)) +  geom_point(size = 2, position = position_dodge(width = 0.5)) + 
    geom_errorbar(aes(ymin=pmax(sparsity-se,0), ymax=pmin(sparsity+se,1)), width=0, linewidth= 1, position = position_dodge(width = 0.5)) +
    scale_color_manual(name = "Method", values = col_scale) + 
    xlab("Band") + ylab("Sparsity") + 
    theme_Publication(base_size = 6) +   theme(legend.background = element_rect(fill='transparent', color = NA), 
                                               legend.box.background = element_rect(fill='transparent', color = NA),
                                               panel.background = element_rect(fill='transparent', color = NA),
                                               plot.background = element_rect(fill='transparent', color = NA))
  print(plt)
  dev.off()
  return(plt)  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#    Finding subjects corresponding to each analysis   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

setwd(dat_folder)
skip_subjs <- c("sub-100", "sub-104")

subjs_w_s1 <- unique(gsub(".*(sub-\\d{3}).*", "\\1", list.dirs()[grepl("ses-t1",list.dirs())]))
subjs_w_s1 <- subjs_w_s1[!subjs_w_s1 %in% skip_subjs]
subjs_w_s2 <- unique(gsub(".*(sub-\\d{3}).*", "\\1", list.dirs()[grepl("ses-t2",list.dirs())]))
subjs_w_s2 <- subjs_w_s2[!subjs_w_s2 %in% skip_subjs]

subj_across_session <- intersect(subjs_w_s1, subjs_w_s2)
subj_within_session <- setdiff(subjs_w_s1, subjs_w_s2)

bands_dt <- data.table(band_nm = c("theta", "beta", "gamma", "high_gamma"), band_min = c(4,12,30,71), band_max = c(8,28,70,150))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#    Across sessions analysis   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

across_sessions <- lapply(subj_across_session, function(sub) analyzeOneSub(sub, bands_dt, dat_folder, res_folder))
across_sessions_sparsity <- rbindlist(lapply(across_sessions, "[[", "sparsity"))
across_sparsity_tables <- make_sparsity_tables(across_sessions_sparsity)
make_sparsity_figures(across_sessions_sparsity, "across_sessions")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#    Within sessions analysis   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

within_sessions <- lapply(subj_within_session, function(sub) analyzeOneSub(sub, bands_dt, dat_folder, res_folder))
within_sessions_sparsity <- rbindlist(lapply(within_sessions, "[[", "sparsity"))
within_sparsity_tables <- make_sparsity_tables(within_sessions_sparsity)
make_sparsity_figures(within_sessions_sparsity, "within_sessions")
