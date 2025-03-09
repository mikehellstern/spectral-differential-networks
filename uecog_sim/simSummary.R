## Library loads
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(data.table)
library(reshape2)
library(ggplot2)
library(ggforce)
library(extrafont)

freqs_coding <-  c(setNames(rep("Theta",5), c(4,5,6,7,8)), # Theta
                   setNames(rep("Beta",5), c(12,16,20,24,28)), # Beta
                   setNames(rep("Gamma",5), c(30,40,50,60,70)), # Gamma
                   setNames(rep("High-gamma",6), c(80,95,110,125,140,150))) # High-gamma #Freq in Hz
## Inputs
control_sessions <- paste0("session", 33:36)
experiment_sessions <- paste0("session", 1:32)
session_keys <- c(control_sessions, experiment_sessions)
block_diff_coding <- c(setNames(rep("RS to RS", 5), c("RecBlock1_RecBlock2", "RecBlock2_RecBlock3", "RecBlock3_RecBlock4", "RecBlock4_RecBlock5", "RecBlock5_RecBlock6")),
                       setNames(rep("RS to SS", 5), c("RecBlock1_CondBlock1", "RecBlock2_CondBlock2", "RecBlock3_CondBlock3", "RecBlock4_CondBlock4", "RecBlock5_CondBlock5")))


## For each session, compute edge dists for all block_diffs / frequencies
get_session_results <- function(session){
  diffnetests_all <- readRDS(paste0("./results/",session))
  params <- readRDS(paste0("./data/",session, "/params.rds"))
  
  laser_from <- params$stimfrom
  laser_to <- params$stimto
  
  ### Computing edge dists for all block_diffs / frequencies
  diffnetests_dists_l <- sapply(diffnetests_all, function(diffnetests){
    
    res_l <- sapply(diffnetests, function(x){
      dtl_est <- x$X_Y_diff$dtl_est
      dtl_est[upper.tri(dtl_est)] <- NA
      dtl_est_m <- melt(dtl_est, na.rm = T, value.name = "edgew", varnames = c("node1", "node2"))
      return(dtl_est_m)
    }, simplify = FALSE, USE.NAMES = TRUE)
    res_dt <- rbindlist(res_l, idcol = "freq")
    
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  diffnetests_dists <- rbindlist(diffnetests_dists_l, idcol = "block_diff")
  diffnetests_dists$laser_to <- ifelse(is.null(laser_to), NA, laser_to)
  diffnetests_dists$laser_from <- ifelse(is.null(laser_from), NA, laser_from)
  diffnetests_dists$delay <- params$delay
  
  return(diffnetests_dists)
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  Summarizing session data  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

sessions_dt <- rbindlist(sapply(session_keys, function(session) {print(session); get_session_results(session)}, simplify = FALSE, USE.NAMES = TRUE), idcol = "session")
sessions_dt[,c("block_diff_type", "band") := .(block_diff_coding[block_diff], freqs_coding[freq])]
sessions_dt$band <- factor(sessions_dt$band, levels = c("Theta", "Beta", "Gamma", "High-gamma"))
control_dt <- sessions_dt[delay %in% c("control")]
experimental_dt <- sessions_dt[delay %in% c("10ms", "100ms")]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#    Probability of edge between stim sites   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

delay_color_map <- c("100ms" = "#DC267F", "10ms" = "#648FFF", "control" = "#FE6100")

experimental_stim_edge <- experimental_dt[(node1 == laser_to & node2 == laser_from) | (node1 == laser_from & node2 == laser_to)]
n_edge_experimental <- nrow(experimental_stim_edge)
experimental_stim_edge <- experimental_stim_edge[,.(prob_edge = mean(abs(edgew)!=0), se_edge = sd(abs(edgew)!=0) / sqrt(.N), .N), by = .(block_diff_type, delay, band)]

control_stim_edge <- control_dt[sample(1:nrow(control_dt), n_edge_experimental)]
control_stim_edge <- control_stim_edge[,.(prob_edge = mean(abs(edgew)!=0), se_edge = sd(abs(edgew)!=0) / sqrt(.N), .N), by = .(block_diff_type, delay, band)]

stim_edges <- rbind(control_stim_edge, experimental_stim_edge)

#------- RS to RS
tikz(file = paste0("./results/prob_edge_stim_sites_rs_to_rs.tex"), width = 4, height = 4, standAlone = TRUE)
plt_rs_to_rs <- ggplot(stim_edges[block_diff_type == "RS to RS"], aes(x = band, y = prob_edge, color = delay)) + 
  geom_point(size = 2, position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin=prob_edge-se_edge, ymax=prob_edge+se_edge), width=0, linewidth=1, position = position_dodge(width = 0.5)) + 
  scale_color_manual(name = "Delay", values = delay_color_map) + 
  xlab("Band") + ylab("Prob. of edge change between stim sites") + 
  theme_Publication(base_size = 6) +   theme(legend.background = element_rect(fill='transparent', color = NA), 
                                             legend.box.background = element_rect(fill='transparent', color = NA),
                                             panel.background = element_rect(fill='transparent', color = NA),
                                             plot.background = element_rect(fill='transparent', color = NA))
plt_rs_to_rs
dev.off()

#------- RS to SS
tikz(file = paste0("./results/prob_edge_stim_sites_rs_to_ss.tex"), width = 4, height = 4, standAlone = TRUE)
plt_rs_to_ss <- ggplot(stim_edges[block_diff_type == "RS to SS"], aes(x = band, y = prob_edge, color = delay)) + 
  geom_point(size = 2, position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin=prob_edge-se_edge, ymax=prob_edge+se_edge), width=0, linewidth=1, position = position_dodge(width = 0.5)) + 
  scale_color_manual(name = "Delay", values = delay_color_map) + 
  xlab("Band") + ylab("Prob. of edge change between stim sites") + 
  theme_Publication(base_size = 6) +   theme(legend.background = element_rect(fill='transparent', color = NA), 
                                             legend.box.background = element_rect(fill='transparent', color = NA),
                                             panel.background = element_rect(fill='transparent', color = NA),
                                             plot.background = element_rect(fill='transparent', color = NA))
plt_rs_to_ss
dev.off()
