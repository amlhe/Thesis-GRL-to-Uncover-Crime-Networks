rm(list=ls())
source("Source/esbm.R")

library(sbm)
library(igraph)
library(aricode)
library(abind)
library(RColorBrewer)
library(greed)
library(igraph)
library(Matrix)
library(blockmodels)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(reshape)
library(gdata)
library(mcclust.ext)
library(gridExtra)
library(grid)
library(cowplot)
library(coda)
library(dummies)
library(randnet)
library(LaplacesDemon)
library(kernlab)
library(clue)
library(psych)
library(Matrix)
library(pROC)
library(cvms)
library(PRROC)
library(caret)
library("gridGraphics")

################################################################################
load("Adjacency matrices/adj_MN.RData")
load("Adjacency matrices/adj_PC.RData")
load("Adjacency matrices/adj_SN.RData")
load("Adjacency matrices/adj_WR.RData")
load("Adjacency matrices/adj_AW.RData")
load("Adjacency matrices/adj_JU.RData")

load("Adjacency matrices/adj_MN_sub.RData")
load("Adjacency matrices/adj_PC_sub.RData")


################################################################################
### Plot initial adjacency matrix
# Plot the reordered adjacency matrix using the image function
par(mai = c(0.1, 0.1, 0.1, 0.1)) 
image(1:nrow(MN_sub), 1:ncol(MN_sub), MN_sub, 
      col = colorRampPalette(brewer.pal(9, "Greys")[c(1, 8)])(30),
      axes = FALSE, xlab = "", ylab = "", asp = 1, xlim = c(0, nrow(MN_sub)+1), ylim = c(0,ncol(MN_sub)+1))
rect(0, 0, ncol(MN_sub)+1, nrow(MN_sub)+1, border = "black", lwd = 2)
fun_grob <- function(){
  grid.echo()
  grid.grab()
}
fun_grob_out <- fun_grob()
pushViewport(viewport(width = 1, height = 1, angle = 270))
grid.draw(fun_grob_out)


par(mai = c(0.1, 0.1, 0.1, 0.1)) 
image(1:nrow(PC_sub), 1:ncol(PC_sub), PC_sub, 
      col = colorRampPalette(brewer.pal(9, "Greys")[c(1, 8)])(30),
      axes = FALSE, xlab = "", ylab = "", asp = 1, xlim = c(0, nrow(PC_sub)+1), ylim = c(0,ncol(PC_sub)+1))
rect(0, 0, ncol(PC_sub)+1, nrow(PC_sub)+1, border = "black", lwd = 2)
fun_grob <- function(){
  grid.echo()
  grid.grab()
}
fun_grob_out <- fun_grob()
pushViewport(viewport(width = 1, height = 1, angle = 270))
grid.draw(fun_grob_out)


################################################################################
num_nodes <- 47
sparse_i_MN_sub <- read.table("Missing links/MN_sub/sparse_i.txt", header = FALSE)$V1
sparse_j_MN_sub <- read.table("Missing links/MN_sub/sparse_j.txt", header = FALSE)$V1
MN_sub_mis <- matrix(0, nrow = num_nodes, ncol = num_nodes)
for (i in 1:length(sparse_i_MN_sub)) {
  MN_sub_mis[sparse_i_MN_sub[i], sparse_j_MN_sub[i]] <- 1
  MN_sub_mis[sparse_j_MN_sub[i], sparse_i_MN_sub[i]] <- 1
}

sparse_i_PC_sub <- read.table("Missing links/PC_sub/sparse_i.txt", header = FALSE)$V1
sparse_j_PC_sub <- read.table("Missing links/PC_sub/sparse_j.txt", header = FALSE)$V1
PC_sub_mis <- matrix(0, nrow = num_nodes, ncol = num_nodes)
for (i in 1:length(sparse_i_PC_sub)) {
  PC_sub_mis[sparse_i_PC_sub[i], sparse_j_PC_sub[i]] <- 1
  PC_sub_mis[sparse_j_PC_sub[i], sparse_i_PC_sub[i]] <- 1
}

N_iter <- 10000
my_seed <- 1


################################################################################
num_nodes <- 182
sparse_i_WR <- read.table("Missing links/WR/sparse_i.txt", header = FALSE)$V1
sparse_j_WR <- read.table("Missing links/WR/sparse_j.txt", header = FALSE)$V1
WR_mis <- matrix(0, nrow = num_nodes, ncol = num_nodes)
for (i in 1:length(sparse_i_WR)) {
  WR_mis[sparse_i_WR[i], sparse_j_WR[i]] <- 1
  WR_mis[sparse_j_WR[i], sparse_i_WR[i]] <- 1
}

sparse_i_AW <- read.table("Missing links/AW/sparse_i.txt", header = FALSE)$V1
sparse_j_AW <- read.table("Missing links/AW/sparse_j.txt", header = FALSE)$V1
AW_mis <- matrix(0, nrow = num_nodes, ncol = num_nodes)
for (i in 1:length(sparse_i_AW)) {
  AW_mis[sparse_i_AW[i], sparse_j_AW[i]] <- 1
  AW_mis[sparse_j_AW[i], sparse_i_AW[i]] <- 1
}

sparse_i_JU <- read.table("Missing links/JU/sparse_i.txt", header = FALSE)$V1
sparse_j_JU <- read.table("Missing links/JU/sparse_j.txt", header = FALSE)$V1
JU_mis <- matrix(0, nrow = num_nodes, ncol = num_nodes)
for (i in 1:length(sparse_i_JU)) {
  JU_mis[sparse_i_JU[i], sparse_j_JU[i]] <- 1
  JU_mis[sparse_j_JU[i], sparse_i_JU[i]] <- 1
}


N_iter <- 10000
my_seed <- 1

################################################################################
# Create multilayer lists
Y_list_Montagna <- abind(MN_sub, PC_sub, along = 3)
Y_list_Montagna_mis <- abind(MN_sub_mis, PC_sub_mis, along = 3)
Y_list_Oversize <- abind(WR, AW, JU, along = 3)
Y_list_Oversize_mis <- abind(WR_mis, AW_mis, JU_mis, along = 3)

################################################################################
# MONTAGNA

############################################################
# Posterior computation using collapsed Gibbs sampler
############################################################

# DIRICHLET MULTINOMIAL
# ------------------------------------
my_prior <- "DM"
M_shared_Z_DM <- multilayer_esbm_shared(Y_list_Montagna, my_seed, N_iter, my_prior, a = 1, b = 1, beta_DM = 12/50, H_DM = 50)
M_indi_Z_DM <- multilayer_esbm_layer_specific(Y_list_Montagna, my_seed, N_iter, my_prior, a = 1, b = 1, beta_DM = 12/50, H_DM = 50)
M_shared_Z_DM_mis <- multilayer_esbm_shared(Y_list_Montagna_mis, my_seed, N_iter, my_prior, a = 1, b = 1, beta_DM = 12/50, H_DM = 50)
M_indi_Z_DM_mis <- multilayer_esbm_layer_specific(Y_list_Montagna_mis, my_seed, N_iter, my_prior, a = 1, b = 1, beta_DM = 12/50, H_DM = 50)

# DIRICHLET PROCESS (CRP)
# ------------------------------------
my_prior <- "DP"
M_shared_Z_DP <- multilayer_esbm_shared(Y_list_Montagna, my_seed, N_iter, my_prior, a = 1, b = 1, alpha_PY = 8, sigma_PY = 0)
M_indi_Z_DP <- multilayer_esbm_layer_specific(Y_list_Montagna, my_seed, N_iter, my_prior, a = 1, b = 1, alpha_PY = 8, sigma_PY = 0)
M_shared_Z_DP_mis <- multilayer_esbm_shared(Y_list_Montagna_mis, my_seed, N_iter, my_prior, a = 1, b = 1, alpha_PY = 8, sigma_PY = 0)
M_indi_Z_DP_mis <- multilayer_esbm_layer_specific(Y_list_Montagna_mis, my_seed, N_iter, my_prior, a = 1, b = 1, alpha_PY = 8, sigma_PY = 0)

# PITMAN-YOR PROCESS
# ------------------------------------
my_prior <- "PY"
M_shared_Z_PY <- multilayer_esbm_shared(Y_list_Montagna, my_seed, N_iter, my_prior, a = 1, b = 1, alpha_PY = -0.350, sigma_PY = 0.725)
M_indi_Z_PY <- multilayer_esbm_layer_specific(Y_list_Montagna, my_seed, N_iter, my_prior, a = 1, b = 1, alpha_PY = -0.350, sigma_PY = 0.725)
M_shared_Z_PY_mis <- multilayer_esbm_shared(Y_list_Montagna_mis, my_seed, N_iter, my_prior, a = 1, b = 1, alpha_PY = -0.350, sigma_PY = 0.725)
M_indi_Z_PY_mis <- multilayer_esbm_layer_specific(Y_list_Montagna_mis, my_seed, N_iter, my_prior, a = 1, b = 1, alpha_PY = -0.350, sigma_PY = 0.725)

# GNEDIN PROCESS
# ------------------------------------
my_prior  <- "GN"
M_shared_Z_GN <- multilayer_esbm_shared(Y_list_Montagna, my_seed, N_iter, my_prior, a = 1, b = 1, gamma_GN = 0.3)
M_indi_Z_GN <- multilayer_esbm_layer_specific(Y_list_Montagna, my_seed, N_iter, my_prior, a = 1, b = 1, gamma_GN = 0.3)
M_shared_Z_GN_mis <- multilayer_esbm_shared(Y_list_Montagna_mis, my_seed, N_iter, my_prior, a = 1, b = 1, gamma_GN = 0.3)
M_indi_Z_GN_mis <- multilayer_esbm_layer_specific(Y_list_Montagna_mis, my_seed, N_iter, my_prior, a = 1, b = 1, gamma_GN = 0.3)


# Save the output
save(M_shared_Z_DM,M_shared_Z_DP,M_shared_Z_PY,M_shared_Z_GN,file="Application/Posterior_Montagna_shared_fixed.RData")
save(M_indi_Z_DM,M_indi_Z_DP,M_indi_Z_PY,M_indi_Z_GN,file="Application/Posterior_Montagna_indi_fixed.RData")
save(M_shared_Z_DM_mis,M_shared_Z_DP_mis,M_shared_Z_PY_mis,M_shared_Z_GN_mis,file="Application/Posterior_Montagna_shared_mis.RData")
save(M_indi_Z_DM_mis,M_indi_Z_DP_mis,M_indi_Z_PY_mis,M_indi_Z_GN_mis,file="Application/Posterior_Montagna_indi_mis.RData")
load("Application/Posterior_Montagna_shared_fixed.RData")
load("Application/Posterior_Montagna_indi_fixed.RData")
load("Application/Posterior_Montagna_shared_mis.RData")
load("Application/Posterior_Montagna_indi_mis.RData")

############################################################
# Posterior inference
############################################################

set.seed(1)
V <- dim(MN_sub)[1]
index_traceplot <- sample(c(1:(V*(V-1)/2)),1)

burn_in <- 10000
a <- b <- 1

### DIRICHLET PROCESS
M_shared_Z_DP_WAIC <- M_shared_Z_DM[,(burn_in+1):N_iter]

for (t in 1:dim(M_shared_Z_DP_WAIC)[2]){
  LL[,t]<-sampleLL(M_shared_Z_DP_WAIC[,t],MN_sub,a,b)
  if (t%%1000 == 0){print(paste("Iteration:", t))}
}
WAIC(LL)$WAIC

# Shared eta
c_DP <- pr_cc(M_shared_Z_DP[,(burn_in+1):N_iter])
c_DM <- pr_cc(M_shared_Z_DM[,(burn_in+1):N_iter])
c_PY <- pr_cc(M_shared_Z_PY[,(burn_in+1):N_iter])
c_GN <- pr_cc(M_shared_Z_GN[,(burn_in+1):N_iter])

memb_DP_VI <- minVI(c_DP,method="avg",max.k=30)
memb_DM_VI <- minVI(c_DM,method="avg",max.k=30)
memb_PY_VI <- minVI(c_PY,method="avg",max.k=30)
memb_GN_VI <- minVI(c_GN,method="avg",max.k=30)

memb_DP_shared <- memb_DP_VI$cl
memb_DM_shared <- memb_DM_VI$cl
memb_PY_shared <- memb_PY_VI$cl
memb_GN_shared <- memb_GN_VI$cl

# Plot reordered adjacency matrix shared
par(mai = c(0.1, 0.1, 0.1, 0.1))
plot_reordered_adjacency(memb_DP_shared, MN_sub)
plot_reordered_adjacency(memb_DM_shared, MN_sub)
plot_reordered_adjacency(memb_PY_shared, MN_sub)
plot_reordered_adjacency(memb_GN_shared, MN_sub)

plot_reordered_adjacency(memb_DP_shared, PC_sub)
plot_reordered_adjacency(memb_DM_shared, PC_sub)
plot_reordered_adjacency(memb_PY_shared, PC_sub)
plot_reordered_adjacency(memb_GN_shared, PC_sub)



# Layer-specific eta
c_DP <- pr_cc(M_indi_Z_DP[,(burn_in+1):N_iter])
c_DM <- pr_cc(M_indi_Z_DM[,(burn_in+1):N_iter])
c_PY <- pr_cc(M_indi_Z_PY[,(burn_in+1):N_iter])
c_GN <- pr_cc(M_indi_Z_GN[,(burn_in+1):N_iter])

memb_DP_VI <- minVI(c_DP,method="avg",max.k=30)
memb_DM_VI <- minVI(c_DM,method="avg",max.k=30)
memb_PY_VI <- minVI(c_PY,method="avg",max.k=30)
memb_GN_VI <- minVI(c_GN,method="avg",max.k=30)

memb_DP_indi <- memb_DP_VI$cl
memb_DM_indi <- memb_DM_VI$cl
memb_PY_indi <- memb_PY_VI$cl
memb_GN_indi <- memb_GN_VI$cl

# Plot reordered adjacency matrix layer-specific
par(mai = c(0.1, 0.1, 0.1, 0.1))
plot_reordered_adjacency(memb_DP_indi, MN_sub)
plot_reordered_adjacency(memb_DM_indi, MN_sub)
plot_reordered_adjacency(memb_PY_indi, MN_sub)
plot_reordered_adjacency(memb_GN_indi, MN_sub)

plot_reordered_adjacency(memb_DP_indi, PC_sub)
plot_reordered_adjacency(memb_DM_indi, PC_sub)
plot_reordered_adjacency(memb_PY_indi, PC_sub)
plot_reordered_adjacency(memb_GN_indi, PC_sub)

# Plot graph with nodes colored according to community membership
block_colors <- rainbow(length(unique(memb_Z_DP_sel)))
node_colors <- block_colors[memb_Z_DP_sel]
plot(graph_from_adjacency_matrix(PC_sub, mode = "undirected", diag = FALSE), vertex.size = 10, vertex.color = node_colors, vertex.label = NA, main = paste("Number of communities:", length(unique(memb_Z_DP))))

############################################################
# Posterior inference - missing links
############################################################
set.seed(1)
index_traceplot <- sample(c(1:(V*(V-1)/2)),1)

V <- dim(MN_sub_mis)[1]
burn_in <- 1000
a <- b <- 1

# Shared eta
c_DP_mis <- pr_cc(M_shared_Z_DP_mis[,(burn_in+1):N_iter])
c_DM_mis <- pr_cc(M_shared_Z_DM_mis[,(burn_in+1):N_iter])
c_PY_mis <- pr_cc(M_shared_Z_PY_mis[,(burn_in+1):N_iter])
c_GN_mis <- pr_cc(M_shared_Z_GN_mis[,(burn_in+1):N_iter])

memb_DP_VI_mis <- minVI(c_DP_mis,method="avg",max.k=30)
memb_DM_VI_mis <- minVI(c_DM_mis,method="avg",max.k=30)
memb_PY_VI_mis <- minVI(c_PY_mis,method="avg",max.k=30)
memb_GN_VI_mis <- minVI(c_GN_mis,method="avg",max.k=30)

memb_DP_shared_mis <- memb_DP_VI_mis$cl
memb_DM_shared_mis <- memb_DM_VI_mis$cl
memb_PY_shared_mis <- memb_PY_VI_mis$cl
memb_GN_shared_mis <- memb_GN_VI_mis$cl

# ------------------------------------
# ROC collection form different priors
# Collect ROC curves and AUC-ROC scores for all priors
roc_curves_mis <- list()
auc_scores_mis <- numeric()

# DIRICHLET PROCESS (CRP) UNSUPERVISED
result_DP_mis <- misclass_roc(memb_DP_shared_mis, PC_sub, a = 1, b = 1)
roc_curves_mis[["DP"]] <- result_DP_mis$perf
auc_scores_mis[["DP"]] <- result_DP_mis$auc_roc

# PITMAN-YOR (PY)
result_PY_mis <- misclass_roc(memb_PY_shared_mis, PC_sub, a = 1, b = 1)
roc_curves_mis[["PY"]] <- result_PY_mis$perf
auc_scores_mis[["PY"]] <- result_PY_mis$auc_roc

# DIRICHLET MULTINOMIAL (DM)
result_DM_mis <- misclass_roc(memb_DM_shared_mis, PC_sub, a = 1, b = 1)
roc_curves_mis[["DM"]] <- result_DM_mis$perf
auc_scores_mis[["DM"]] <- result_DM_mis$auc_roc

# GNEDIN PROCESS (GN)
result_GN_mis <- misclass_roc(memb_GN_shared_mis, PC_sub, a = 1, b = 1)
roc_curves_mis[["GN"]] <- result_GN_mis$perf
auc_scores_mis[["GN"]] <- result_GN_mis$auc_roc

# Plot ROC curves with different colors
colors <- c("red", "green", "blue", "purple")  # Add more colors if needed
legend_text <- c("DP", "PY", "DM", "GN")  # Legend labels

# Create an empty plot
plot(roc_curves_mis[["DP"]], col = colors[1], lwd = 1)

# Add diagonal line for reference (color set to grey)
abline(a = 0, b = 1, col = "grey", lwd = 1, lty = 3)

# Iterate through the list and add each ROC curve to the plot
for (i in seq_along(roc_curves_mis)) {
  lines(roc_curves_mis[[i]]@x.values[[1]], roc_curves_mis[[i]]@y.values[[1]], col = colors[i], lwd = 1)
}

# Add a legend with AUC-ROC scores
legend("bottomright", legend = sprintf("%s AUC-ROC = %.2f", legend_text, auc_scores_mis),
       col = colors, lty = 1, cex = 0.8)


write.csv(memb_DP_shared_mis, "../ESBM ROC/memb_Z_DP_mis.csv", row.names = FALSE)
write.csv(memb_DM_shared_mis, "../ESBM ROC/memb_Z_DM_mis.csv", row.names = FALSE)
write.csv(memb_PY_shared_mis, "../ESBM ROC/memb_Z_PY_mis.csv", row.names = FALSE)
write.csv(memb_GN_shared_mis, "../ESBM ROC/memb_Z_GN_mis.csv", row.names = FALSE)
write.csv(PC_sub, "../ESBM ROC/Y.csv", row.names = FALSE)



# Layer-specific eta
c_DP_mis <- pr_cc(M_indi_Z_DP_mis[,(burn_in+1):N_iter])
c_DM_mis <- pr_cc(M_indi_Z_DM_mis[,(burn_in+1):N_iter])
c_PY_mis <- pr_cc(M_indi_Z_PY_mis[,(burn_in+1):N_iter])
c_GN_mis <- pr_cc(M_indi_Z_GN_mis[,(burn_in+1):N_iter])

memb_DP_VI_mis <- minVI(c_DP_mis,method="avg",max.k=30)
memb_DM_VI_mis <- minVI(c_DM_mis,method="avg",max.k=30)
memb_PY_VI_mis <- minVI(c_PY_mis,method="avg",max.k=30)
memb_GN_VI_mis <- minVI(c_GN_mis,method="avg",max.k=30)

memb_DP_indi_mis <- memb_DP_VI_mis$cl
memb_DM_indi_mis <- memb_DM_VI_mis$cl
memb_PY_indi_mis <- memb_PY_VI_mis$cl
memb_GN_indi_mis <- memb_GN_VI_mis$cl

# ------------------------------------
# ROC collection form different priors
# Collect ROC curves and AUC-ROC scores for all priors
roc_curves_mis <- list()
auc_scores_mis <- numeric()

# DIRICHLET PROCESS (CRP) UNSUPERVISED
result_DP_mis <- misclass_roc(memb_DP_indi_mis, PC_sub, a = 1, b = 1)
roc_curves_mis[["DP"]] <- result_DP_mis$perf
auc_scores_mis[["DP"]] <- result_DP_mis$auc_roc

# PITMAN-YOR (PY)
result_PY_mis <- misclass_roc(memb_PY_indi_mis, PC_sub, a = 1, b = 1)
roc_curves_mis[["PY"]] <- result_PY_mis$perf
auc_scores_mis[["PY"]] <- result_PY_mis$auc_roc

# DIRICHLET MULTINOMIAL (DM)
result_DM_mis <- misclass_roc(memb_DM_indi_mis, PC_sub, a = 1, b = 1)
roc_curves_mis[["DM"]] <- result_DM_mis$perf
auc_scores_mis[["DM"]] <- result_DM_mis$auc_roc

# GNEDIN PROCESS (GN)
result_GN_mis <- misclass_roc(memb_GN_indi_mis, PC_sub, a = 1, b = 1)
roc_curves_mis[["GN"]] <- result_GN_mis$perf
auc_scores_mis[["GN"]] <- result_GN_mis$auc_roc

# Plot ROC curves with different colors
colors <- c("red", "green", "blue", "purple")  # Add more colors if needed
legend_text <- c("DP", "PY", "DM", "GN")  # Legend labels

# Create an empty plot
plot(roc_curves_mis[["DP"]], col = colors[1], lwd = 1)

# Add diagonal line for reference (color set to grey)
abline(a = 0, b = 1, col = "grey", lwd = 1, lty = 3)

# Iterate through the list and add each ROC curve to the plot
for (i in seq_along(roc_curves_mis)) {
  lines(roc_curves_mis[[i]]@x.values[[1]], roc_curves_mis[[i]]@y.values[[1]], col = colors[i], lwd = 1)
}

# Add a legend with AUC-ROC scores
legend("bottomright", legend = sprintf("%s AUC-ROC = %.2f", legend_text, auc_scores_mis),
       col = colors, lty = 1, cex = 0.8)

write.csv(memb_DP_indi_mis, "../ESBM ROC/memb_Z_DP_mis.csv", row.names = FALSE)
write.csv(memb_DM_indi_mis, "../ESBM ROC/memb_Z_DM_mis.csv", row.names = FALSE)
write.csv(memb_PY_indi_mis, "../ESBM ROC/memb_Z_PY_mis.csv", row.names = FALSE)
write.csv(memb_GN_indi_mis, "../ESBM ROC/memb_Z_GN_mis.csv", row.names = FALSE)
write.csv(PC_sub, "../ESBM ROC/Y.csv", row.names = FALSE)



################################################################################
# Oversize

############################################################
# Posterior computation using collapsed Gibbs sampler
############################################################

# DIRICHLET MULTINOMIAL
# ------------------------------------
my_prior <- "DM"
O_shared_Z_DM <- multilayer_esbm_shared(Y_list_Oversize, my_seed, N_iter, my_prior, a = 1, b = 1, beta_DM = 12/50, H_DM = 50)
O_indi_Z_DM <- multilayer_esbm_layer_specific(Y_list_Oversize, my_seed, N_iter, my_prior, a = 1, b = 1, beta_DM = 12/50, H_DM = 50)
O_shared_Z_DM_mis <- multilayer_esbm_shared(Y_list_Oversize_mis, my_seed, N_iter, my_prior, a = 1, b = 1, beta_DM = 12/50, H_DM = 50)
O_indi_Z_DM_mis <- multilayer_esbm_layer_specific(Y_list_Oversize_mis, my_seed, N_iter, my_prior, a = 1, b = 1, beta_DM = 12/50, H_DM = 50)

# DIRICHLET PROCESS (CRP)
# ------------------------------------
my_prior <- "DP"
O_shared_Z_DP <- multilayer_esbm_shared(Y_list_Oversize, my_seed, N_iter, my_prior, a = 1, b = 1, alpha_PY = 8, sigma_PY = 0)
O_indi_Z_DP <- multilayer_esbm_layer_specific(Y_list_Oversize, my_seed, N_iter, my_prior, a = 1, b = 1, alpha_PY = 8, sigma_PY = 0)
O_shared_Z_DP_mis <- multilayer_esbm_shared(Y_list_Oversize_mis, my_seed, N_iter, my_prior, a = 1, b = 1, alpha_PY = 8, sigma_PY = 0)
O_indi_Z_DP_mis <- multilayer_esbm_layer_specific(Y_list_Oversize_mis, my_seed, N_iter, my_prior, a = 1, b = 1, alpha_PY = 8, sigma_PY = 0)

# PITMAN-YOR PROCESS
# ------------------------------------
my_prior <- "PY"
O_shared_Z_PY <- multilayer_esbm_shared(Y_list_Oversize, my_seed, N_iter, my_prior, a = 1, b = 1, alpha_PY = -0.350, sigma_PY = 0.725)
O_indi_Z_PY <- multilayer_esbm_layer_specific(Y_list_Oversize, my_seed, N_iter, my_prior, a = 1, b = 1, alpha_PY = -0.350, sigma_PY = 0.725)
O_shared_Z_PY_mis <- multilayer_esbm_shared(Y_list_Oversize_mis, my_seed, N_iter, my_prior, a = 1, b = 1, alpha_PY = -0.350, sigma_PY = 0.725)
O_indi_Z_PY_mis <- multilayer_esbm_layer_specific(Y_list_Oversize_mis, my_seed, N_iter, my_prior, a = 1, b = 1, alpha_PY = -0.350, sigma_PY = 0.725)

# GNEDIN PROCESS
# ------------------------------------
my_prior  <- "GN"
O_shared_Z_GN <- multilayer_esbm_shared(Y_list_Oversize, my_seed, N_iter, my_prior, a = 1, b = 1, gamma_GN = 0.3)
O_indi_Z_GN <- multilayer_esbm_layer_specific(Y_list_Oversize, my_seed, N_iter, my_prior, a = 1, b = 1, gamma_GN = 0.3)
O_shared_Z_GN_mis <- multilayer_esbm_shared(Y_list_Oversize_mis, my_seed, N_iter, my_prior, a = 1, b = 1, gamma_GN = 0.3)
O_indi_Z_GN_mis <- multilayer_esbm_layer_specific(Y_list_Oversize_mis, my_seed, N_iter, my_prior, a = 1, b = 1, gamma_GN = 0.3)


# Save the output
save(O_shared_Z_DM,O_shared_Z_DP,O_shared_Z_PY,O_shared_Z_GN,file="Application/Posterior_Oversize_shared_fixed.RData")
save(O_indi_Z_DM,O_indi_Z_DP,O_indi_Z_PY,O_indi_Z_GN,file="Application/Posterior_Oversize_indi_fixed.RData")
save(O_shared_Z_DM_mis,O_shared_Z_DP_mis,O_shared_Z_PY_mis,O_shared_Z_GN_mis,file="Application/Posterior_Oversize_shared_mis.RData")
save(O_indi_Z_DM_mis,O_indi_Z_DP_mis,O_indi_Z_PY_mis,O_indi_Z_GN_mis,file="Application/Posterior_Oversize_indi_mis.RData")
load("Application/Posterior_Oversize_shared_fixed.RData")
load("Application/Posterior_Oversize_indi_fixed.RData")
load("Application/Posterior_Oversize_shared_mis.RData")
load("Application/Posterior_Oversize_indi_mis.RData")


############################################################
# Posterior inference
############################################################

set.seed(1)
V <- dim(WR)[1]
index_traceplot <- sample(c(1:(V*(V-1)/2)),1)

burn_in <- 10000
a <- b <- 1
LL <- matrix(nrow=V*(V-1)/2,ncol=40000)

### DIRICHLET PROCESS
M_shared_Z_DP_WAIC <- M_shared_Z_DM[,(burn_in+1):N_iter]

for (t in 1:dim(M_shared_Z_DP_WAIC)[2]){
  LL[,t]<-sampleLL(M_shared_Z_DP_WAIC[,t],MN_sub,a,b)
  if (t%%1000 == 0){print(paste("Iteration:", t))}
}
WAIC(LL)$WAIC

# Shared eta
c_DP <- pr_cc(O_shared_Z_DP[,(burn_in+1):N_iter])
c_DM <- pr_cc(O_shared_Z_DM[,(burn_in+1):N_iter])
c_PY <- pr_cc(O_shared_Z_PY[,(burn_in+1):N_iter])
c_GN <- pr_cc(O_shared_Z_GN[,(burn_in+1):N_iter])

memb_DP_VI <- minVI(c_DP,method="avg",max.k=30)
memb_DM_VI <- minVI(c_DM,method="avg",max.k=30)
memb_PY_VI <- minVI(c_PY,method="avg",max.k=30)
memb_GN_VI <- minVI(c_GN,method="avg",max.k=30)

memb_DP_shared <- memb_DP_VI$cl
memb_DM_shared <- memb_DM_VI$cl
memb_PY_shared <- memb_PY_VI$cl
memb_GN_shared <- memb_GN_VI$cl

# Plot reordered adjacency matrix shared
par(mai = c(0.1, 0.1, 0.1, 0.1))
plot_reordered_adjacency(memb_DP_shared, WR)
plot_reordered_adjacency(memb_DM_shared, WR)
plot_reordered_adjacency(memb_PY_shared, WR)
plot_reordered_adjacency(memb_GN_shared, WR)

plot_reordered_adjacency(memb_DP_shared, AW)
plot_reordered_adjacency(memb_DM_shared, AW)
plot_reordered_adjacency(memb_PY_shared, AW)
plot_reordered_adjacency(memb_GN_shared, AW)

plot_reordered_adjacency(memb_DP_shared, JU)
plot_reordered_adjacency(memb_DM_shared, JU)
plot_reordered_adjacency(memb_PY_shared, JU)
plot_reordered_adjacency(memb_GN_shared, JU)


# Layer-specific eta
c_DP <- pr_cc(O_indi_Z_DP[,(burn_in+1):N_iter])
c_DM <- pr_cc(O_indi_Z_DM[,(burn_in+1):N_iter])
c_PY <- pr_cc(O_indi_Z_PY[,(burn_in+1):N_iter])
c_GN <- pr_cc(O_indi_Z_GN[,(burn_in+1):N_iter])

memb_DP_VI <- minVI(c_DP,method="avg",max.k=30)
memb_DM_VI <- minVI(c_DM,method="avg",max.k=30)
memb_PY_VI <- minVI(c_PY,method="avg",max.k=30)
memb_GN_VI <- minVI(c_GN,method="avg",max.k=30)

memb_DP_indi <- memb_DP_VI$cl
memb_DM_indi <- memb_DM_VI$cl
memb_PY_indi <- memb_PY_VI$cl
memb_GN_indi <- memb_GN_VI$cl

# Plot reordered adjacency matrix layer-specific
par(mai = c(0.1, 0.1, 0.1, 0.1))
plot_reordered_adjacency(memb_DP_indi, WR)
plot_reordered_adjacency(memb_DM_indi, WR)
plot_reordered_adjacency(memb_PY_indi, WR)
plot_reordered_adjacency(memb_GN_indi, WR)

plot_reordered_adjacency(memb_DP_indi, AW)
plot_reordered_adjacency(memb_DM_indi, AW)
plot_reordered_adjacency(memb_PY_indi, AW)
plot_reordered_adjacency(memb_GN_indi, AW)

plot_reordered_adjacency(memb_DP_indi, JU)
plot_reordered_adjacency(memb_DM_indi, JU)
plot_reordered_adjacency(memb_PY_indi, JU)
plot_reordered_adjacency(memb_GN_indi, JU)

# Plot graph with nodes colored according to community membership
block_colors <- rainbow(length(unique(memb_Z_DP_sel)))
node_colors <- block_colors[memb_Z_DP_sel]
plot(graph_from_adjacency_matrix(PC_sub, mode = "undirected", diag = FALSE), vertex.size = 10, vertex.color = node_colors, vertex.label = NA, main = paste("Number of communities:", length(unique(memb_Z_DP))))



############################################################
# Posterior inference - missing links
############################################################
set.seed(1)
index_traceplot <- sample(c(1:(V*(V-1)/2)),1)

V <- dim(WR_mis)[1]
burn_in <- 1000
a <- b <- 1

# Shared eta
c_DP_mis <- pr_cc(O_shared_Z_DP_mis[,(burn_in+1):N_iter])
c_DM_mis <- pr_cc(O_shared_Z_DM_mis[,(burn_in+1):N_iter])
c_PY_mis <- pr_cc(O_shared_Z_PY_mis[,(burn_in+1):N_iter])
c_GN_mis <- pr_cc(O_shared_Z_GN_mis[,(burn_in+1):N_iter])

memb_DP_VI_mis <- minVI(c_DP_mis,method="avg",max.k=30)
memb_DM_VI_mis <- minVI(c_DM_mis,method="avg",max.k=30)
memb_PY_VI_mis <- minVI(c_PY_mis,method="avg",max.k=30)
memb_GN_VI_mis <- minVI(c_GN_mis,method="avg",max.k=30)

memb_DP_shared_mis <- memb_DP_VI_mis$cl
memb_DM_shared_mis <- memb_DM_VI_mis$cl
memb_PY_shared_mis <- memb_PY_VI_mis$cl
memb_GN_shared_mis <- memb_GN_VI_mis$cl

# ------------------------------------
# ROC collection form different priors
# Collect ROC curves and AUC-ROC scores for all priors
roc_curves_mis <- list()
auc_scores_mis <- numeric()

# DIRICHLET PROCESS (CRP) UNSUPERVISED
result_DP_mis <- misclass_roc(memb_DP_shared_mis, JU, a = 1, b = 1)
roc_curves_mis[["DP"]] <- result_DP_mis$perf
auc_scores_mis[["DP"]] <- result_DP_mis$auc_roc

# PITMAN-YOR (PY)
result_PY_mis <- misclass_roc(memb_PY_shared_mis, JU, a = 1, b = 1)
roc_curves_mis[["PY"]] <- result_PY_mis$perf
auc_scores_mis[["PY"]] <- result_PY_mis$auc_roc

# DIRICHLET MULTINOMIAL (DM)
result_DM_mis <- misclass_roc(memb_DM_shared_mis, JU, a = 1, b = 1)
roc_curves_mis[["DM"]] <- result_DM_mis$perf
auc_scores_mis[["DM"]] <- result_DM_mis$auc_roc

# GNEDIN PROCESS (GN)
result_GN_mis <- misclass_roc(memb_GN_shared_mis, JU, a = 1, b = 1)
roc_curves_mis[["GN"]] <- result_GN_mis$perf
auc_scores_mis[["GN"]] <- result_GN_mis$auc_roc

# Plot ROC curves with different colors
colors <- c("red", "green", "blue", "purple")  # Add more colors if needed
legend_text <- c("DP", "PY", "DM", "GN")  # Legend labels

# Create an empty plot
plot(roc_curves_mis[["DP"]], col = colors[1], lwd = 1)

# Add diagonal line for reference (color set to grey)
abline(a = 0, b = 1, col = "grey", lwd = 1, lty = 3)

# Iterate through the list and add each ROC curve to the plot
for (i in seq_along(roc_curves_mis)) {
  lines(roc_curves_mis[[i]]@x.values[[1]], roc_curves_mis[[i]]@y.values[[1]], col = colors[i], lwd = 1)
}

# Add a legend with AUC-ROC scores
legend("bottomright", legend = sprintf("%s AUC-ROC = %.2f", legend_text, auc_scores_mis),
       col = colors, lty = 1, cex = 0.8)



write.csv(memb_DP_shared_mis, "../ESBM ROC/memb_Z_DP_mis.csv", row.names = FALSE)
write.csv(memb_DM_shared_mis, "../ESBM ROC/memb_Z_DM_mis.csv", row.names = FALSE)
write.csv(memb_PY_shared_mis, "../ESBM ROC/memb_Z_PY_mis.csv", row.names = FALSE)
write.csv(memb_GN_shared_mis, "../ESBM ROC/memb_Z_GN_mis.csv", row.names = FALSE)
write.csv(JU, "../ESBM ROC/Y.csv", row.names = FALSE)



# Layer-specific eta
c_DP_mis <- pr_cc(O_indi_Z_DP_mis[,(burn_in+1):N_iter])
c_DM_mis <- pr_cc(O_indi_Z_DM_mis[,(burn_in+1):N_iter])
c_PY_mis <- pr_cc(O_indi_Z_PY_mis[,(burn_in+1):N_iter])
c_GN_mis <- pr_cc(O_indi_Z_GN_mis[,(burn_in+1):N_iter])

memb_DP_VI_mis <- minVI(c_DP_mis,method="avg",max.k=30)
memb_DM_VI_mis <- minVI(c_DM_mis,method="avg",max.k=30)
memb_PY_VI_mis <- minVI(c_PY_mis,method="avg",max.k=30)
memb_GN_VI_mis <- minVI(c_GN_mis,method="avg",max.k=30)

memb_DP_indi_mis <- memb_DP_VI_mis$cl
memb_DM_indi_mis <- memb_DM_VI_mis$cl
memb_PY_indi_mis <- memb_PY_VI_mis$cl
memb_GN_indi_mis <- memb_GN_VI_mis$cl

# ------------------------------------
# ROC collection form different priors
# Collect ROC curves and AUC-ROC scores for all priors
roc_curves_mis <- list()
roc_curves_mis <- list()
auc_scores_mis <- numeric()

# DIRICHLET PROCESS (CRP) UNSUPERVISED
result_DP_mis <- misclass_roc(memb_DP_indi_mis, JU, a = 1, b = 1)
roc_curves_mis[["DP"]] <- result_DP_mis$perf
auc_scores_mis[["DP"]] <- result_DP_mis$auc_roc

# PITMAN-YOR (PY)
result_PY_mis <- misclass_roc(memb_PY_indi_mis, JU, a = 1, b = 1)
roc_curves_mis[["PY"]] <- result_PY_mis$perf
auc_scores_mis[["PY"]] <- result_PY_mis$auc_roc

# DIRICHLET MULTINOMIAL (DM)
result_DM_mis <- misclass_roc(memb_DM_indi_mis, JU, a = 1, b = 1)
roc_curves_mis[["DM"]] <- result_DM_mis$perf
auc_scores_mis[["DM"]] <- result_DM_mis$auc_roc

# GNEDIN PROCESS (GN)
result_GN_mis <- misclass_roc(memb_GN_indi_mis, JU, a = 1, b = 1)
roc_curves_mis[["GN"]] <- result_GN_mis$perf
auc_scores_mis[["GN"]] <- result_GN_mis$auc_roc

# Plot ROC curves with different colors
colors <- c("red", "green", "blue", "purple")  # Add more colors if needed
legend_text <- c("DP", "PY", "DM", "GN")  # Legend labels

# Create an empty plot
plot(roc_curves_mis[["DP"]], col = colors[1], lwd = 1)

# Add diagonal line for reference (color set to grey)
abline(a = 0, b = 1, col = "grey", lwd = 1, lty = 3)

# Iterate through the list and add each ROC curve to the plot
for (i in seq_along(roc_curves_mis)) {
  lines(roc_curves_mis[[i]]@x.values[[1]], roc_curves_mis[[i]]@y.values[[1]], col = colors[i], lwd = 1)
}

# Add a legend with AUC-ROC scores
legend("bottomright", legend = sprintf("%s AUC-ROC = %.2f", legend_text, auc_scores_mis),
       col = colors, lty = 1, cex = 0.8)

write.csv(memb_DP_indi_mis, "../ESBM ROC/memb_Z_DP_mis.csv", row.names = FALSE)
write.csv(memb_DM_indi_mis, "../ESBM ROC/memb_Z_DM_mis.csv", row.names = FALSE)
write.csv(memb_PY_indi_mis, "../ESBM ROC/memb_Z_PY_mis.csv", row.names = FALSE)
write.csv(memb_GN_indi_mis, "../ESBM ROC/memb_Z_GN_mis.csv", row.names = FALSE)
write.csv(JU, "../ESBM ROC/Y.csv", row.names = FALSE)

