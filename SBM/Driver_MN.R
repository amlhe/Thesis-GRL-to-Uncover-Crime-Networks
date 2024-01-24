rm(list=ls())
source("Source/esbm.R")

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
library(RColorBrewer)
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


load("Adjacency matrices/adj_MN.RData")
load("Adjacency matrices/adj_PC.RData")
load("Adjacency matrices/adj_SN.RData")
load("Adjacency matrices/adj_WR.RData")
load("Adjacency matrices/adj_AW.RData")
load("Adjacency matrices/adj_JU.RData")


### Plot initial adjacency matrix
# Plot the reordered adjacency matrix using the image function
par(mai = c(0.1, 0.1, 0.1, 0.1)) 
image(1:nrow(MN), 1:ncol(MN), MN, 
      col = colorRampPalette(brewer.pal(9, "Greys")[c(1, 8)])(30),
      axes = FALSE, xlab = "", ylab = "", asp = 1, xlim = c(0, nrow(MN)+1), ylim = c(0,ncol(MN)+1))
rect(0, 0, ncol(MN)+1, nrow(MN)+1, border = "black", lwd = 2)
fun_grob <- function(){
  grid.echo()
  grid.grab()
}
fun_grob_out <- fun_grob()
pushViewport(viewport(width = 1, height = 1, angle = 270))
grid.draw(fun_grob_out)

# Missing links for link prediction
sparse_i_MN <- read.table("Missing links/MN/sparse_i.txt", header = FALSE)$V1
sparse_j_MN <- read.table("Missing links/MN/sparse_j.txt", header = FALSE)$V1
num_nodes <- 101
MN_mis <- matrix(0, nrow = num_nodes, ncol = num_nodes)
for (i in 1:length(sparse_i_MN)) {
  MN_mis[sparse_i_MN[i], sparse_j_MN[i]] <- 1
  MN_mis[sparse_j_MN[i], sparse_i_MN[i]] <- 1
}


set.seed(1)  # Set seed for reproducibility

### Create igraph object from adjacency matrix
G_MN <- graph.adjacency(MN, mode="undirected", diag=FALSE)
G_MN_mis <- graph.adjacency(MN_mis, mode="undirected", diag=FALSE)

### ESBM
V <- dim(MN)[1]
#V <- dim(PC)[1]
#V <- dim(SN)[1]
#V <- dim(WR)[1]
#V <- dim(AW)[1]
#V <- dim(JU)[1]

Y <- MN

# ------------------------------------
# DIRICHLET PROCESS (CRP)
sigma_dp <- 0  
H_dp <- Inf 
alpha_dp <- 8
round(expected_cl_py(V, sigma = sigma_dp, theta = alpha_dp, H = H_dp))


# ------------------------------------
# PITMAN-YOR (PY)
sigma_py <- 0.7
H_py <- Inf 
alpha_py <- -0.350
round(expected_cl_py(V, sigma = sigma_py, theta = alpha_py, H = H_py))


# DIRICHLET MULTINOMIAL
# ------------------------------------
sigma_dm <- 0   
H_dm <- 50 # Conservative upper bound 
beta_dm <- 12/H_dm 
round(expected_cl_py(V, sigma = sigma_dm, theta = beta_dm*H_dm, H = H_dm))


# GNEDIN PROCESS
# ------------------------------------
gamma <- 0.3
probs_gnedin <- HGnedin(V, 1:V, gamma = gamma)
round(sum(1:V*probs_gnedin))



############################################################
# Posterior computation using collapsed Gibbs sampler
############################################################

N_iter <- 50000
V <- dim(Y)[1]
my_seed <- 1
my_z <- c(1:V)


# DIRICHLET MULTINOMIAL
# ------------------------------------
my_prior <- "DM"
Z_DM <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, beta_DM = 12/50, H_DM = 50)
Z_DM_mis <- esbm(MN_mis, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, beta_DM = 12/50, H_DM = 50)

# DIRICHLET PROCESS (CRP)
# ------------------------------------
my_prior <- "DP"
Z_DP <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, beta_DM = 12/50, H_DM = 50)
Z_DP_mis <- esbm(MN_mis, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, alpha_PY = 8, sigma_PY = 0)

# PITMAN-YOR PROCESS
# ------------------------------------
my_prior <- "PY"
Z_PY <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, beta_DM = 12/50, H_DM = 50)
Z_PY_mis <- esbm(MN_mis, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, alpha_PY = -0.350, sigma_PY = 0.725)

# GNEDIN PROCESS
# ------------------------------------
my_prior  <- "GN"
Z_GN <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, beta_DM = 12/50, H_DM = 50)
Z_GN_mis <- esbm(MN_mis, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, gamma_GN = 0.3)


# Save the output
save(Z_DP,Z_PY,Z_GN,Z_DM,file="Application/Posterior_No_Attributes_MN_fixed.RData")
load("Application/Posterior_No_Attributes_MN_fixed.RData")

save(Z_DP_mis,Z_PY_mis,Z_GN_mis,Z_DM_mis,file="Application/Posterior_No_Attributes_MN_mis.RData")
load("Application/Posterior_No_Attributes_MN_mis.RData")


############################################################
# Posterior inference
############################################################
set.seed(1)
index_traceplot <- sample(c(1:(V*(V-1)/2)),1)


V <- dim(Y)[1]
burn_in <- 10000
a <- b <- 1
LL <- matrix(nrow=V*(V-1)/2,ncol=40000)



### DIRICHLET PROCESS
Z_DP_WAIC <- Z_DP[,(burn_in+1):N_iter]
for (t in 1:dim(Z_DP_WAIC)[2]){
  LL[,t]<-sampleLL(Z_DP_WAIC[,t],Y,a,b)
  if (t%%1000 == 0){print(paste("Iteration:", t))}
}
WAIC(LL)$WAIC
# Traceplot
plot(ts(LL[index_traceplot,]),xlab="",ylab="",ylim=c(-0.2,0.1), main = "Dirichlet process")



### DIRICHLET MULTINOMIAL
Z_DM_WAIC <- Z_DM[,(burn_in+1):N_iter]
for (t in 1:dim(Z_DM_WAIC)[2]){
  LL[,t]<-sampleLL(Z_DM_WAIC[,t],Y,a,b)
  if (t%%1000 == 0){print(paste("Iteration:", t))}
}
WAIC(LL)$WAIC
# Traceplot
plot(ts(LL[index_traceplot,]),xlab="",ylab="",ylim=c(-0.2,0.1), main = "Dirichlet multinomial")



### PITMAN-YOR
Z_PY_WAIC <- Z_PY[,(burn_in+1):N_iter]
for (t in 1:dim(Z_PY_WAIC)[2]){
  LL[,t]<-sampleLL(Z_PY_WAIC[,t],Y,a,b)
  if (t%%1000 == 0){print(paste("Iteration:", t))}
}
WAIC(LL)$WAIC
# Traceplot
plot(ts(LL[index_traceplot,]),xlab="",ylab="",ylim=c(-0.2,0.1), main = "Pitman-Yor process")



### GNEDIN PROCESS
Z_GN_WAIC <- Z_GN[,(burn_in+1):N_iter]
for (t in 1:dim(Z_GN_WAIC)[2]){
  LL[,t]<-sampleLL(Z_GN_WAIC[,t],Y,a,b)
  if (t%%1000 == 0){print(paste("Iteration:", t))}
}
WAIC(LL)$WAIC
# Traceplot
plot(ts(LL[index_traceplot,]),xlab="",ylab="",ylim=c(-0.2,0.1), main = "Gnedin process")

### QUANTILES FROM DIFFERENT PRIORS
quantile(apply(Z_DP[,(burn_in+1):N_iter],2,max))[c(2:4)]
quantile(apply(Z_DM[,(burn_in+1):N_iter],2,max))[c(2:4)]
quantile(apply(Z_PY[,(burn_in+1):N_iter],2,max))[c(2:4)]
quantile(apply(Z_GN[,(burn_in+1):N_iter],2,max))[c(2:4)]

# ------------------------------------
### Point estimates and credible balls
### DIRICHLET PROCESS
c_Z_DP <- pr_cc(Z_DP[,(burn_in+1):N_iter])

# Point estimate
memb_Z_DP_VI <- minVI(c_Z_DP,method="avg",max.k=21)
memb_Z_DP <- memb_Z_DP_VI$cl

# Credible ball
crediblecall_DP <- credibleball(memb_Z_DP_VI$cl,t(Z_DP[,(burn_in+1):N_iter]))[[5]]

# Misclassification error 
misclass_DP <- misclass(memb_Z_DP,Y,a=1,b=1)

# Deviance
-2*log_pY_z(Y,memb_Z_DP,1,1)

### DIRICHLET MULTINOMIAL
c_Z_DM <- pr_cc(Z_DM[,(burn_in+1):N_iter])

# Point estimate
memb_Z_DM_VI <- minVI(c_Z_DM,method="avg",max.k=21)
memb_Z_DM <- memb_Z_DM_VI$cl

# Credible ball
crediblecall_DM <- credibleball(memb_Z_DM_VI$cl,t(Z_DM[,(burn_in+1):N_iter]))[[5]]

# Misclassification error
misclass_DM <- misclass(memb_Z_DM,Y,a=1,b=1)

# Deviance
-2*log_pY_z(Y,memb_Z_DM,1,1)


### PITMAN-YOR 
c_Z_PY <- pr_cc(Z_PY[,(burn_in+1):N_iter])

# Point estimate
memb_Z_PY_VI <- minVI(c_Z_PY,method="avg",max.k=20)
memb_Z_PY <- memb_Z_PY_VI$cl

# Credible ball
crediblecall_PY <- credibleball(memb_Z_PY_VI$cl,t(Z_PY[,(burn_in+1):N_iter]))[[5]]

# Misclassification error
misclass_PY <- misclass(memb_Z_PY,Y,a=1,b=1)

# Deviance
-2*log_pY_z(Y,memb_Z_PY,1,1)


### GNEDIN PROCESS
c_Z_GN <- pr_cc(Z_GN[,(burn_in+1):N_iter])

# Point estimate
memb_Z_GN_VI <- minVI(c_Z_GN,method="avg",max.k=23)
memb_Z_GN <- memb_Z_GN_VI$cl

# Credible ball
crediblecall_GN <- credibleball(memb_Z_GN_VI$cl,t(Z_GN[,(burn_in+1):N_iter]))[[5]]

# Misclassification error
misclass_GN <- misclass(memb_Z_GN,Y,a=1,b=1)

# Deviance
-2*log_pY_z(Y,memb_Z_GN,1,1)

# Plot reordered adjacency matrix
par(mai = c(0.1, 0.1, 0.1, 0.1))
plot_reordered_adjacency(memb_Z_DP, Y)
plot_reordered_adjacency(memb_Z_DM, Y)
plot_reordered_adjacency(memb_Z_PY, Y)
plot_reordered_adjacency(memb_Z_GN, Y)


# ------------------------------------
# ROC collection form different priors
# Collect ROC curves and AUC-ROC scores for all priors
roc_curves <- list()
auc_scores <- numeric()

# DIRICHLET PROCESS (CRP) UNSUPERVISED
result_DP <- misclass_roc(memb_Z_DP, Y, a = 1, b = 1)
roc_curves[["DP"]] <- result_DP$perf
auc_scores[["DP"]] <- result_DP$auc_roc

# PITMAN-YOR (PY)
result_PY <- misclass_roc(memb_Z_PY, Y, a = 1, b = 1)
roc_curves[["PY"]] <- result_PY$perf
auc_scores[["PY"]] <- result_PY$auc_roc

# DIRICHLET MULTINOMIAL (DM)
result_DM <- misclass_roc(memb_Z_DM, Y, a = 1, b = 1)
roc_curves[["DM"]] <- result_DM$perf
auc_scores[["DM"]] <- result_DM$auc_roc

# GNEDIN PROCESS (GN)
result_GN <- misclass_roc(memb_Z_GN, Y, a = 1, b = 1)
roc_curves[["GN"]] <- result_GN$perf
auc_scores[["GN"]] <- result_GN$auc_roc

# Plot ROC curves with different colors
colors <- c("red", "green", "blue", "purple")  # Add more colors if needed
legend_text <- c("DP", "PY", "DM", "GN")  # Legend labels

# Create an empty plot
plot(roc_curves[["DP"]], col = colors[1], lwd = 1)

# Add diagonal line for reference (color set to grey)
abline(a = 0, b = 1, col = "grey", lwd = 1, lty = 3)

# Iterate through the list and add each ROC curve to the plot
for (i in seq_along(roc_curves)) {
  lines(roc_curves[[i]]@x.values[[1]], roc_curves[[i]]@y.values[[1]], col = colors[i], lwd = 1)
}

# Add a legend with AUC-ROC scores
legend("bottomright", legend = sprintf("%s AUC-ROC = %.2f", legend_text, auc_scores),
       col = colors, lty = 1, cex = 0.8)



############################################################
# Posterior inference - missing links
############################################################
set.seed(1)
index_traceplot <- sample(c(1:(V*(V-1)/2)),1)

V <- dim(Y)[1]
burn_in <- 10000
a <- b <- 1
LL <- matrix(nrow=V*(V-1)/2,ncol=40000)


### DIRICHLET PROCESS
Z_DP_WAIC_mis <- Z_DP_mis[,(burn_in+1):N_iter]
for (t in 1:dim(Z_DP_WAIC_mis)[2]){
  LL[,t]<-sampleLL(Z_DP_WAIC_mis[,t],MN_mis,a,b)
  if (t%%1000 == 0){print(paste("Iteration:", t))}
}

### DIRICHLET MULTINOMIAL
Z_DM_WAIC_mis <- Z_DM_mis[,(burn_in+1):N_iter]
for (t in 1:dim(Z_DM_WAIC_mis)[2]){
  LL[,t]<-sampleLL(Z_DM_WAIC_mis[,t],MN_mis,a,b)
  if (t%%1000 == 0){print(paste("Iteration:", t))}
}

### PITMAN-YOR
Z_PY_WAIC_mis <- Z_PY_mis[,(burn_in+1):N_iter]
for (t in 1:dim(Z_PY_WAIC_mis)[2]){
  LL[,t]<-sampleLL(Z_PY_WAIC_mis[,t],MN_mis,a,b)
  if (t%%1000 == 0){print(paste("Iteration:", t))}
}

### GNEDIN PROCESS
Z_GN_WAIC_mis <- Z_GN_mis[,(burn_in+1):N_iter]
for (t in 1:dim(Z_GN_WAIC_mis)[2]){
  LL[,t]<-sampleLL(Z_GN_WAIC_mis[,t],MN_mis,a,b)
  if (t%%1000 == 0){print(paste("Iteration:", t))}
}

# ------------------------------------
### Point estimates and credible balls
### DIRICHLET PROCESS
c_Z_DP_mis <- pr_cc(Z_DP_mis[,(burn_in+1):N_iter])

# Point estimate
memb_Z_DP_VI_mis <- minVI(c_Z_DP_mis,method="avg",max.k=21)
memb_Z_DP_mis <- memb_Z_DP_VI_mis$cl

### DIRICHLET MULTINOMIAL
c_Z_DM_mis <- pr_cc(Z_DM_mis[,(burn_in+1):N_iter])

# Point estimate
memb_Z_DM_VI_mis <- minVI(c_Z_DM_mis,method="avg",max.k=21)
memb_Z_DM_mis <- memb_Z_DM_VI_mis$cl

### PITMAN-YOR 
c_Z_PY_mis <- pr_cc(Z_PY_mis[,(burn_in+1):N_iter])

# Point estimate
memb_Z_PY_VI_mis <- minVI(c_Z_PY_mis,method="avg",max.k=20)
memb_Z_PY_mis <- memb_Z_PY_VI_mis$cl

### GNEDIN PROCESS
c_Z_GN_mis <- pr_cc(Z_GN_mis[,(burn_in+1):N_iter])

# Point estimate
memb_Z_GN_VI_mis <- minVI(c_Z_GN_mis,method="avg",max.k=23)
memb_Z_GN_mis <- memb_Z_GN_VI_mis$cl


# ------------------------------------
# ROC collection form different priors
# Collect ROC curves and AUC-ROC scores for all priors
roc_curves_mis <- list()
pr_curves_mis <- list()
auc_roc_scores_mis <- numeric()
auc_pr_scores_mis <- numeric()

# DIRICHLET PROCESS (CRP) UNSUPERVISED
result_DP_mis <- misclass_roc(memb_Z_DP_mis, Y, a = 1, b = 1)
roc_curves_mis[["DP"]] <- result_DP_mis$perf
pr_curves_mis[["DP"]] <- result_DP_mis$pr_curve
auc_roc_scores_mis[["DP"]] <- result_DP_mis$auc_roc
auc_pr_scores_mis[["DP"]] <- result_DP_mis$auc_pr

# PITMAN-YOR (PY)
result_PY_mis <- misclass_roc(memb_Z_PY_mis, Y, a = 1, b = 1)
roc_curves_mis[["PY"]] <- result_PY_mis$perf
pr_curves_mis[["PY"]] <- result_PY_mis$pr_curve
auc_roc_scores_mis[["PY"]] <- result_PY_mis$auc_roc
auc_pr_scores_mis[["PY"]] <- result_PY_mis$auc_pr

# DIRICHLET MULTINOMIAL (DM)
result_DM_mis <- misclass_roc(memb_Z_DM_mis, Y, a = 1, b = 1)
roc_curves_mis[["DM"]] <- result_DM_mis$perf
pr_curves_mis[["DM"]] <- result_DM_mis$pr_curve
auc_roc_scores_mis[["DM"]] <- result_DM_mis$auc_roc
auc_pr_scores_mis[["DM"]] <- result_DM_mis$auc_pr

# GNEDIN PROCESS (GN)
result_GN_mis <- misclass_roc(memb_Z_GN_mis, Y, a = 1, b = 1)
roc_curves_mis[["GN"]] <- result_GN_mis$perf
pr_curves_mis[["GN"]] <- result_GN_mis$pr_curve
auc_roc_scores_mis[["GN"]] <- result_GN_mis$auc_roc
auc_pr_scores_mis[["GN"]] <- result_GN_mis$auc_pr

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
legend("bottomright", legend = sprintf("%s AUC-ROC = %.2f", legend_text, auc_roc_scores_mis),
       col = colors, lty = 1, cex = 0.8)


# Plot PR curves with different colors
colors <- c("red", "green", "blue", "purple")  # Add more colors if needed
legend_text <- c("DP", "PY", "DM", "GN")  # Legend labels

# Create an empty plot
plot(pr_curves_mis[["DP"]], col = colors[1], lwd = 1)

# Add diagonal line for reference (color set to grey)
abline(a = 0, b = 1, col = "grey", lwd = 1, lty = 3)

# Iterate through the list and add each ROC curve to the plot
for (i in seq_along(pr_curves_mis)) {
  lines(pr_curves_mis[[i]]@x.values[[1]], pr_curves_mis[[i]]@y.values[[1]], col = colors[i], lwd = 1)
}

# Add a legend with AUC-ROC scores
legend("bottomright", legend = sprintf("%s AUC-PR = %.2f", legend_text, auc_pr_scores_mis),
       col = colors, lty = 1, cex = 0.8)


# ------------------------------------------------------------
# ADJACENCY AND GRAPH VISUALIZATIONS OF ESBM WITH DP PRIOR
# ------------------------------------------------------------
# Get cluster labels
memb_Z_DP <- memb_Z_DP_VI$cl
group_order <- unique(memb_Z_DP)

# Re-order the rows and columns of Y
sel <- which(memb_Z_DP==group_order[1])
for (k in 2:length(group_order)){
  sel <- c(sel,which(memb_Z_DP==group_order[k]))	
}

memb_Z_DP_sel <- memb_Z_DP[sel]
Y_esbm_DP <- Y[sel,sel]

# plot the adjacency with the grouping structure defined by ESBM with unsupervised DP prior
Adj_esbm <- pheatmap(Y_esbm_DP,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F, annotation_names_row=F,show_rownames=F, show_colnames=F, legend=F ,border_color=FALSE, annotation_legend=F,gaps_row=c(which(diff(memb_Z_DP)!=0)),gaps_col=c(which(diff(memb_Z_DP)!=0)))

block_colors <- rainbow(length(unique(memb_Z_DP_sel)))
node_colors <- block_colors[memb_Z_DP_sel]
plot(G_MN, vertex.size = 10, vertex.color = node_colors, vertex.label = NA, main = paste("Number of communities:", length(unique(memb_Z_DP))))



# ------------------------------------------------------------
# ADJACENCY AND GRAPH VISUALIZATIONS OF ESBM WITH PY PRIOR
# ------------------------------------------------------------
# Get cluster labels
memb_Z_PY <- memb_Z_PY_VI$cl
group_order <- unique(memb_Z_PY)

# Re-order the rows and columns of Y
sel <- which(memb_Z_PY==group_order[1])
for (k in 2:length(group_order)){
  sel <- c(sel,which(memb_Z_PY==group_order[k]))	
}

memb_Z_PY_sel <- memb_Z_PY[sel]
Y_esbm_PY <- Y[sel,sel]

# plot the adjacency with the grouping structure defined by ESBM with unsupervised DP prior
Adj_esbm <- pheatmap(Y_esbm_PY,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F, annotation_names_row=F,show_rownames=F, show_colnames=F, legend=F ,border_color=FALSE, annotation_legend=F,gaps_row=c(which(diff(memb_Z_PY)!=0)),gaps_col=c(which(diff(memb_Z_PY)!=0)))

block_colors <- rainbow(length(unique(memb_Z_PY_sel)))
node_colors <- block_colors[memb_Z_PY_sel]
plot(G_MN, vertex.size = 10, vertex.color = node_colors, vertex.label = NA, main = paste("Number of communities:", length(unique(memb_Z_PY))))

# ------------------------------------------------------------
# ADJACENCY AND GRAPH VISUALIZATIONS OF ESBM WITH DM PRIOR
# ------------------------------------------------------------
# Get cluster labels
memb_Z_DM <- memb_Z_DM_VI$cl
group_order <- unique(memb_Z_DM)

# Re-order the rows and columns of Y
sel <- which(memb_Z_DM==group_order[1])
for (k in 2:length(group_order)){
  sel <- c(sel,which(memb_Z_DM==group_order[k]))	
}

memb_Z_DM_sel <- memb_Z_DM[sel]
Y_esbm_DM <- Y[sel,sel]

# plot the adjacency with the grouping structure defined by ESBM with unsupervised DP prior
Adj_esbm <- pheatmap(Y_esbm_DM,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F, annotation_names_row=F,show_rownames=F, show_colnames=F, legend=F ,border_color=FALSE, annotation_legend=F,gaps_row=c(which(diff(memb_Z_DM)!=0)),gaps_col=c(which(diff(memb_Z_DM)!=0)))

block_colors <- rainbow(length(unique(memb_Z_DM_sel)))
node_colors <- block_colors[memb_Z_DM_sel]
plot(G_MN, vertex.size = 10, vertex.color = node_colors, vertex.label = NA, main = paste("Number of communities:", length(unique(memb_Z_DM))))


# ------------------------------------------------------------
# ADJACENCY AND GRAPH VISUALIZATIONS OF ESBM WITH GN PRIOR
# ------------------------------------------------------------
# Get cluster labels
memb_Z_GN <- memb_Z_GN_VI$cl
group_order <- unique(memb_Z_GN)

# Re-order the rows and columns of Y
sel <- which(memb_Z_GN==group_order[1])
for (k in 2:length(group_order)){
  sel <- c(sel,which(memb_Z_GN==group_order[k]))	
}

memb_Z_GN_sel <- memb_Z_GN[sel]
Y_esbm_GN <- Y[sel,sel]

# plot the adjacency with the grouping structure defined by ESBM with unsupervised GN prior
color_palette <- colorRampPalette(brewer.pal(9, "Greys")[c(1, 8)])(30)
Adj_esbm <- pheatmap(Y_esbm_GN,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F, annotation_names_row=F,show_rownames=F, show_colnames=F, legend=F ,border_color=FALSE, annotation_legend=F,gaps_row=c(which(diff(memb_Z_GN)!=0)),gaps_col=c(which(diff(memb_Z_GN)!=0)))




block_colors <- rainbow(length(unique(memb_Z_GN_sel)))
node_colors <- block_colors[memb_Z_GN_sel]
plot(G_MN, vertex.size = 10, vertex.color = node_colors, vertex.label = NA, main = paste("Number of communities:", length(unique(memb_Z_GN))))

