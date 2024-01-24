rm(list=ls())
source("Source/esbm.R")

load("Adjacency matrices/adj_MN.RData")
load("Adjacency matrices/adj_PC.RData")
load("Adjacency matrices/adj_SN.RData")
load("Adjacency matrices/adj_WR.RData")
load("Adjacency matrices/adj_AW.RData")
load("Adjacency matrices/adj_JU.RData")

Y <- SN


# ------------------------------------
# LOUVAIN ALGORITHM
# ------------------------------------
net <- graph_from_adjacency_matrix(Y, mode=c("undirected"), weighted=NULL, diag=FALSE)
# point estimate
Louv <- cluster_louvain(net)$membership
# estimated H
H_louv <- length(table(Louv))
# deviance (D)
dev_lou <- -2*log_pY_z(Y,Louv,1,1)
# Plot reordered adjacency matrix
par(mai = c(0.1, 0.1, 0.1, 0.1)) 
plot_reordered_adjacency(Louv, Y)


############################################################
# Find number of groups for spectral clustering
############################################################

set.seed(1)
H_select <- rep(0,8)

# Le and Levina (2015)
bhmc <- BHMC.estimate(Y,K.max=20)
H_select[1] <- bhmc$K

# Wang and Bickel (2017)
lrbic <- LRBIC(Y,Kmax=20)
H_select[2] <- lrbic$SBM.K

# Chen and Lei (2018)
ncv <- NCV.select(Y,max.K=20)
H_select[3] <- which.min(ncv$l2)
H_select[4] <- which.min(ncv$dev)

# Li et al. (2020)
ecv <- ECV.block(Y,max.K=20)
H_select[5] <- which.min(ecv$l2)
H_select[6] <- which.min(ecv$dev)

# Li et al. (2020)
ecv.R <- ECV.Rank(Y,20,weighted=FALSE,mode="undirected")
H_select[7] <- ecv.R$sse.rank
H_select[8] <- ecv.R$auc.rank

sel_H <- round(median(H_select))


# ------------------------------------
# SPECTRAL CLUSTERING
# ------------------------------------
set.seed(1)
# point estimate
sc <- reg.SP(Y,K=sel_H,lap=TRUE,tau=1)$cluster
# estimated H
H_sc <- length(table(sc))
# deviance (D)
dev_sc <- -2*log_pY_z(Y,sc,1,1)
# Plot reordered adjacency matrix
plot_reordered_adjacency(sc, Y)


# ------------------------------------
# GREED SBM
# ------------------------------------
set.seed(1)
greed_out <- greed(Y, model=Sbm(alpha = 8, a0 = 1, b0 = 1, type = "undirected"))
# point estimate
g_sbm <- greed_out@cl
# estimated H
H_greed <- length(table(g_sbm))
# deviance (D)
dev_greed <- -2*log_pY_z(Y,g_sbm,1,1)
# Plot reordered adjacency matrix
plot_reordered_adjacency(g_sbm, Y)


# ------------------------------------
# GREED DC-SBM
# ------------------------------------
set.seed(1)
greed_dc <- greed(Y, model=DcSbm(alpha = 8, type = "undirected"))
# point estimate
g_dcsbm <- greed_dc@cl
# estimated H
H_greed_dc <- length(table(g_dcsbm))
# deviance (D)
dev_greed_sc <- -2*log_pY_z(Y,g_dcsbm,1,1)
# Plot reordered adjacency matrix
plot_reordered_adjacency(g_dcsbm, Y)


############################################################
# Plotting
############################################################
plot_reordered_adjacency(Louv, Y)
plot_reordered_adjacency(sc, Y)
plot_reordered_adjacency(g_sbm, Y)
plot_reordered_adjacency(g_dcsbm, Y)


layout_matrix <- matrix(c(1, 2, 3, 4), nrow = 1, byrow = TRUE)

# Set up the layout
layout(layout_matrix)

plot_reordered_adjacency(Louv, Y)
plot_reordered_adjacency(sc, Y)
plot_reordered_adjacency(g_sbm, Y)
plot_reordered_adjacency(g_dcsbm, Y)



############################################################
# Plot adj from HM-LDM
############################################################
hmldm <- read.csv("Community from HM-LDM/community_JU_D10_p2.csv")
par(mai = c(0.15, 0.15, 0.15, 0.15))
idx = which(hmldm$f_z==0)
plot_reordered_adjacency_hmldm(hmldm$z_idx[-idx], hmldm$f_z[-idx], JU)
