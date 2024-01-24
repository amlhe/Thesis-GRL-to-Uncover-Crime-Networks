rm(list=ls())
source("Source/esbm.R")

library(greed)
library(igraph)
library(Matrix)
library(blockmodels)
library(pheatmap)

# Read in the data
WR <- read.csv(file="Application/netWR.csv", header=FALSE, sep = ";")
AW <- read.csv(file="Application/netAW.csv", header=FALSE, sep = ";")
JU <- read.csv(file="Application/netJU.csv", header=FALSE, sep = ";")

WR <- as.matrix(WR)
AW <- as.matrix(AW)
JU <- as.matrix(JU)

WR_sym <- WR | t(WR)
WR_num <- as.numeric(WR_sym)
WR <- matrix(WR_num, nrow = nrow(WR_sym), ncol = ncol(WR_sym))
isSymmetric(WR)

AW_sym <- AW | t(AW)
AW_num <- as.numeric(AW_sym)
AW <- matrix(AW_num, nrow = nrow(AW_sym), ncol = ncol(AW_sym))
isSymmetric(AW)

JU_sym <- JU | t(JU)
JU_num <- as.numeric(JU_sym)
JU <- matrix(JU_num, nrow = nrow(JU_sym), ncol = ncol(JU_sym))
isSymmetric(JU)

heatmap(WR, 
        col = colorRampPalette(c("white", "black"))(256),  # Choose a color palette
        scale = "none",  # Use "none" to scale the colors
        Rowv = NA, Colv = NA,  # Turn off row and column clustering
        labRow = "", labCol = "",
        main = "Adjacency Matrix of WR")

heatmap(AW, 
        col = colorRampPalette(c("white", "black"))(256),  # Choose a color palette
        scale = "none",  # Use "none" to scale the colors
        Rowv = NA, Colv = NA,  # Turn off row and column clustering
        labRow = "", labCol = "",
        main = "Adjacency Matrix of AW")

heatmap(JU, 
        col = colorRampPalette(c("white", "black"))(256),  # Choose a color palette
        scale = "none",  # Use "none" to scale the colors
        Rowv = NA, Colv = NA,  # Turn off row and column clustering
        labRow = "", labCol = "",
        main = "Adjacency Matrix of JU")

save(WR, file="Adjacency matrices/adj_WR.RData")
save(AW, file="Adjacency matrices/adj_AW.RData")
save(JU, file="Adjacency matrices/adj_JU.RData")

write.csv(WR, file="Adjacency matrices/adj_WR.csv", row.names = FALSE)
write.csv(AW, file="Adjacency matrices/adj_AW.csv", row.names = FALSE)
write.csv(JU, file="Adjacency matrices/adj_JU.csv", row.names = FALSE)




# Transform adjacency matrix into an igraph object
net_WR <- graph.adjacency(WR, mode=c("undirected"), weighted=NULL, diag=FALSE)
net_AW <- graph.adjacency(AW, mode=c("undirected"), weighted=NULL, diag=FALSE)
net_JU <- graph.adjacency(AW, mode=c("undirected"), weighted=NULL, diag=FALSE)


# Point estimate
Louv_WR <- cluster_louvain(net_WR)$membership
lou_WR <- cluster_louvain(net_WR)
Louv_AW <- cluster_louvain(net_AW)$membership
lou_AW <- cluster_louvain(net_AW)
Louv_JU <- cluster_louvain(net_JU)$membership
lou_JU <- cluster_louvain(net_JU)

# Estimated H
length(table(Louv_WR))
length(table(Louv_AW))
length(table(Louv_JU))

# Plot communities
plot(lou_WR,net_WR,vertex.label = NA)
plot(lou_AW,net_AW,vertex.label = NA)
plot(lou_JU,net_JU,vertex.label = NA)


Y <- WR
V <- dim(Y)[1]


# ------------------------------------
# DIRICHLET MULTINOMIAL
# ------------------------------------
sigma_dm <- 0   
H_dm <- 50 # Conservative upper bound 
beta_dm <- 12/H_dm 
round(expected_cl_py(V, sigma = sigma_dm, theta = beta_dm*H_dm, H = H_dm))



N_iter <- 50000
V <- dim(Y)[1]
my_seed <- 1
my_z <- c(1:V)

# ------------------------------------
# DIRICHLET MULTINOMIAL
# ------------------------------------

my_prior <- "DM"
Z_DM <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, beta_DM = 12/50, H_DM = 50)




set.seed(1)
index_traceplot <- sample(c(1:(V*(V-1)/2)),1)

set.seed(1)
V <- dim(Y)[1]
burn_in <- 10000
a <- b <- 1
LL <- matrix(nrow=V*(V-1)/2,ncol=40000)

# ------------------------------------
# DIRICHLET MULTINOMIAL UNSUPERVISED
# ------------------------------------
Z_DM_WAIC <- Z_DM[,(burn_in+1):N_iter]


for (t in 1:dim(Z_DM_WAIC)[2]){
  LL[,t]<-sampleLL(Z_DM_WAIC[,t],Y,a,b)
  if (t%%10000 == 0){print(paste("Iteration:", t))}
}
WAIC(LL)$WAIC
# Selected traceplot
plot(ts(LL[index_traceplot,]),xlab="",ylab="",ylim=c(-4,0))



# ------------------------------------
# DIRICHLET MULTINOMIAL
# ------------------------------------
quantile(apply(Z_DM[,(burn_in+1):N_iter],2,max))[c(2:4)]




# ------------------------------------
# DIRICHLET MULTINOMIAL
# ------------------------------------

c_Z_DM <- pr_cc(Z_DM[,(burn_in+1):N_iter])

# point estimate
memb_Z_DM_VI <- minVI(c_Z_DM,method="avg",max.k=20)
memb_Z_DM <- memb_Z_DM_VI$cl

# horizontal bound of the credible ball
credibleball(memb_Z_DM_VI$cl,t(Z_DM[,(burn_in+1):N_iter]))[[5]]

# misclassification error [in-sample for edges]
misclass(memb_Z_DM,Y,a=1,b=1)



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
sc <- reg.SP(Y,K=sel_H,lap=TRUE,tau=0)$cluster

# estimated H
length(table(sc))

# deviance (D)
-2*log_pY_z(Y,sc,1,1)

# ------------------------------------
# REGULARIZED SPECTRAL CLUSTERING
# ------------------------------------

set.seed(1)

# point estimate
r_sc <- reg.SP(Y,K=sel_H,lap=TRUE,tau=1)$cluster

# estimated H
length(table(r_sc))

# deviance (D)
-2*log_pY_z(Y,r_sc,1,1)


# ------------------------------------
# GREED SBM
# ------------------------------------

set.seed(1)

graph <- igraph::graph.adjacency(Y, mode = "undirected")
preprocessed_graph <- igraph::preprocess(graph)
sbm_object <- sbm(alpha = 12/50)
result <- greed(graph, K = sel_H, model = sbm_object, alg = methods::new("hybrid"), verbose = FALSE)

greed_out <- greed(Y,K=sel_H,model=new("sbm",alpha=12/50,type="undirected"),alg=methods::new("hybrid"),verbose=FALSE)

# point estimate
g_sbm <- greed_out@cl

# estimated H
length(table(g_sbm))

# deviance (D)
-2*log_pY_z(Y,g_sbm,1,1)




# transform the adjacency matrix into an igraph object
net_Y <- graph.adjacency(Y, mode=c("undirected"), weighted=NULL, diag=FALSE)

# compute the node betweenness to be used for the size of the nodes
betw <- betweenness(net_Y)
# node sizes are proportional to their betweenness 
# Note: for graphical purposes, we consider a monotone transformation of such a measure
V(net_Y)$size <- sqrt(betw/1.5+mean(betw))*2

# additional graphical settings
V(net_Y)$frame.color <- "black"
V(net_Y)$label <- "" 
E(net_Y)$color <- brewer.pal(9,"Greys")[3]

# node positions are obtained via force–directed placement
set.seed(12)
l <- layout_with_fr(net_Y)
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1.5, xmax=1.5)

# plot Figure 1
plot(net_Y, rescale=F, layout=l*1.5,edge.curved=.3,edge.width=0.5)




# ------------------------------------
# LOUVAIN ALGORITHM
# ------------------------------------

# transform the adjacency matrix into an igraph object
net <- graph.adjacency(Y, mode=c("undirected"), weighted=NULL, diag=FALSE)

# point estimate
Louv <- cluster_louvain(net)$membership

# to display the block structures, re-order the rows and columns of Y, and the elements 
# in RoleLocale according to the groupings estimated by the Louvain algorithm
sel <- order(Louv_WR)
Louv <- Louv_WR[sel]
Y_Louv <- Y[sel,sel]
RoleLocale_Louv <- RoleLocale[sel]

# plot the adjacency with the grouping structure defined by the Louvain algorithm
row_plotLouv <- as.data.frame(as.factor(matrix(RoleLocale_Louv,V,1)))
names(row_plotLouv) <- "RoleLocale_Louv"
rownames(Y_Louv) <- rownames(row_plotLouv)
names(mycolors) <- sort(unique(row_plotLouv$RoleLocale_Louv))

Adj_Louv <- pheatmap(Y_Louv,cluster_cols = F, cluster_rows= F,annotation_names_row=F, show_rownames=F,show_colnames=F,legend=F,border_color=FALSE, annotation_legend=F, annotation_colors=list(RoleLocale_Louv = mycolors),gaps_row=c(which(diff(Louv)!=0)),gaps_col=c(which(diff(Louv)!=0)))

