rm(list=ls())
source("Source/esbm.R")

library(igraph)

# Read in the data
A <- read.csv(file="Application/Montagna_Meetings_Edgelist.csv",header=FALSE)
B <- read.csv(file="Application/Montagna_Phonecalls_Edgelist.csv",header=FALSE)

# Split the Factor_Column into three columns using space as a delimiter
A <- strsplit(A[,1], '\t')
B <- strsplit(B[,1], '\t')

# Convert the result to a data frame
A <- data.frame(do.call(rbind, A))
B <- data.frame(do.call(rbind, B))

# Rename the columns if needed
colnames(A) <- c('V1', 'V2', 'W')
colnames(B) <- c('V1', 'V2', 'W')

# Get unique persons
persons_A <- sort(unique(as.numeric(c(A$V1, A$V2))))
persons_B <- sort(unique(as.numeric(c(B$V1, B$V2))))

# Create an empty matrix
adj_A <- matrix(0, nrow = length(persons_A), ncol = length(persons_A), dimnames = list(persons_A, persons_A))
adj_B <- matrix(0, nrow = length(persons_B), ncol = length(persons_B), dimnames = list(persons_B, persons_B))

# Fill in the values and order
adj_A[cbind(A$V1, A$V2)] <- 1
adj_A <- adj_A[order(as.numeric(rownames(adj_A))), order(as.numeric(colnames(adj_A)))]
adj_B[cbind(B$V1, B$V2)] <- 1
adj_B <- adj_B[order(as.numeric(rownames(adj_B))), order(as.numeric(colnames(adj_B)))]

MN <- as.matrix(adj_A)
PC <- as.matrix(adj_B)

MN_sym <- MN | t(MN)
MN_num <- as.numeric(MN_sym)
MN <- matrix(MN_num, nrow = nrow(MN_sym), ncol = ncol(MN_sym))
isSymmetric(MN)

PC_sym <- PC | t(PC)
PC_num <- as.numeric(PC_sym)
PC <- matrix(PC_num, nrow = nrow(PC_sym), ncol = ncol(PC_sym))
isSymmetric(PC)

heatmap(MN, 
        col = colorRampPalette(c("white", "black"))(256),  # Choose a color palette
        scale = "none",  # Use "none" to scale the colors
        Rowv = NA, Colv = NA,  # Turn off row and column clustering
        labRow = "", labCol = "",
        main = "Adjacency Matrix of MN")

heatmap(PC, 
        col = colorRampPalette(c("white", "black"))(256),  # Choose a color palette
        scale = "none",  # Use "none" to scale the colors
        Rowv = NA, Colv = NA,  # Turn off row and column clustering
        labRow = "", labCol = "",
        main = "Adjacency Matrix of PC")


write.csv(MN, file="Adjacency matrices/adj_MN.csv", row.names = FALSE)
write.csv(PC, file="Adjacency matrices/adj_PC.csv", row.names = FALSE)

save(MN,file="Adjacency matrices/adj_MN.RData")
save(PC,file="Adjacency matrices/adj_PC.RData")

### Find common nodes
common_MN_PC <- intersect(persons_A, persons_B)
G_MN <- graph_from_adjacency_matrix(MN, mode = "undirected", diag = FALSE)
G_PC <- graph_from_adjacency_matrix(PC, mode = "undirected", diag = FALSE)
G_MN_sub <- induced_subgraph(G_MN, common_MN_PC)
G_PC_sub <- induced_subgraph(G_PC, common_MN_PC)

MN_sub <- MN[common_MN_PC, common_MN_PC]
PC_sub <- PC[common_MN_PC, common_MN_PC]

write.csv(MN_sub, file="Adjacency matrices/adj_MN_sub.csv", row.names = FALSE)
write.csv(PC_sub, file="Adjacency matrices/adj_PC_sub.csv", row.names = FALSE)

save(MN_sub, file="Adjacency matrices/adj_MN_sub.RData")
save(PC_sub, file="Adjacency matrices/adj_PC_sub.RData")

# Create a heatmap of the adjacency matrix
heatmap(MN, 
        col = colorRampPalette(c("white", "blue"))(256),  # Choose a color palette
        scale = "none",  # Use "none" to scale the colors
        Rowv = NA, Colv = NA,  # Turn off row and column clustering
        labRow = "", labCol = "",
        main = "Adjacency Matrix of MN")

heatmap(PC, 
        col = colorRampPalette(c("white", "blue"))(256),  # Choose a color palette
        scale = "none",  # Use "none" to scale the colors
        Rowv = NA, Colv = NA,  # Turn off row and column clustering
        labRow = "", labCol = "",
        main = "Adjacency Matrix of PC")


# Transform adjacency matrix into an igraph object
net_MN <- graph.adjacency(adj_MN_num, mode=c("undirected"), weighted=NULL, diag=FALSE)
net_PC <- graph.adjacency(adj_PC_num, mode=c("undirected"), weighted=NULL, diag=FALSE)

# Point estimate
Louv_MN <- cluster_louvain(net_MN)$membership
lou_MN <- cluster_louvain(net_MN)
Louv_PC <- cluster_louvain(net_PC)$membership
lou_PC <- cluster_louvain(net_PC)

# Estimated H
length(table(Louv_MN))
length(table(Louv_PC))

# Plot communities
plot(lou_MN,net_MN,vertex.label = NA)
plot(lou_PC,net_PC,vertex.label = NA)


Y <- adj_MN_num


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


