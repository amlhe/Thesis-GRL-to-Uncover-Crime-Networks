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

# Ensure symmetry
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

# Make graphical visualizations
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

# Save output
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
