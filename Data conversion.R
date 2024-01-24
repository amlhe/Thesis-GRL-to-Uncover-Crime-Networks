#### Script to convert data
library(igraph)
library(Matrix)

load("Adjacency matrices/adj_MN.RData")
load("Adjacency matrices/adj_PC.RData")
load("Adjacency matrices/adj_SN.RData")
load("Adjacency matrices/adj_WR.RData")
load("Adjacency matrices/adj_AW.RData")
load("Adjacency matrices/adj_JU.RData")

load("Adjacency matrices/adj_MN_sub.RData")
load("Adjacency matrices/adj_PC_sub.RData")

G_MN <- graph_from_adjacency_matrix(MN, mode = "undirected", diag = FALSE)
MN_edge <- as_edgelist(G_MN)
write.table(MN_edge, "../HM-LDM-main/MN_edge.txt", col.names = FALSE, row.names = FALSE)

G_PC <- graph_from_adjacency_matrix(PC, mode = "undirected", diag = FALSE)
PC_edge <- as_edgelist(G_PC)
write.table(PC_edge, "../HM-LDM-main/PC_edge.txt", col.names = FALSE, row.names = FALSE)

G_SN <- graph_from_adjacency_matrix(SN, mode = "undirected", diag = FALSE)
SN_edge <- as_edgelist(G_SN)
write.table(SN_edge, "../HM-LDM-main/SN_edge.txt", col.names = FALSE, row.names = FALSE)

G_WR <- graph_from_adjacency_matrix(WR, mode = "undirected", diag = FALSE)
WR_edge <- as_edgelist(G_WR)
write.table(WR_edge, "../HM-LDM-main/WR_edge.txt", col.names = FALSE, row.names = FALSE)

G_AW <- graph_from_adjacency_matrix(AW, mode = "undirected", diag = FALSE)
AW_edge <- as_edgelist(G_AW)
write.table(AW_edge, "../HM-LDM-main/AW_edge.txt", col.names = FALSE, row.names = FALSE)

G_JU <- graph_from_adjacency_matrix(JU, mode = "undirected", diag = FALSE)
JU_edge <- as_edgelist(G_JU)
write.table(JU_edge, "../HM-LDM-main/JU_edge.txt", col.names = FALSE, row.names = FALSE)

G_MN_sub <- graph_from_adjacency_matrix(MN_sub, mode = "undirected", diag = FALSE)
MN_sub_edge <- as_edgelist(G_MN_sub)
write.table(MN_sub_edge, "../HM-LDM-main/MN_sub_edge.txt", col.names = FALSE, row.names = FALSE)

G_PC_sub <- graph_from_adjacency_matrix(PC_sub, mode = "undirected", diag = FALSE)
PC_sub_edge <- as_edgelist(G_PC_sub)
write.table(PC_sub_edge, "../HM-LDM-main/PC_sub_edge.txt", col.names = FALSE, row.names = FALSE)
