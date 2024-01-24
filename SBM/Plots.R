# Load the igraph and poweRlaw libraries
library(igraph)
library(poweRlaw)

# Load data
load("Adjacency matrices/adj_MN.RData")
load("Adjacency matrices/adj_PC.RData")
load("Adjacency matrices/adj_SN.RData")
load("Adjacency matrices/adj_WR.RData")
load("Adjacency matrices/adj_AW.RData")
load("Adjacency matrices/adj_JU.RData")

# Create graph objects
G_MN <- graph.adjacency(MN, mode="undirected", diag=FALSE)
G_PC <- graph.adjacency(PC, mode="undirected", diag=FALSE)
G_SN <- graph.adjacency(SN, mode="undirected", diag=FALSE)
G_WR <- graph.adjacency(WR, mode="undirected", diag=FALSE)
G_AW <- graph.adjacency(AW, mode="undirected", diag=FALSE)
G_JU <- graph.adjacency(JU, mode="undirected", diag=FALSE)

# Calculate the degree of each node and determine size for plotting
degrees_MN <- degree(G_MN)
node_sizes_MN <- degrees_MN * 1  
plot(G_MN, layout = layout_nicely, vertex.size = node_sizes_MN, 
     vertex.label = NA, edge.width = 1)

degrees_PC <- degree(G_PC)
node_sizes_PC <- degrees_PC * 1  
plot(G_PC, layout = layout_nicely, vertex.size = node_sizes_PC, 
     vertex.label = NA, edge.width = 1)

degrees_SN <- degree(G_SN)
node_sizes_SN <- degrees_SN * 1  
plot(G_SN, layout = layout_nicely, vertex.size = node_sizes_SN, 
     vertex.label = NA, edge.width = 1)

degrees_WR <- degree(G_WR)
node_sizes_WR <- degrees_WR * 1  
plot(G_WR, layout = layout_nicely, vertex.size = node_sizes_WR, 
     vertex.label = NA, edge.width = 1)

degrees_AW <- degree(G_AW)
node_sizes_AW <- degrees_AW * 1  
plot(G_AW, layout = layout_nicely, vertex.size = node_sizes_AW, 
     vertex.label = NA, edge.width = 1)

degrees_JU <- degree(G_JU)
node_sizes_JU <- degrees_JU * 1  
plot(G_JU, layout = layout_nicely, vertex.size = node_sizes_JU, 
     vertex.label = NA, edge.width = 1)



# Merge the two graphs into an aggregate network
G_Montagna <- graph.union(G_MN, G_PC)
adj_Montagna <- as_adjacency_matrix(G_Montagna, sparse = FALSE)
G_Oversize <- graph.union(G_WR, G_AW, G_JU)
adj_Oversize <- as_adjacency_matrix(G_Oversize, sparse = FALSE)

write.csv(adj_Montagna, file="Adjacency matrices/adj_Montagna.csv", row.names = FALSE)
write.csv(adj_Oversize, file="Adjacency matrices/adj_Oversize.csv", row.names = FALSE)
save(adj_Montagna,file="Adjacency matrices/adj_Montagna.RData")
save(adj_Oversize,file="Adjacency matrices/adj_Oversize.RData")

# Set the layout for better visualization
layout_Montagna <- layout_with_fr(G_Montagna)
layout_Oversize <- layout_with_fr(G_Oversize)

# Plot the aggregate networks of the two multi-layer operations
degrees_Montagna <- degree(G_Montagna)
node_sizes_Montagna <- degrees_Montagna * 1  
plot(G_Montagna, 
     vertex.label.color = "black", 
     vertex.size = 5,    
     vertex.label = NA, 
     edge.width = 1, 
     layout = layout_Montagna)

degrees_Oversize <- degree(G_Oversize)
node_sizes_Oversize <- degrees_Oversize * 1  
plot(G_Oversize,
     vertex.label.color = "black", 
     vertex.size = 5,   
     vertex.label = NA, 
     edge.width = 1, 
     layout = layout_Oversize)


### Make histogram
hist(degrees_SN, main = "Degree Distribution Histogram", xlab = "Degree", ylab = "Density", col = "lightblue", border = "black", probability = TRUE, breaks = unique(degrees_SN))
avg_degree <- mean(degrees_SN)
abline(v = avg_degree, col = "red", lty = 2, lw = 2)
legend("topright", legend = paste("Mean Degree: ", round(avg_degree, 2)), col = c("red"), lty = c(2), cex = 0.8)




# Function to plot CCDF with power law fit
plot_CCDF_with_powerlaw <- function(adj_matrix, threshold, color) {
  g <- graph.adjacency(adj_matrix, mode = "undirected", weighted = NULL)
  degrees <- degree(g)
  degrees <- degrees[degrees > 0]
  max_degree <- max(degrees)
  
  # Calculate the CCDF
  k_values <- 0:max_degree
  ccdf_values <- 1 - cumsum(table(factor(degrees, levels = k_values))) / vcount(g)
  
  # Fit the power law
  mpl <- displ$new(degrees)
  estimate <- estimate_xmin(mpl)
  mpl$setXmin(estimate)
  
  # Plot the CCDF and power law fit on a log-log scale
  lines(k_values, ccdf_values, type = "l", col = color)
  lines(mpl, col = color)
}

# Example adjacency matrices
adj_matrices <- list(MN, PC, SN, WR, AW, JU)

# Example thresholds
thresholds <- c(0.5, 0.6)

# Plot multiple CCDFs with power law fit
plot(NULL, xlim = c(1, 1000), ylim = c(0.0001, 1), 
     xlab = "k", ylab = "P(X > k)", log = "xy",
     main = "CCDF with Power Law Fit for Multiple Matrices")

colors <- c("red", "blue", "green", "orange", "purple", "yellow") # Add more colors if necessary
for (i in 1:length(adj_matrices)) {
  plot_CCDF_with_powerlaw(adj_matrices[[i]], thresholds[i], colors[i])
}
legend("topright", legend = c("Power Law Fit", "CCDF"), col = colors, lty = 1)
