library(Matrix)
library(ggplot2)
library(MASS)
library(tidyr)
library(pROC)
library(igraph)

load("Adjacency matrices/adj_MN.RData")
load("Adjacency matrices/adj_PC.RData")
load("Adjacency matrices/adj_SN.RData")
load("Adjacency matrices/adj_WR.RData")
load("Adjacency matrices/adj_AW.RData")
load("Adjacency matrices/adj_JU.RData")

###########################################################################
# COMPUTE LOCAL SIMILARITY SCORES ##########################################
############################################################################

# Function to calculate Common Neighbors (CN)
common_neighbors <- function(graph, node_i, node_j) {
  neighbors_i <- which(graph[node_i, ] == 1)
  neighbors_j <- which(graph[node_j, ] == 1)
  common <- length(intersect(neighbors_i, neighbors_j))
  return(common)
}

# Function to calculate Jaccard Coefficient (JC)
jaccard_coefficient <- function(graph, node_i, node_j) {
  neighbors_i <- which(graph[node_i, ] == 1)
  neighbors_j <- which(graph[node_j, ] == 1)
  intersection_size <- length(intersect(neighbors_i, neighbors_j))
  union_size <- length(union(neighbors_i, neighbors_j))
  jaccard <- intersection_size / union_size
  return(jaccard)
}

# Function to calculate Preferential Attachment (PA)
preferential_attachment <- function(graph, node_i, node_j) {
  neighbors_i <- which(graph[node_i, ] == 1)
  neighbors_j <- which(graph[node_j, ] == 1)
  attachment <- length(neighbors_i) * length(neighbors_j)
  return(attachment)
}

# Function to calculate Adamic-Adar coefficient (AA)
academic_adar <- function(graph, node_i, node_j) {
  neighbors_i <- which(graph[node_i, ] == 1)
  neighbors_j <- which(graph[node_j, ] == 1)
  common_neighbors <- intersect(neighbors_i, neighbors_j)
  adar_sum <- sum(1 / log(length(which(graph[common_neighbors, ] == 1))))
  return(adar_sum)
}

############################################################################
# COMPUTE GLOBAL SIMILARITY SCORES #########################################
############################################################################

# Function to calculate the Katz score (K)
katz_coefficient <- function(adjacency_matrix, alpha) {
  # Check alpha
  maxEigenvalue <- (1 / abs(eigen(adjacency_matrix)$values[1]))
  if (alpha <= 0 || alpha >= maxEigenvalue) {
    stop("Invalid alpha value.", call. = FALSE)
  }
  
  I <- diag(nrow(adjacency_matrix))
  K_a <- ginv(I - alpha * adjacency_matrix) - I
  diag(K_a) <- 0
  
  # Normalize
  katz_score <- (K_a - min(K_a)) / (max(K_a) - min(K_a))

  return(katz_score)
}


# Function to calculate the Personalized PageRank score (PPR)
personalized_pagerank_formula <- function(adjacency_matrix, alpha) {
  num_nodes <- nrow(adjacency_matrix)
  
  # Calculate the row-stochastic matrix P = D^(-1) * A
  D <- diag(colSums(adjacency_matrix))
  D_inv <- ginv(D)
  P <- D_inv %*% adjacency_matrix
  
  # Identity matrix
  I <- diag(num_nodes)
  
  # Calculate the Personalized PageRank scores using the formula PPR_alpha = (I - alpha * P)^-1
  pagerank_scores <- ginv(I - alpha * P)
  
  return(pagerank_scores)
}


true_edges <- function(graph){
  edges <- which(graph != 0 & row(graph) > col(graph), arr.ind = TRUE)
  return(edges)
}

false_edges <- function(graph){
  edges <- which(graph == 0 & row(graph) > col(graph), arr.ind = TRUE)
}


# Function to get likelihood scores using different approaches
get_likelihood_scores <- function(graph, approach, alpha = NULL) {
  num_nodes <- nrow(graph)
  likelihood_scores <- matrix(0, nrow = num_nodes, ncol = num_nodes)
  
  for (i in 1:(num_nodes - 1)) {
    for (j in (i + 1):num_nodes) {
      if (graph[i, j] == 0) {  # For non-edges only
        score <- switch(
          approach,
          "CN" = common_neighbors(graph, i, j),
          "JC" = jaccard_coefficient(graph, i, j),
          "PA" = preferential_attachment(graph, i, j),
          "AA" = academic_adar(graph, i, j),
          "Random" = runif(1),
          "K" = katz_coefficient(graph, alpha),
          "PPR" = {
            if (is.null(alpha)) stop("Alpha parameter missing for PPR.")
            personalized_pagerank_formula(graph, alpha)[i, j]
          }
        )
        likelihood_scores[i, j] <- score
        likelihood_scores[j, i] <- score  # Graph is undirected, so we set both entries
      }
    }
  }
  
  return(likelihood_scores)
}

compute_betweenness <- function(adj){
  graph <- graph.adjacency(adj, mode = "undirected", weighted = NULL)
  betweenness_centrality <- betweenness(graph)
  avg_betweenness <- mean(betweenness_centrality)
  hist(betweenness_centrality, main = "", xlab = "Betweenness Centrality", col = "skyblue", border = "black")
  abline(v = avg_betweenness, col = "red", lwd = 2, lty = 2)
  legend("topright", legend = c("Betweenness Centrality", paste("Average: ", round(avg_betweenness, 2))), col = c("skyblue", "red"), lty = c(1, 2), lwd = c(1, 2))
}


###########################################################
# FUNCTIONS TO COMPUTE TPR AND TNR ########################
###########################################################

# Function to calculate TPR and TNR for Katz scores
calculate_katz_metrics <- function(adjacency_matrix, alpha_values, true_edges, false_edges, threshold_value) {
  tpr_values <- c()
  tnr_values <- c()
  fpr_values <- c()
  fnr_values <- c()
  
  for (alpha in alpha_values) {
    katz_scores <- katz_coefficient(adjacency_matrix, alpha)
    predicted_edges <- which(katz_scores > threshold_value & row(katz_scores) > col(katz_scores), arr.ind = TRUE)
    
    # Convert true edge list and predicted edge list to dataframes
    df_true <- data.frame(start = true_edges[, 1], end = true_edges[, 2])
    df_pred <- data.frame(start = predicted_edges[, 1], end = predicted_edges[, 2])
    
    common_edges <- as.matrix(merge(df_true, df_pred, by = c("start", "end")))
    
    true_positive <- nrow(common_edges)
    false_negative <- nrow(df_true) - true_positive
    false_positive <- nrow(df_pred) - true_positive
    true_negative <- nrow(true_edges) - true_positive - false_negative - false_positive
    
    # Check for division by zero or NaN/Inf
    tpr_denominator <- true_positive + false_negative
    tnr_denominator <- true_negative + false_positive
    
    tpr <- ifelse(tpr_denominator > 0, true_positive / tpr_denominator, 0)
    tnr <- ifelse(tnr_denominator > 0, true_negative / tnr_denominator, 0)
    
    fpr <- ifelse(tnr_denominator > 0, false_positive / tnr_denominator, NaN)
    fnr <- ifelse(tpr_denominator > 0, false_negative / tpr_denominator, NaN)
    
    # Replace Inf and -Inf with NaN
    fpr[is.infinite(fpr)] <- NaN
    fnr[is.infinite(fnr)] <- NaN
    
    tpr_values <- c(tpr_values, tpr)
    tnr_values <- c(tnr_values, tnr)
    fpr_values <- c(fpr_values, fpr)
    fnr_values <- c(fnr_values, fnr)
  }
  
  return(data.frame(N = length(predicted_edges), Alpha = alpha_values, TPR = tpr_values, TNR = tnr_values, FPR = fpr_values, FNR = fnr_values))
}

# Function to calculate TPR and TNR for Personalized PageRank scores
calculate_ppr_metrics <- function(adjacency_matrix, alpha_values, true_edges, false_edges, threshold_value) {
  tpr_values <- c()
  tnr_values <- c()
  
  for (alpha in alpha_values) {
    pagerank_scores <- personalized_pagerank_formula(adjacency_matrix, alpha)
    predicted_edges <- sort(which(pagerank_scores > threshold_value), decreasing = TRUE)
    
    true_positive <- length(intersect(predicted_edges, true_edges))
    false_negative <- length(setdiff(true_edges, predicted_edges))
    false_positive <- length(setdiff(predicted_edges, true_edges))
    true_negative <- length(intersect(false_edges, predicted_edges))
    
    tpr <- true_positive / (true_positive + false_negative)
    tnr <- true_negative / (true_negative + false_positive)
    
    tpr_values <- c(tpr_values, tpr)
    tnr_values <- c(tnr_values, tnr)
  }
  
  return(data.frame(Alpha = alpha_values, TPR = tpr_values, TNR = tnr_values))
}

# Example usage
alpha_values <- seq(0.001, (1 / eigen(MN)$values[1]), 0.001)

# Identify true edges from the adjacency matrix
true_edges_MN <- true_edges(MN)
true_edges_PC <- true_edges(PC)

# Identify false edges from the adjacency matrix
false_edges_MN <- false_edges(MN)
false_edges_PC <- false_edges(PC)

# Set a threshold value for classification (adjust as needed)
threshold_value <- 0.5

# Calculate TPR and TNR for Katz scores
katz_metrics_MN <- calculate_katz_metrics(MN, alpha_values, true_edges_MN, false_edges_MN, threshold_value)
katz_metrics_PC <- calculate_katz_metrics(PC, alpha_values, true_edges_PC, false_edges_PC, threshold_value)

# Calculate TPR and TNR for Personalized PageRank scores
#ppr_metrics_MN <- calculate_ppr_metrics(MN, alpha_values, true_edges_MN, false_edges_MN, threshold_value)
#ppr_metrics_PC <- calculate_ppr_metrics(PC, alpha_values, true_edges_PC, false_edges_PC, threshold_value)

# Combine metrics for the different networks
combined_katz <- rbind(cbind(Method = "MN", katz_metrics_MN),
                       cbind(Method = "PC", katz_metrics_PC))

# Plotting TPR
ggplot(combined_katz, aes(x = Alpha, y = TPR, color = Method, linetype = Method)) +
  geom_line(size = 1) +
  labs(title = "TPR vs. Alpha for Katz Score",
       x = "Alpha",
       y = "Rate") +
  theme_minimal() +
  scale_linetype_manual(values = c("MN" = "solid", "PC" = "dashed"))

# Plotting TNR
ggplot(combined_katz, aes(x = Alpha, y = TNR, color = Method, linetype = Method)) +
  geom_line(size = 1) +
  labs(title = "TNR vs. Alpha for Katz Score",
       x = "Alpha",
       y = "Rate") +
  theme_minimal() +
  scale_linetype_manual(values = c("MN" = "solid", "PC" = "dashed"))

###########################################################
# EXAMPLE USAGES ##########################################
###########################################################

# Example usage
adjacency_matrix <- WR
approach <- "AA"  # Change to "JC" or "Random" as needed
scores <- get_likelihood_scores(adjacency_matrix, approach)
print(scores)


# Step 1
link_likelihood_scores <- get_likelihood_scores(MN)  # Replace with your function to calculate likelihood scores
sorted_non_edges <- sort_non_edges(link_likelihood_scores)
top_50_percent <- sorted_non_edges[1:round(0.5 * length(sorted_non_edges))]
C <- top_50_percent

# Step 2
p_values <- c(0.01, 0.05, 0.1, 0.15)  # Set your desired p values
for (p in p_values) {
  R_p <- sample(C, size = round(p * length(C)))
  
  # Step 3
  MN0_p <- add_edges_to_graph(MN, R_p)  # Replace with your function to add edges to graph
  # Similarly, update PC0_p if needed
  
  # Step 4
  r <- calculate_relative_variation(MN, MN0_p)  # Replace with your function to calculate relative variation
  # Similarly, update for PC if needed
  
  # Display or store the results as needed
  cat(sprintf("For p = %f, relative variation r = %f\n", p, r))
}

roc_obj <- roc(katz_metrics_MN$FPR, katz_metrics_MN$TPR)
