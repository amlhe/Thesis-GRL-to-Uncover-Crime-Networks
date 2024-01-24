
# Function to put cluster labels into binary form
vec2mat <- function(clust_lab){
  V <- length(clust_lab)
  H <- max(clust_lab)
  M <- matrix(0,V,H)
  for (v in 1:V){
    M[v,clust_lab[v]] <- 1
  }
  return(M)
}

# Stochastic Block Model for single layer binary undirected network
sbm <- function(Y, N_iter, seed, z_init=c(1:nrow(Y)), a = 1, b = 1, alpha = 8){
  
  # Initialization
  set.seed(seed)
  V <- nrow(Y)
  Z <- vec2mat(z_init) # binary VxH matrix with 1 if node v is in cluster h where input z_init is a vector of length V
  z_post <- matrix(NA,V,N_iter)
  
  # Matrix with block connections
  temp   <- Y%*%Z
  m_full <- t(Z)%*%temp - diag(0.5*colSums(temp*Z),ncol(Z))
  
  # Gibbs sampler
  for (t in 1:N_iter){ # iterate over number of Gibbs sweeps
    for (v in 1:V){ # iterate over rows in adjacency matrix 
      
      # Discard empty and singularly populated clusters
      if(ncol(Z) > 1){
        nonempty_v <- which(colSums(Z[-v,]) > 0) # find columns not including v that are populated
        Z <- Z[, nonempty_v] # update Z with populated clusters
        if (length(nonempty_v) == 1){
          Z <- matrix(Z,V,1) 
        }
        
        # Reduce the dimensions of the m_full matrix
        m_full <- matrix(m_full[nonempty_v,nonempty_v],ncol(Z),ncol(Z))
      }
      
      H <- ncol(Z) # number of active clusters
      Z_v <- Z[-v,] # Z matrix excluding node v
      
      if (H == 1){
        v_minus <- sum(Z[-v]) # number of nodes in each cluster excluding node v
      } else {
        v_minus <- colSums(Z_v)
      }
      
      r_v <- crossprod(Z_v, Y[-v,v]) # number of edges between node v and each cluster
      h_v <- which(Z[v,] > 0) # indices of the clusters to which vertex v is assigned
      
      # Compute the m matrix by difference
      if(length(h_v) == 1){
        resid1 <- matrix(0,H,H)
        resid1[,h_v] <- r_v # values in columns set to r_v
        resid1[h_v,] <- r_v # values in rows set to r_v
        m <- m_full - resid1 # update expected connection between clusters
      } else {
        m <- m_full
      } 
      
      m_bar <- matrix(v_minus,H,1) %*% matrix(v_minus,1,H) - diag(0.5*v_minus*(v_minus+1),H) - m # number of non-edges between cluster h and cluster k, excluding node v 
      V_minus <- matrix(1,H,1) %*% matrix(v_minus,1,H)
      R <- matrix(1,H,1) %*% matrix(r_v,1,H)
      
      # Compute probabilities
      log_lhds_old <- rowSums(lbeta(m + R + a, m_bar + V_minus - R + b) - lbeta(m + a, m_bar + b))
      log_lhd_new  <- sum(lbeta(r_v + a, v_minus - r_v + b) - lbeta(a, b))
      
      log_p <- log(c(v_minus,alpha)) + c(log_lhds_old, log_lhd_new)
      max_log_p <- max(log_p)
      p <- exp(log_p - max_log_p)
      p <- p / sum(p)
      
      new_sample <- which(rmultinom(1,1,p) > 0)
      
      if(length(h_v) == 1){
        Z[v,h_v] <- 0
      }
      
      if(new_sample == H+1){
        Z <- cbind(Z,rep(0,V))
        Z[v,new_sample] <- 1
        m <- rbind(cbind(m,0),0)
        H <- H + 1
        r_v <- crossprod(Z[-v,], Y[-v,v])
      } else {
        Z[v, new_sample] <- 1
      }

      resid2  <- matrix(0,H,H)
      resid2[,new_sample] <- r_v
      resid2[new_sample,] <- r_v
      m_full <- m + resid2
    }
    
    z_post[,t] <- Z %*% c(1:ncol(Z))
  }
  
  return(z_post[,N_iter])
}


# Stochastic Block Model for multi-layer binary undirected network with layer-specific link probabilities
multilayer_sbm_layer_specific <- function(Y_list, N_iter, seed, z_init = c(1:nrow(Y_list[,,1])), a = 1, b = 1, alpha = 8) {
  
  # Initialization
  set.seed(seed)
  V <- nrow(Y_list[,,1])
  H <- length(z_init)
  L <- length(Y_list[1,1,])
  Z <- vec2mat(z_init)
  z_post <- matrix(NA, V, N_iter)
  log_p_vec <- list()
  
  # Gibbs sampler
  for (t in 1:N_iter) {  # iterate over number of Gibbs sweeps
    for (v in 1:V) {  # iterate over rows in adjacency matrix
      
      # Iterate over layers
      for (l in 1:L) {
        Y <- Y_list[,,l]
        
        # Update model parameters
        temp <- Y %*% Z
        m_full <- t(Z) %*% temp - diag(0.5 * colSums(temp * Z), ncol(Z))
        
        # Discard empty and singularly populated clusters
        if (ncol(Z) > 1) {
          nonempty_v <- which(colSums(Z[-v, ]) > 0)  # find columns not including v that are populated
          Z <- Z[, nonempty_v]  # update Z with populated clusters
          if (length(nonempty_v) == 1) {
            Z <- matrix(Z, V, 1)
          }
          # Reduce the dimensions of the m_full matrix
          m_full <- matrix(m_full[nonempty_v,nonempty_v],ncol(Z),ncol(Z))
        }
        
        H <- ncol(Z)  # number of active clusters
        Z_v <- Z[-v, ]  # Z matrix excluding node v
        
        if (H == 1) {
          v_minus <- sum(Z[-v])  # number of nodes in each cluster excluding node v
        } else {
          v_minus <- colSums(Z_v)
        }
        
        r_v <- crossprod(Z_v, Y[-v, v])  # number of edges between node v and each cluster in each layer
        h_v <- which(Z[v, ] > 0)  # indices of the clusters to which vertex v is assigned
        
        # Compute the m matrix by difference
        if (length(h_v) == 1) {
          resid1 <- matrix(0, H, H)
          resid1[, h_v] <- r_v # values in columns set to r_v
          resid1[h_v, ] <- r_v # values in rows set to r_v
          m <- m_full - resid1 # update expected connection between clusters
        } else {
          m <- m_full
        }
        
        m_bar <- matrix(v_minus, H, 1) %*% matrix(v_minus, 1, H) - diag(0.5 * v_minus * (v_minus + 1), H) - m  # number of non-edges between cluster h and cluster k, excluding node v
        V_minus <- matrix(1, H, 1) %*% matrix(v_minus, 1, H)
        R <- matrix(1, H, 1) %*% matrix(r_v, 1, H)
        
        # Compute probabilities
        log_lhds_old <- apply(lbeta(m + R + a, m_bar + V_minus - R + b) - lbeta(m + a, m_bar + b), 1, function(row) sum(ifelse(is.finite(row), row, 0)))
        log_lhd_new <- sum(lbeta(r_v + a, v_minus - r_v + b) - lbeta(a, b))
        log_p <- log(c(v_minus, alpha)) + c(log_lhds_old, log_lhd_new)
        
        # Normalize probabilities within each layer
        max_log_p <- max(log_p)
        p <- exp(log_p - max_log_p)
        p <- p / sum(p)
        
        # Accumulate probabilities across layers
        log_p_vec[[l]] <- p
      }
      
      log_p_mat <- do.call(cbind, log_p_vec)
      
      # Likelihood as product
      log_p_layer <- apply(log_p_mat, 1, prod)
      
      # Sample new assignment
      new_sample <- which(rmultinom(1, 1, log_p_layer) > 0)
      
      if(length(h_v) == 1) {
        Z[v, h_v] <- 0
      }
      
      if (new_sample == H + 1) {
        Z <- cbind(Z, rep(0, V))
        Z[v, new_sample] <- 1
        m <- rbind(cbind(m, 0), 0)
        H <- H + 1
        r_v <- crossprod(Z[-v, ], Y[-v, v])  # Update r_v for the new layer
      } else {
        Z[v, new_sample] <- 1
      }
      
      resid2 <- matrix(0, H, H)
      resid2[, new_sample] <- r_v
      resid2[new_sample, ] <- r_v
      m_full <- m + resid2
    }
    
    z_post[, t] <- Z %*% c(1:ncol(Z))
  }
  
  return(z_post)
}

# Multilayer Stochastic Block Model with shared community structure and link probabilities
multilayer_sbm_shared <- function(Y_list, N_iter, seed, z_init = c(1:nrow(Y_list[,,1])), a = 1, b = 1, alpha = 8) {
  
  # Initialization
  set.seed(seed)
  V <- nrow(Y_list[,,1])
  H <- length(z_init)
  L <- length(Y_list[1,1,])
  Z <- vec2mat(z_init)
  z_post <- matrix(NA, V, N_iter)
  
  # Gibbs sampler
  for (t in 1:N_iter) {  # iterate over the number of Gibbs sweeps
    for (v in 1:V) {  # iterate over rows in the adjacency matrix
      
      # Accumulate statistics for each layer
      m_sum <- R_sum <- a_sum <- m_bar_sum <- V_minus_sum <- b_sum <- v_minus_sum <- r_v_sum <- 0
      
      # Iterate over layers
      for (l in 1:L) {
        Y <- Y_list[,,l]
        
        # Update model parameters
        temp <- Y %*% Z
        m_full <- t(Z) %*% temp - diag(0.5 * colSums(temp * Z), ncol(Z))
        
        # Discard empty and singularly populated clusters
        if (ncol(Z) > 1) {
          nonempty_v <- which(colSums(Z[-v, ]) > 0)  # find columns not including v that are populated
          Z <- Z[, nonempty_v]  # update Z with populated clusters
          if (length(nonempty_v) == 1) {
            Z <- matrix(Z, V, 1)
          }
          # Reduce the dimensions of the m_full matrix
          m_full <- matrix(m_full[nonempty_v,nonempty_v],ncol(Z),ncol(Z))
        }
        
        H <- ncol(Z)  # number of active clusters
        Z_v <- Z[-v, ]  # Z matrix excluding node v
        
        if (H == 1) {
          v_minus <- sum(Z[-v])  # number of nodes in each cluster excluding node v
        } else {
          v_minus <- colSums(Z_v)
        }
        
        r_v <- crossprod(Z_v, Y[-v, v])  # number of edges between node v and each cluster in each layer
        h_v <- which(Z[v, ] > 0)  # indices of the clusters to which vertex v is assigned
        
        # Compute the m matrix by difference
        if (length(h_v) == 1) {
          resid1 <- matrix(0, H, H)
          resid1[, h_v] <- r_v # values in columns set to r_v
          resid1[h_v, ] <- r_v # values in rows set to r_v
          m <- m_full - resid1 # update expected connection between clusters
        } else {
          m <- m_full
        }
        
        m_bar <- matrix(v_minus, H, 1) %*% matrix(v_minus, 1, H) - diag(0.5 * v_minus * (v_minus + 1), H) - m  # number of non-edges between cluster h and cluster k, excluding node v
        V_minus <- matrix(1, H, 1) %*% matrix(v_minus, 1, H)
        R <- matrix(1, H, 1) %*% matrix(r_v, 1, H)
        
        # Accumulate statistics
        m_sum <- m_sum + sum(m)
        R_sum <- R_sum + sum(R)
        a_sum <- a_sum + a
        m_bar_sum <- m_bar_sum + sum(m_bar)
        V_minus_sum <- V_minus_sum + sum(V_minus)
        b_sum <- b_sum + b
        v_minus_sum <- v_minus_sum + sum(v_minus)
        r_v_sum <- r_v_sum + sum(r_v)
      }
      
      # Compute probabilities
      log_lhds_old <- apply(lbeta(m + R + a, m_bar + V_minus - R + b) - lbeta(m + a, m_bar + b), 1, function(row) sum(ifelse(is.finite(row), row, 0)))
      log_lhd_new <- sum(lbeta(r_v + a, v_minus - r_v + b) - lbeta(a, b))
      log_p_layer <- log(c(v_minus, alpha)) + c(log_lhds_old, log_lhd_new)
      
      max_log_p_layer <- max(log_p_layer)
      p_layer <- exp(log_p_layer - max_log_p_layer)
      p_layer <- p_layer / sum(p_layer)
      
      # Sample new assignment
      new_sample <- which(rmultinom(1, 1, p_layer) > 0)
      
      if (length(h_v) == 1) {
        Z[v, h_v] <- 0
      }
      
      if (new_sample == H + 1) {
        Z <- cbind(Z, rep(0, V))
        Z[v, new_sample] <- 1
        m <- rbind(cbind(m, 0), 0)
        H <- H + 1
        r_v <- crossprod(Z[-v, ], Y[-v, v])  # Update r_v for the new layer
      } else {
        Z[v, new_sample] <- 1
      }
      
      resid2 <- matrix(0, H, H)
      resid2[, new_sample] <- r_v
      resid2[new_sample, ] <- r_v
      m_full <- m + resid2
    }
    
    z_post[, t] <- Z %*% c(1:ncol(Z))
  }
  
  return(z_post)
}
