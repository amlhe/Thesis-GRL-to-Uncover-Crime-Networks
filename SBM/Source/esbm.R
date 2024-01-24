####################################################################################
# FUNCTIONS TO COMPUTE THE DIFFERENT URN SCHEMES FOR THE PRIORS IN THE ARTICLE #####
####################################################################################

urn_DP <- function(v_minus, alpha_PY) {
  return(c(v_minus,alpha_PY))
}

urn_PY <- function(v_minus,alpha_PY,sigma_PY){
  H<-length(v_minus)
  return(c(v_minus-sigma_PY,alpha_PY+H*sigma_PY))
}

urn_DM <- function(v_minus,beta_DM,H_DM){
  H<-length(v_minus)
  return(c(v_minus+beta_DM,beta_DM*(H_DM-H)*(H_DM>H)))
}

urn_GN <- function(v_minus,gamma_GN){
  H<-length(v_minus)
  return(c((v_minus+1)*(sum(v_minus)-H+gamma_GN),H^2-H*gamma_GN))
}

####################################################################################
# GIBBS SAMPLER FOR THE EXTENDED STOCHASTIC BLOCK MODEL  ###########################
####################################################################################

# Inputs:
# Y = VxV symmetric adjacency matrix
# z_init = V-vector of initialization assignment for each node (default = one cluster for each node)
# a,b = parameters of the Beta prior on the thetas
# prior = string in c("DP","PY","DM","GN")
# if prior=="DP", alpha_PY must be set
# if prior=="PY", alpha_PY and sigma_PY must be set
# if prior=="DM", beta_DM and H_DM must be set
# if prior=="GN", gamma_GN must be set
# x = V-vector of categorical covariates
# if x is specified, also alpha_xi (a C-vector of parameters for a Dirichlet distribution) must be set

# Output:
# Posterior samples of the community labels for each node v=1,...,V

esbm <- function(Y, seed, N_iter, prior, z_init=c(1:nrow(Y)), a=1, b=1,
                 alpha_PY=NA, sigma_PY=NA, beta_DM=NA, H_DM=NA, gamma_GN=NA, x=NULL, alpha_xi=NULL){
  
  # ----------------------------------------------
  # Selection of the prior distribution to be used
  # ----------------------------------------------
  
  if (prior=="DP"){
    urn<-function(v_minus){return(urn_DP(v_minus,alpha_PY))}
  } else if (prior=="PY"){
    urn<-function(v_minus){return(urn_PY(v_minus,alpha_PY,sigma_PY))}
  } else if (prior=="DM"){
    urn<-function(v_minus){return(urn_DM(v_minus,beta_DM,H_DM))}
  } else if (prior=="GN"){
    urn<-function(v_minus){return(urn_GN(v_minus,gamma_GN))}
  } else { 
    stop("Invalid value for prior")  
  }
  
  # ----------------------------------------------
  # Pre-processing of the node attributes
  # ----------------------------------------------
  
  if (!is.null(x)){
    print("Covariates have been provided")
    x <- as.numeric(as.factor(x))
    X <- vec2mat(x)
    if (!is.null(alpha_xi)){
      alpha0 <- sum(alpha_xi)
    } else {
      stop("If covariates x are given, then alpha_xi must be set as well")
    }
  }
  
  # ----------------------------------------------
  # Initialization
  # ----------------------------------------------
  
  set.seed(seed)
  
  V <- nrow(Y)
  
  # cluster assignments are encoded in two equivalent ways:
  # [i] a VxH matrix Z, s.t. Z[v,h]=1{node v in cluster h}, faster to use within each iteration
  Z <- vec2mat(z_init)
  
  # [ii] a vector of length V containing the cluster label for each node, more compact to store;
  # such vectors for all iterations are packed in a VxN_iter matrix z_post, 
  # s.t. z_post[v,t]=h if node v is in cluster h at iteration t
  # Note: matrix z_post takes less memory than a list of N_iter matrices Z
  z_post <- matrix(NA,V,N_iter)
  
  
  # Create the matrix with block connections
  temp   <- Y%*%Z
  m_full <- t(Z)%*%temp - diag(0.5*colSums(temp*Z),ncol(Z))
  
  # ----------------------------------------------
  # Beginning of the Gibbs sampler
  # ----------------------------------------------
  
  for (t in 1:N_iter){
    for (v in 1:V){
      
      # remove empty clusters and
      # if the cluster containing node v has no other node, discard it as well
      if(ncol(Z) > 1){
        nonempty_v <- which(colSums(Z[-v,]) > 0)  
        Z <- Z[, nonempty_v]
        if (length(nonempty_v)==1){Z <- matrix(Z,V,1)}
        
        # Reduce the dimensions of the m_full matrix
        m_full <- matrix(m_full[nonempty_v,nonempty_v],ncol(Z),ncol(Z))
      } 
      
      # H = number of active clusters
      H   <- ncol(Z)
      Z_v <- Z[-v,]
      
      # v_minus = number of nodes in each cluster, excluding node v
      if (H==1){v_minus <- sum(Z[-v])} else {v_minus <- colSums(Z_v)}
      
      # r_v = number of edges from node v to each cluster (no self-loops)
      r_v         <- crossprod(Z_v, Y[-v,v])
      h_v         <- which(Z[v,] > 0)
      
      # Compute the m matrix by difference
      if(length(h_v) == 1){
        resid1       <- matrix(0,H,H)
        resid1[,h_v] <- r_v; resid1[h_v,] <- r_v
        m            <- m_full - resid1
      } else {m <- m_full} # No need to update m in this case
      
      # m_bar = number of non-edges between cluster h and cluster k, excluding node v 
      m_bar   <- matrix(v_minus,H,1)%*%matrix(v_minus,1,H) - diag(0.5*v_minus*(v_minus+1),H) - m
      V_minus <- matrix(1,H,1)%*%matrix(v_minus,1,H)
      R       <- matrix(1,H,1)%*%matrix(r_v,1,H)
      
      # ----------------------------------------------
      # Computing the probabilities
      # ----------------------------------------------
      
      #log_lhds_old <- rowSums(lbeta(m + R + a, m_bar + V_minus - R + b) - lbeta(m + a, m_bar + b)) # vector of length H
      # Calculate the log likelihoods with handling NaN or Inf
      log_lhds_old <- apply(lbeta(m + R + a, m_bar + V_minus - R + b) - lbeta(m + a, m_bar + b), 1, function(row) sum(ifelse(is.finite(row), row, 0)))
      log_lhd_new  <- sum(lbeta(r_v + a, v_minus - r_v + b) - lbeta(a, b)) # scalar
      log_addit    <- 0
      
      if(!is.null(x)){
        Vx        <- crossprod(Z_v, X[-v,])
        addit_old <- (Vx[,x[v]] + alpha_xi[x[v]]) / (v_minus+alpha0)
        addit_new <- alpha_xi[x[v]] / alpha0
        log_addit <- log(c(addit_old, addit_new))
      }
      
      # Probabilities
      log_p <- log_addit + log(urn(v_minus)) + c(log_lhds_old, log_lhd_new)
      max_log_p <- max(log_p)
      p <- exp(log_p - max_log_p)
      p <- p / sum(p)
      #p     <- exp(log_p - max(log_p)); #p <- p / sum(p)
      
      # ----------------------------------------------
      # Sampling the indicator
      # ----------------------------------------------
      
      new_sample <- which(rmultinom(1,1,p) > 0)
      
      # ----------------------------------------------
      # Adjusting Z, H, r_v and m
      # ----------------------------------------------
      
      if(length(h_v) == 1){Z[v,h_v] <- 0}
      
      # If a new value is sampled, increase the dimension of m_full
      if(new_sample== H+1){
        Z               <- cbind(Z,rep(0,V))
        Z[v,new_sample] <- 1
        m               <- rbind(cbind(m,0),0)
        H               <- H + 1
        r_v             <- crossprod(Z[-v,], Y[-v,v])
      } else {Z[v, new_sample] <- 1}
      
      # Updating m_full
      resid2              <- matrix(0,H,H)
      resid2[,new_sample] <- r_v; resid2[new_sample,] <- r_v
      m_full              <- m + resid2
    }
    
    # store cluster assignments at time t in matrix z_post s.t.
    # z_post[v,t]=h if node v is in cluster h at iteration t
    z_post[,t] <- Z %*% c(1:ncol(Z))
    
    #print(table(z_post[,t])) 
    if (t%%1000 == 0){print(paste("Iteration:", t))}
  }
  return(z_post)
}


####################################################################################
# EXTENTIONS OF THE ESBM TO ACCOUNT FOR MULTILAYER NETWORKS  #######################
####################################################################################

# Extended Stochastic Block Model for multi-layer undirected network with layer-specific link probabilities
multilayer_esbm_layer_specific <- function(Y_list, seed, N_iter, prior, z_init=c(1:nrow(Y_list[,,1])), a=1, b=1,
                 alpha_PY=NA, sigma_PY=NA, beta_DM=NA, H_DM=NA, gamma_GN=NA, x=NULL, alpha_xi=NULL){
  
  # ----------------------------------------------
  # Selection of the prior distribution to be used
  # ----------------------------------------------
  
  if (prior=="DP"){
    urn<-function(v_minus){return(urn_DP(v_minus,alpha_PY))}
  } else if (prior=="PY"){
    urn<-function(v_minus){return(urn_PY(v_minus,alpha_PY,sigma_PY))}
  } else if (prior=="DM"){
    urn<-function(v_minus){return(urn_DM(v_minus,beta_DM,H_DM))}
  } else if (prior=="GN"){
    urn<-function(v_minus){return(urn_GN(v_minus,gamma_GN))}
  } else { 
    stop("Invalid value for prior")  
  }
  
  # Initialization
  set.seed(seed)
  V <- nrow(Y_list[,,1])
  H <- length(z_init)
  L <- length(Y_list[1,1,])
  Z <- vec2mat(z_init)
  z_post <- matrix(NA, V, N_iter)
  log_p_vec <- list()
  
  # Gibbs sampler
  for (t in 1:N_iter) {  # iterate over the number of Gibbs sweeps
    for (v in 1:V) {  # iterate over rows in the adjacency matrix
      
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
        
        # Computing the probabilities
        log_lhds_old <- apply(lbeta(m + R + a, m_bar + V_minus - R + b) - lbeta(m + a, m_bar + b), 1, function(row) sum(ifelse(is.finite(row), row, 0)))
        log_lhd_new <- sum(lbeta(r_v + a, v_minus - r_v + b) - lbeta(a, b))
        log_p <- log(urn(v_minus)) + c(log_lhds_old, log_lhd_new)
        
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
      
      
      # Adjusting Z, H, r_v and m
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
    
    # store cluster assignments at time t in matrix z_post s.t.
    # z_post[v,t]=h if node v is in cluster h at iteration t
    z_post[,t] <- Z %*% c(1:ncol(Z))
    
    if (t%%1000 == 0){print(paste("Iteration:", t))}
  }
  return(z_post)
}


# Stochastic Block Model for multi-layer binary undirected network with layer-specific link probabilities
multilayer_esbm_shared <- function(Y_list, seed, N_iter, prior, z_init=c(1:nrow(Y_list[,,1])), a=1, b=1,
                 alpha_PY=NA, sigma_PY=NA, beta_DM=NA, H_DM=NA, gamma_GN=NA, x=NULL, alpha_xi=NULL) {
  
  # ----------------------------------------------
  # Selection of the prior distribution to be used
  # ----------------------------------------------
  
  if (prior=="DP"){
    urn<-function(v_minus){return(urn_DP(v_minus,alpha_PY))}
  } else if (prior=="PY"){
    urn<-function(v_minus){return(urn_PY(v_minus,alpha_PY,sigma_PY))}
  } else if (prior=="DM"){
    urn<-function(v_minus){return(urn_DM(v_minus,beta_DM,H_DM))}
  } else if (prior=="GN"){
    urn<-function(v_minus){return(urn_GN(v_minus,gamma_GN))}
  } else { 
    stop("Invalid value for prior")  
  }
  
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
    
        # H = number of active clusters
        H   <- ncol(Z)
        Z_v <- Z[-v,]
        
        # v_minus = number of nodes in each cluster, excluding node v
        if (H==1){v_minus <- sum(Z[-v])} else {v_minus <- colSums(Z_v)}
        
        # r_v = number of edges from node v to each cluster (no self-loops)
        r_v         <- crossprod(Z_v, Y[-v,v])
        h_v         <- which(Z[v,] > 0)
        
        # Compute the m matrix by difference
        if(length(h_v) == 1){
          resid1       <- matrix(0,H,H)
          resid1[,h_v] <- r_v; resid1[h_v,] <- r_v
          m            <- m_full - resid1
        } else {m <- m_full} # No need to update m in this case
        
        # m_bar = number of non-edges between cluster h and cluster k, excluding node v 
        m_bar   <- matrix(v_minus,H,1)%*%matrix(v_minus,1,H) - diag(0.5*v_minus*(v_minus+1),H) - m
        V_minus <- matrix(1,H,1)%*%matrix(v_minus,1,H)
        R       <- matrix(1,H,1)%*%matrix(r_v,1,H)
        
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
      
      # Computing the probabilities
      log_lhds_old <- apply(lbeta(m + R + a, m_bar + V_minus - R + b) - lbeta(m + a, m_bar + b), 1, function(row) sum(ifelse(is.finite(row), row, 0)))
      log_lhd_new  <- sum(lbeta(r_v + a, v_minus - r_v + b) - lbeta(a, b)) 
      log_p_layer <- log(urn(v_minus)) + c(log_lhds_old, log_lhd_new)
      
      max_log_p_layer <- max(log_p_layer)
      p_layer <- exp(log_p_layer - max_log_p_layer)
      p_layer <- p_layer / sum(p_layer)
      
      # Sample new assignment
      new_sample <- which(rmultinom(1, 1, p_layer) > 0)
      
      # ----------------------------------------------
      # Adjusting Z, H, r_v and m
      # ----------------------------------------------
      
      if(length(h_v) == 1){Z[v,h_v] <- 0}
      
      # If a new value is sampled, increase the dimension of m_full
      if(new_sample== H+1){
        Z               <- cbind(Z,rep(0,V))
        Z[v,new_sample] <- 1
        m               <- rbind(cbind(m,0),0)
        H               <- H + 1
        r_v             <- crossprod(Z[-v,], Y[-v,v])
      } else {Z[v, new_sample] <- 1}
      
      # Updating m_full
      resid2              <- matrix(0,H,H)
      resid2[,new_sample] <- r_v; resid2[new_sample,] <- r_v
      m_full              <- m + resid2
    }
    
    # store cluster assignments at time t in matrix z_post s.t.
    # z_post[v,t]=h if node v is in cluster h at iteration t
    z_post[,t] <- Z %*% c(1:ncol(Z))
    
    #print(table(z_post[,t])) 
    if (t%%1000 == 0){print(paste("Iteration:", t))}
  }
  return(z_post)
}


####################################################################################
# PUT CLUSTER LABELS IN BINARY MATRIX FORM  ########################################
####################################################################################

vec2mat <- function(clust_lab){
  # in: vector clust_lab of length V s.t. clust_lab[v]=h if node v is in cluster h
  # out: binary VxH matrix M s.t. M[v,h]=1{node v is in cluster h}
  V <- length(clust_lab)
  H <- max(clust_lab)
  M <- matrix(0,V,H)
  for (v in 1:V){
    M[v,clust_lab[v]] <- 1
  }
  return(M)
}


####################################################################################
# COMPUTE POSTERIOR CO-CLUSTERING MATRIX  ##########################################
####################################################################################

pr_cc <- function(z_post){
  # in: posterior sample of assignments (VxN_iter matrix)
  # out: VxV matrix c with elements c[vu]=fraction of iterations in which v and u are in the same cluster
  V <- nrow(z_post)    
  N_iter <- ncol(z_post)
  c <- matrix(0,V,V)
  for (t in 1:N_iter){
    Z <- vec2mat(z_post[,t])
    c <- c + Z%*%t(Z)
  }
  return(c/N_iter)
}

####################################################################################
# COMPUTE MISCLASSIFICATION ERROR  #################################################
####################################################################################

misclass <- function(memb,Y,a,b){
 # in: vector of cluster labels (memb), VxV adjancency matrix (Y) and hyperparameters beta priors (a,b)
 # out: misclassification error when predicting the edges with the estimated block probabilities
z <- dummy(memb)
H <- ncol(z)
V <- dim(Y)[1]
Abs_Freq <- t(z)%*%Y%*%z
diag(Abs_Freq) <- diag(Abs_Freq)/2
Tot <- t(z)%*%matrix(1,V,V)%*%z
diag(Tot) <- (diag(Tot)-table(memb))/2
Rel_Freq <- (a+Abs_Freq)/(a+b+Tot)
pred <- z%*%Rel_Freq%*%t(z)
return(1-sum(diag(table(lowerTriangle(Y),lowerTriangle(pred>0.5))))/length(lowerTriangle(Y)))
}


####################################################################################
# COMPUTE AUC-ROC AND ROC CURVES ##################################################
####################################################################################

misclass_roc <- function(memb, Y, a=1, b=1) {
  z <- dummy(memb)
  H <- ncol(z)
  V <- dim(Y)[1]
  Abs_Freq <- t(z) %*% Y %*% z
  diag(Abs_Freq) <- diag(Abs_Freq) / 2
  Tot <- t(z) %*% matrix(1, V, V) %*% z
  diag(Tot) <- (diag(Tot) - table(memb)) / 2
  Rel_Freq <- (a + Abs_Freq) / (a + b + Tot)
  pred <- z %*% Rel_Freq %*% t(z)
  
  true_labels <- as.numeric(lowerTriangle(Y))
  predicted_scores <- as.numeric(lowerTriangle(pred))
  
  # Create a prediction object
  prediction <- prediction(predicted_scores, true_labels)
  
  # Create a ROC curve
  roc_curve <- performance(prediction, "tpr", "fpr")
  
  # Create a PR curve
  pr_curve <- performance(prediction, "prec", "rec")
  
  # Calculate AUC score
  auc_roc <- performance(prediction, "auc")@y.values[[1]]
  auc_pr <- performance(prediction, "aucpr")@y.values[[1]]
  acc <- performance(prediction, "acc")@y.values[[1]]
  
  return(list(perf = roc_curve, auc_roc = auc_roc, pr_curve = pr_curve, auc_pr = auc_pr, acc = acc))
}


####################################################################################
# COMPUTE MATRIX OF ESTIMATED EDGE PROBABILITIES  ##################################
####################################################################################

edge_est <- function(memb,Y,a,b){
 # in: vector of cluster labels (memb), VxV adjancency matrix (Y) and hyperparameters beta priors (a,b)
 # out: matrix of estimated edge probabilities
z <- dummy(memb)
H <- ncol(z)
V <- dim(Y)[1]
Abs_Freq <- t(z)%*%Y%*%z
diag(Abs_Freq) <- diag(Abs_Freq)/2
Tot <- t(z)%*%matrix(1,V,V)%*%z
diag(Tot) <- (diag(Tot)-table(memb))/2
Rel_Freq <- (a+Abs_Freq)/(a+b+Tot)
edge_matr <- z%*%Rel_Freq%*%t(z)
diag(edge_matr)<-0
return(edge_matr)
}

####################################################################################
# COMPUTE LOG MARGINAL LIKELIHOOD  #################################################
####################################################################################

log_pY_z <- function(Y,z,a,b){
# in: Adjacency matrix Y, one vector of node labels z, hyperparameters (a,b) of Beta priors for block probabilities
# out: logarithm of the marginal likelihood in eq. [1] evaluated at z.

H <- length(unique(z))
colnames(Y) <- rownames(Y) <- z

edge_counts <- melt(Y)

Y_c <- 1 - Y
diag(Y_c) <- 0
non_edge_counts <- melt(Y_c)
	
Edge <- matrix(aggregate(edge_counts[,3],by=list(edge_counts[,1],edge_counts[,2]),sum,na.rm=TRUE)[,3],H,H)
diag(Edge) <- diag(Edge)/2

No_Edge <- matrix(aggregate(non_edge_counts[,3],by=list(non_edge_counts[,1],non_edge_counts[,2]),sum,na.rm=TRUE)[,3],H,H)
diag(No_Edge) <- diag(No_Edge)/2

a_n <- lowerTriangle(Edge,diag=TRUE)+a
b_bar_n <- lowerTriangle(No_Edge,diag=TRUE)+b

return(sum(lbeta(a_n,b_bar_n))-(H*(H+1)/2)*lbeta(a,b))
}

####################################################################################
# COMPUTE THE LOG-LIKELIHOOD OF THE EDGES ##########################################
####################################################################################

sampleLL <- function(memb,Y,a,b){
  # in: vector of cluster labels (memb), VxV adjancency matrix (Y) and hyperparameters beta priors (a,b)
  # out: vector of Bernoulli log-likelihoods for the edges under ESBM (conditioned on memb and on block-probabilities)
  
  z <- dummy(memb)
  H <- ncol(z)
  V <- dim(Y)[1]
  
  M <- t(z)%*%Y%*%z # block-wise edge counts
  diag(M) <- diag(M)/2
  Tot <- t(z)%*%matrix(1,V,V)%*%z # number of edges for each community
  diag(Tot) <- (diag(Tot)-table(memb))/2
  Mbar <- Tot-M
  a_n <- lowerTriangle(M,diag=TRUE)+a # shape parameter
  b_bar_n <- lowerTriangle(Mbar,diag=TRUE)+b # scale parameter 
  
  theta <- rbeta(length(a_n),a_n,b_bar_n) # samples from beta distribution
  Theta <- matrix(0,H,H) 
  Theta[lower.tri(Theta,diag=TRUE)] <- theta
  Theta <- Theta+t(Theta)
  diag(Theta) <- diag(Theta)/2
  edge_prob <- z%*%Theta%*%t(z) # edge probability matrix 
  
  LL <- dbinom(lowerTriangle(Y,diag=FALSE), size=1, prob=lowerTriangle(edge_prob,diag=FALSE),log=TRUE) # log-likelihood using binomial
  
  return(LL)
}

####################################################################################
# COMPUTE OUT-OF-SAMPLE CLUSTER PROBABILITIES  #####################################
####################################################################################

pred_esbm <- function(Y, prior, z_hat, a=1, b=1,alpha_PY=NA, sigma_PY=NA, beta_DM=NA, H_DM=NA, gamma_GN=NA){

# in: Adjacency matrix Y whose last row and column contain the observed edges for the new node to be classified, one vector of node labels z_hat for the already observed nodes, hyperparameters (a,b) of Beta priors for block probabilities, and prior for the partition process
# out: full conditional probabilities for the membership of the new node to the different clusters.

	
  # ----------------------------------------------
  # Selection of the prior distribution to be used
  # ----------------------------------------------
  
  if (prior=="DP"){
    urn<-function(v_minus){return(urn_DP(v_minus,alpha_PY))}
  } else if (prior=="PY"){
    urn<-function(v_minus){return(urn_PY(v_minus,alpha_PY,sigma_PY))}
  } else if (prior=="DM"){
    urn<-function(v_minus){return(urn_DM(v_minus,beta_DM,H_DM))}
  } else if (prior=="GN"){
    urn<-function(v_minus){return(urn_GN(v_minus,gamma_GN))}
  } else { 
    stop("Invalid value for prior")  
  }
  
  # ----------------------------------------------
  # Initialization
  # ----------------------------------------------
    
  V <- nrow(Y)
  z_init <- c(z_hat,max(z_hat)+1)
   
  # cluster assignments are encoded in two equivalent ways:
  # [i] a VxH matrix Z, s.t. Z[v,h]=1{node v in cluster h}, faster to use within each iteration
  Z <- vec2mat(z_init)  
  
  # Create the matrix with block connections
  temp   <- Y%*%Z
  m_full <- t(Z)%*%temp - diag(0.5*colSums(temp*Z),ncol(Z))
  
  # ----------------------------------------------
  # Predictive
  # ----------------------------------------------
    
  v <- V
  
   # remove empty clusters and
      # if the cluster containing node v has no other node, discard it as well
      if(ncol(Z) > 1){
        nonempty_v <- which(colSums(Z[-v,]) > 0)  
        Z <- Z[, nonempty_v]
        if (length(nonempty_v)==1){Z <- matrix(Z,V,1)}
        
        # Reduce the dimensions of the m_full matrix
        m_full <- matrix(m_full[nonempty_v,nonempty_v],ncol(Z),ncol(Z))
      } 
      
      # H = number of active clusters
      H   <- ncol(Z)
      Z_v <- Z[-v,]
      
      # v_minus = number of nodes in each cluster, excluding node v
      if (H==1){v_minus <- sum(Z[-v])} else {v_minus <- colSums(Z_v)}
      
      # r_v = number of edges from node v to each cluster (no self-loops)
      r_v         <- crossprod(Z_v, Y[-v,v])
      h_v         <- which(Z[v,] > 0)
      
      # Compute the m matrix by difference
      if(length(h_v) == 1){
        resid1       <- matrix(0,H,H)
        resid1[,h_v] <- r_v; resid1[h_v,] <- r_v
        m            <- m_full - resid1
      } else {m <- m_full} # No need to update m in this case
      
      # m_bar = number of non-edges between cluster h and cluster k, excluding node v 
      m_bar   <- matrix(v_minus,H,1)%*%matrix(v_minus,1,H) - diag(0.5*v_minus*(v_minus+1),H) - m
      V_minus <- matrix(1,H,1)%*%matrix(v_minus,1,H)
      R       <- matrix(1,H,1)%*%matrix(r_v,1,H)
      
      # ----------------------------------------------
      # Computing the probabilities
      # ----------------------------------------------
      
      log_lhds_old <- rowSums(lbeta(m + R + a, m_bar + V_minus - R + b) - lbeta(m + a, m_bar + b)) # vector of length H
      log_lhd_new  <- sum(lbeta(r_v + a, v_minus - r_v + b) - lbeta(a, b)) # scalar
      log_addit    <- 0  
  
      log_p <- log_addit + log(urn(v_minus)) + c(log_lhds_old, log_lhd_new)
      p     <- exp(log_p - max(log_p)); 
            
      return(p)
}



####################################################################################
# COMPUTE THE DISTRIBUTION OF H AND THE EXPECTED VALUE UNDER VARIOUS GIBBS TYPE ####
####################################################################################

expected_cl_py <- function(n, sigma, theta, H){
  n <- as.integer(n)
  stopifnot(sigma >= 0, sigma < 1, theta > - sigma, n > 0, H > 1)
  if(H == Inf){
    if(sigma==0) {
    out <- theta * sum(1/(theta - 1 + 1:n))
    } else {
    out <- 1/sigma*exp(lgamma(theta + sigma + n) - lgamma(theta + sigma) - lgamma(theta + n) + lgamma(theta + 1)) - theta/sigma
    }
  } else if(H < Inf){
    if(sigma==0) {
      index <- 0:(n-1)
      out <- H - H*exp(sum(log(index + theta*(1 - 1/H)) - log(theta+ index)))
    }
  }
  return(out)
}


# ----------------------------------------------
# INTERNAL USAGE ONLY
# ----------------------------------------------

lV_py <- function(n, k, sigma, theta){
  if(k==1) return(lgamma(theta + 1) - lgamma(theta + n))
  index <- 1:(k-1)
  sum(log(theta + index*sigma)) + lgamma(theta + 1) - lgamma(theta + n)
}

# ----------------------------------------------
# INTERNAL USAGE ONLY
# ----------------------------------------------

lV_py <- Vectorize(lV_py, vectorize.args = "k")

dclust_py <- function(n, sigma, theta, H, log_scale=FALSE){
  
  probs <- numeric(n) - Inf
  if(H == Inf){
    index <- 1:n
    if(sigma > 0) {
      probs  <- c(lgfactorial(n, sigma)) - index*log(sigma) + lV_py(n,index,sigma,theta)
    } else {
      probs  <- c(lastirling1(n)) + lV_py(n,index,sigma,theta)
    }
  }
  if(H < Inf){
    index <- 1:min(n,H)
    if(sigma == 0){
      probs[index]  <- lfactorial(H) - lfactorial(H-index) + c(lgfactorial_ns(n,-theta/H))[index] + lgamma(theta) - lgamma(theta + n)
    }
  }
  if(log_scale) return(probs)
  return(exp(probs))
}


# ----------------------------------------------
# GNEDIN
# ----------------------------------------------

HGnedin <- function(V, h, gamma=0.5){
  exp(lchoose(V, h) + lgamma(h-gamma) - lgamma(1-gamma) + log(gamma) + lgamma(V+ gamma - h) - lgamma(V +gamma))
}



# ------------------------------------------------------------------------------
#   FUNCTION TO PLOT REPRDERED ADJACENCY MATRICES BASED ON COMMUNITY ASSIGNMENT
# ------------------------------------------------------------------------------
plot_reordered_adjacency <- function(membership, adjacency_matrix) {
  unique_memberships <- unique(membership)
  order_indices <- order(membership)
  reordered_matrix <- adjacency_matrix[order_indices, order_indices]
  #reordered_matrix <- reordered_matrix[,nrow(reordered_matrix):1]
  
  # Determine the aspect ratio based on the dimensions of the reordered matrix
  asp_ratio <- ncol(reordered_matrix) / nrow(reordered_matrix)
  
  # Plot the reordered adjacency matrix using the image function
  image(1:nrow(reordered_matrix), 1:ncol(reordered_matrix), reordered_matrix, 
        col = colorRampPalette(brewer.pal(9, "Greys")[c(1, 8)])(30),
        axes = FALSE, xlab = "", ylab = "", asp = 1, ylim = c(0,dim(adjacency_matrix)[1]))
  
  # Calculate adjusted positions for vertical grid lines
  v_lines <- c(0.5, cumsum(table(membership)) + 0.5)
  v_lines_adj <- v_lines * asp_ratio
  
  # Calculate adjusted positions for horizontal grid lines
  h_lines <- c(0.5, cumsum(table(membership)) + 0.5)
  h_lines_adj <- h_lines / asp_ratio
  
  # Add vertical grid lines at community boundaries
  for (v in v_lines_adj) {
    segments(v, 0.5, v, ncol(reordered_matrix) + 0.5, col = "blue", lty = 1)
  }
  
  # Add horizontal grid lines at community boundaries
  for (h in h_lines_adj) {
    segments(0.5, h, nrow(reordered_matrix) + 0.5, h, col = "blue", lty = 1)
  }
  
  library("gridGraphics")
  fun_grob <- function(){
    grid.echo()
    grid.grab()
  }
  fun_grob_out <- fun_grob()
  pushViewport(viewport(width = 1, height = 1, angle = 270))
  grid.draw(fun_grob_out)
}



plot_reordered_adjacency_hmldm <- function(membership, order_indices, adjacency_matrix) {
  unique_memberships <- unique(membership)
  #order_indices <- order_indices+1
  reordered_matrix <- adjacency_matrix[order_indices, order_indices]
  #reordered_matrix <- reordered_matrix[,nrow(reordered_matrix):1]
  
  # Determine the aspect ratio based on the dimensions of the reordered matrix
  asp_ratio <- ncol(reordered_matrix) / nrow(reordered_matrix)
  
  # Plot the reordered adjacency matrix using the image function
  image(1:nrow(reordered_matrix), 1:ncol(reordered_matrix), reordered_matrix, 
        col = colorRampPalette(brewer.pal(9, "Greys")[c(1, 8)])(30),
        axes = FALSE, xlab = "", ylab = "", asp = 1)
  
  # Calculate adjusted positions for vertical grid lines
  v_lines <- c(0.5, cumsum(table(membership)) + 0.5)
  v_lines_adj <- v_lines * asp_ratio
  
  # Calculate adjusted positions for horizontal grid lines
  h_lines <- c(0.5, cumsum(table(membership)) + 0.5)
  h_lines_adj <- h_lines / asp_ratio
  
  # Add vertical grid lines at community boundaries
  for (v in v_lines_adj[-length(v_lines_adj)]) {
    segments(v, 0.5, v, ncol(reordered_matrix) + 0.5, col = "blue", lty = 1)
  }
  
  # Add horizontal grid lines at community boundaries
  for (h in h_lines_adj[-length(h_lines_adj)]) {
    segments(0.5, h, nrow(reordered_matrix) + 0.5, h, col = "blue", lty = 1)
  }
  
  segments(nrow(reordered_matrix) + 0.5, 0.5, nrow(reordered_matrix) + 0.5, ncol(reordered_matrix) + 0.5, col = "blue", lty = 1)
  segments(0.5, ncol(reordered_matrix) + 0.5, nrow(reordered_matrix) + 0.5, ncol(reordered_matrix) + 0.5, col = "blue", lty = 1)
  
  library("gridGraphics")
  fun_grob <- function(){
    grid.echo()
    grid.grab()
  }
  fun_grob_out <- fun_grob()
  pushViewport(viewport(width = 1, height = 1, angle = 270))
  grid.draw(fun_grob_out)
}

