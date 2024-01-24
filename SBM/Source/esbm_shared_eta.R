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
# Y_list = List of VxV symmetric adjacency matrices for each layer
# seed = Random seed for reproducibility
# N_iter = Number of iterations
# prior = String in c("DP","PY","DM","GN")
# z_init_list = List of V-vector initialization assignments for each layer
# a, b = Parameters of the Beta prior on the thetas
# alpha_PY, sigma_PY = Parameters for PY prior (if prior=="PY")
# beta_DM, H_DM = Parameters for DM prior (if prior=="DM")
# gamma_GN = Parameter for GN prior (if prior=="GN")
# x = V-vector of categorical covariates
# alpha_xi = C-vector of parameters for a Dirichlet distribution (if x is specified)
# eta = Shared link probability for the shared eta parameter approach
# eta_list = List of layer-specific link probabilities for the layer-specific eta parameter approach

esbm_shared_eta <- function(Y_list, seed, N_iter, prior, z_init_list, a=1, b=1,
                            alpha_PY=NA, sigma_PY=NA, beta_DM=NA, H_DM=NA, gamma_GN=NA, x=NULL, alpha_xi=NULL, eta){
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
  l <- length(Y_list)  # Number of layers
  V <- nrow(Y_list[[1]])
  
  # Initialization for each layer
  Z_list <- lapply(z_init_list, vec2mat)
  
  # [ii] a vector of length V containing the cluster label for each node, more compact to store;
  # such vectors for all iterations are packed in a VxN_iter matrix z_post, 
  # s.t. z_post[v,t]=h if node v is in cluster h at iteration t
  # Note: matrix z_post takes less memory than a list of N_iter matrices Z
  z_post <- matrix(NA,V,N_iter)
  
  
  # Create the matrix with block connections
  temp   <- Y_list%*%Z_lsit
  m_full <- t(Z_list)%*%temp - diag(0.5*colSums(temp*Z_list),ncol(Z_list))
  
  
  # ----------------------------------------------
  # Beginning of the Gibbs sampler
  # ----------------------------------------------
  
  for (t in 1:N_iter){
    for (i in 1:l){  # Loop over layers
      Y <- Y_list[[i]]
      z_init <- z_init_list[[i]]
      
      # remove empty clusters and
      # if the cluster containing node v has no other node, discard it as well
      if(ncol(Z_list) > 1){
        nonempty_v <- which(colSums(Z_list[-v,]) > 0)  
        Z_list <- Z_list[, nonempty_v]
        if (length(nonempty_v)==1){Z <- matrix(Z_list,V,1)}
        
        # Reduce the dimensions of the m_full matrix
        m_full <- matrix(m_full[nonempty_v,nonempty_v],ncol(Z_list),ncol(Z_list))
      } 
      
      # H = number of active clusters
      H   <- ncol(Z_list)
      Z_v <- Z_list[-v,]
      
      # v_minus = number of nodes in each cluster, excluding node v
      if (H==1){v_minus <- sum(Z_list[-v])} else {v_minus <- colSums(Z_v)}
      
      # r_v = number of edges from node v to each cluster (no self-loops)
      r_v         <- crossprod(Z_v, Y[-v,v])
      h_v         <- which(Z_list[v,] > 0)
      
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
      log_p <- log_addit + log(urn(v_minus)) + c(log_lhds_old, log_lhd_new) + log(eta)
      max_log_p <- max(log_p)
      p <- exp(log_p - max_log_p)
      p <- p / sum(p)
      
      # ----------------------------------------------
      # Sampling the indicator
      # ----------------------------------------------
      
      new_sample <- which(rmultinom(1,1,p) > 0)
      
      # ----------------------------------------------
      # Adjusting Z, H, r_v and m
      # ----------------------------------------------
      
      if(length(h_v) == 1){
        Z_list[v,h_v] <- 0
      }
      
      # If a new value is sampled, increase the dimension of m_full
      if(new_sample== H+1){
        Z_list         <- cbind(Z_list,rep(0,V))
        Z[v,new_sample] <- 1
        m               <- rbind(cbind(m,0),0)
        H               <- H + 1
        r_v             <- crossprod(Z_list[-v,], Y[-v,v])
      } 
      else {
        Z_list[v, new_sample] <- 1
        }
      
      # Updating m_full
      resid2              <- matrix(0,H,H)
      resid2[,new_sample] <- r_v; 
      resid2[new_sample,] <- r_v
      m_full              <- m + resid2
    }
    
    # store cluster assignments at time t in matrix z_post s.t.
    # z_post[v,t]=h if node v is in cluster h at iteration t
    z_post[,t] <- Z_list %*% c(1:ncol(Z_list))
    
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
  
  M <- t(z)%*%Y%*%z
  diag(M) <- diag(M)/2
  Tot <- t(z)%*%matrix(1,V,V)%*%z
  diag(Tot) <- (diag(Tot)-table(memb))/2
  Mbar <- Tot-M
  a_n <- lowerTriangle(M,diag=TRUE)+a
  b_bar_n <- lowerTriangle(Mbar,diag=TRUE)+b
  
  theta <- rbeta(length(a_n),a_n,b_bar_n)
  Theta <- matrix(0,H,H)
  Theta[lower.tri(Theta,diag=TRUE)] <- theta
  Theta <- Theta+t(Theta)
  diag(Theta) <- diag(Theta)/2
  edge_prob <- z%*%Theta%*%t(z)
  
  LL <- dbinom(lowerTriangle(Y,diag=FALSE), size=1, prob=lowerTriangle(edge_prob,diag=FALSE),log=TRUE)
  
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
# GNEDIN
# ----------------------------------------------

HGnedin <- function(V, h, gamma=0.5){
  exp(lchoose(V, h) + lgamma(h-gamma) - lgamma(1-gamma) + log(gamma) + lgamma(V+ gamma - h) - lgamma(V +gamma))
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

# ... (rest of the code)
