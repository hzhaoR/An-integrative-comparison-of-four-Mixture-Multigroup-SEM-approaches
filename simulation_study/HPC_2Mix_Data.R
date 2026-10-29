library(MASS)
library(matrixcalc)
library(lavaan)

Args <- commandArgs(TRUE)
i <- as.numeric(Args[1])
load("Design_cluster_0.3.RData")

#' DataGlavaan#' DataGeneration
#'
#' Generates simulated data according to the manipulated conditions. 
#'
#' INPUT: Arguments required by the function
#' @param model: full model in lavaan syntax (string).
#' @param SM_nclus: number of clusters of the data. 
#' @param ngroups: number of groups for the generated data.
#' @param N_g: number of observations per group.
#' @param SM_coeff: regression coefficients of cluster 1. Must a numeric vector of length 3.
#' @param balance: can take values "balanced" and "unbalanced". Determines if the clusters will have the same
#'        number of groups.
#' @param reliability: reliability per item (manipulates ratio between loading and unique variance).
#' @param NonInvSize: size of the loading non-invariance (for the simulation: 0.2 and 0.4).
#' @param NonInvItems: number of items affected by the non-invariance (for now, only 2).
#' @param NonInvG: proportion of groups affected by the non-invariance (e.g., 0.25, 0.5).
#' @param randomVarX: define the variance of the exogenous variable. If TRUE, then random values will be defined
#'        for all groups. If FALSE, the variance will be 1 for all groups.
#'
#' OUTPUT
#' @return SimData: generated data
#' 

#   # exog1 -> endog2 = B1
#   # exog1 -> endog1 = B2
#   # endog1 -> endog2 = B3
#   # endog1 -> endog2 = B4


DataGeneration <- function(model, ngroups, N_g,
                           smallG_proportion, smallG_samplesize, 
                           MM_nclus, MM_balance,
                           SM_nclus, SM_balance,
                           SM_coeff, MM_loading, MM_uniqvar,
                           reliability = "low", emp){
  
  # Get number of variables
  par_table <- lavaanify(model)
  lat_var <- lavNames(par_table, "lv")
  obs_var <- lavNames(par_table, "ov")
  m <- length(lat_var)
  p <- length(obs_var)
  
  # Identify type of latent variable
  endog1 <- lat_var[(lat_var %in% par_table$rhs[which(par_table$op == "~")]) &
                      (lat_var %in% par_table$lhs[which(par_table$op == "~")])]
  endog2 <- lat_var[!c(lat_var %in% par_table$rhs[which(par_table$op == "~")]) &
                      (lat_var %in% par_table$lhs[which(par_table$op == "~")])]
  endog <- c(endog1, endog2)
  exog <- lat_var[!c(lat_var %in% endog)]
  
  # Reorder latent variables
  lat_var <- c(exog, endog)
  
  # Give names to the regression parameters (only for better understanding):
  #   # exog1 -> endog2 = B1
  #   # exog1 -> endog1 = B2
  #   # endog1 -> endog2 = B3
  #   # endog1 -> endog2 = B4
  B1 <- numeric(SM_nclus)
  B2 <- numeric(SM_nclus)
  B3 <- numeric(SM_nclus)
  B4 <- numeric(SM_nclus)
  
  # Get regression parameters for each cluster
  # Cluster 1
  B1[1] <- 0
  B2[1] <- SM_coeff
  B3[1] <- SM_coeff
  B4[1] <- SM_coeff
  
  # Cluster 2
  B1[2] <- SM_coeff
  B2[2] <- 0
  B3[2] <- SM_coeff
  B4[2] <- SM_coeff
  
  if (SM_nclus == 4){
    # Cluster 3
    B1[3] <- SM_coeff
    B2[3] <- SM_coeff
    B3[3] <- 0
    B4[3] <- SM_coeff
    
    # Cluster 4
    B1[4] <- SM_coeff
    B2[4] <- SM_coeff
    B3[4] <- SM_coeff
    B4[4] <- 0
  }
  # Generate the structural parameters
  # beta is cluster-specific. Only SM_nclus matrices needed
  beta <- array(data = 0, dim = c(m, m, SM_nclus), dimnames = list(lat_var, lat_var))
  colnames(beta) <- rownames(beta) <- lat_var
  
  # Initialize
  for(k in 1:SM_nclus){
    # beta
    beta[endog2, exog[1], k] <- B1[k]
    beta[endog1, exog[1], k] <- B2[k]
    beta[endog1, exog[2], k] <- B3[k]
    beta[endog2, endog1, k] <- B4[k]
    
  }
  # psi is group- and cluster-specific. ngroups matrices are needed.
  # First, let's define the cluster-specific part
  # psi_k <- array(data = diag(m), dim = c(m, m, SM_nclus), dimnames = list(lat_var, lat_var))
  # Generate enough psi matrices
  psi_g <- array(data = diag(m), dim = c(m, m, ngroups), dimnames = list(lat_var, lat_var))
  # Now, let's define the group-specific part
  # Manipulate variance any way you like
  # randomVarX == T
    exog_var1 <- runif(n = ngroups, min = 0.75, max = 1.25)
    exog_var2 <- runif(n = ngroups, min = 0.75, max = 1.25)
    exog_cov <- runif(n = ngroups, min = -0.30, max = 0.30)
    endo_var1 <- runif(n = ngroups, min = 0.75, max = 1.25)
    endo_var2 <- runif(n = ngroups, min = 0.75, max = 1.25)

    # Get a SM_cluster label for each group (randomized_SM_GperK)
    # Example: If balanced, G = 6, K = 2. randomized_SM_GperK = 111222
    if (SM_balance == "balanced"){
      SM_GperK <- rep(x = 1:SM_nclus, each = (ngroups/SM_nclus))
      randomized_SM_GperK <- sample(SM_GperK)
    } else if (SM_balance == "unbalanced"){
      smallest <- ngroups/(2+SM_nclus) #3:1; 3:1:1:1
      largest <- 3*smallest # largest and smaller clusters
      largest_cluster <- sample(1:SM_nclus, 1)  # Randomly choose which cluster will be largest
      
      SM_nclusk <- 1:SM_nclus
      SM_GperK <- c(rep(x = largest_cluster, times = largest), rep(x = SM_nclusk[-largest_cluster], each = smallest))
      #GperK <- c(rep(x = 1, times = largest), rep(x = 2:nclus, each = smallest))
      randomized_SM_GperK <- sample(SM_GperK)
    }
    
    # Insert the corresponding group- and cluster-specific parts of psi
    for(g in 1:ngroups){
      # Insert group-specific parts
      psi_g[exog[1], exog[1], g] <- exog_var1[g] # Group-specific
      psi_g[exog[2], exog[2], g] <- exog_var2[g] # Group-specific
      psi_g[exog[1], exog[2], g] <- exog_cov[g] # Group-specific
      psi_g[exog[2], exog[1], g] <- exog_cov[g] # Group-specific
      
      # Insert the group-and-cluster-specific parts (they use the cluster-specific betas)
      # For the endogenous part, start from the total var (endog_var) and subtract the explained by the regression
      psi_g[endog1, endog1, g] <- endo_var1[g] - ((B2[randomized_SM_GperK[g]]^2 * exog_var1[g]) + 
                                                    (B3[randomized_SM_GperK[g]]^2 * exog_var2[g]) + 
                                                    (2 * B2[randomized_SM_GperK[g]] * B3[randomized_SM_GperK[g]] * exog_cov[g])) # Group-specific
      
      psi_g[endog2, endog2, g] <- endo_var2[g] - ((B1[randomized_SM_GperK[g]]^2 * exog_var1[g]) + 
                                                    (B4[randomized_SM_GperK[g]]^2 * endo_var1[g]) + 
                                                    (2 * B1[randomized_SM_GperK[g]] * B4[randomized_SM_GperK[g]] * ((B2[randomized_SM_GperK[g]] * exog_var1[g]) + (B3[randomized_SM_GperK[g]] * exog_cov[g])))) # Group-specific
      
    }  #'psi_g' It contains both total (for exogenous) and residual (for endogenous) variances and covariances of 'eta'.
    
    # Create cov matrix as: ----
    I <- diag(m) # Identity matrix
    cov_eta <- array(data = 0, dim = c(m, m, ngroups), dimnames = list(lat_var, lat_var))
    # browser()
    for(g in 1:ngroups){
      cov_eta[, , g] <- solve(I - beta[, , randomized_SM_GperK[g]]) %*% psi_g[, , g] %*% solve(t(I - beta[, , randomized_SM_GperK[g]]))
    }
    
    # Make sure that the average covariance per cluster is 1
    # Add (or remove) the necessary values
    endo_variances <- as.data.frame(cbind(cov_eta[endog1, endog1, ], cov_eta[endog2, endog2, ], randomized_SM_GperK))
    colnames(endo_variances) <- c("F3", "F4", "Cluster")
    
    for(k in 1:SM_nclus){
      this_k <- endo_variances$Cluster == k
      rsd_F3 <- endo_variances$F3[this_k]
      rsd_F4 <- endo_variances$F4[this_k]
      
      needed_F3 <- mean(1 - rsd_F3)
      needed_F4 <- mean(1 - rsd_F4)
      
      endo_variances$F3[this_k] <- endo_variances$F3[this_k] + needed_F3
      endo_variances$F4[this_k] <- endo_variances$F4[this_k] + needed_F4
      
      cov_eta[endog1, endog1, this_k] <- endo_variances$F3[this_k]
      cov_eta[endog2, endog2, this_k] <- endo_variances$F4[this_k]
    }
    
  # Generate the measurement model parameters
  # Lambda (depends on reliability lvl)
  if (reliability == "low"){
    load <- .4
  } else if (reliability == "high"){
    load <- .6
  }
  
  loadings <- sqrt(load) #high: 0.7745967; low: 0.6324555
  
  # Non-invariance
  # Create Invariant Lambda
  Lambda <- array(matrix(data = rep(x = c(1, rep(loadings, 4), rep(0, p)), times = m)[1:(p*m)], nrow = p, ncol = m), dim = c(p, m, MM_nclus))

  if (MM_nclus==2) {
    Lambda[2,1,1] <- Lambda[7,2,1] <- Lambda[3,1,2] <- Lambda[8,2,2] <- loadings - MM_loading 
    Lambda[3,1,1] <- Lambda[8,2,1] <- Lambda[2,1,2] <- Lambda[7,2,2] <- loadings + MM_loading  
  }
  
  if (MM_nclus==4){
    Lambda[2,1,1] <- Lambda[7,2,1] <- Lambda[3,1,2] <- Lambda[8,2,2] <- loadings - MM_loading 
    Lambda[3,1,1] <- Lambda[8,2,1] <- Lambda[2,1,2] <- Lambda[7,2,2] <- loadings + MM_loading  
    Lambda[12,3,3] <- Lambda[17,4,3] <-  Lambda[13,3,4] <- Lambda[18,4,4] <- loadings - MM_loading 
    Lambda[13,3,3] <- Lambda[18,4,3] <-  Lambda[12,3,4] <- Lambda[17,4,4] <- loadings + MM_loading 
  }
  
  # Get a MM_cluster label for each group (randomized_MM_GperK)
  # Example: If balanced, G = 6, K = 2. randomized_MM_GperK = 111222
  if (MM_balance == "balanced"){
    MM_GperK <- rep(x = 1:MM_nclus, each = (ngroups/MM_nclus))
    randomized_MM_GperK <- sample(MM_GperK)
  } else if (MM_balance == "unbalanced"){
    smallest <- ngroups/(2+MM_nclus) #3:1; 3:1:1:1
    largest <- 3*smallest # largest and smaller clusters
    largest_cluster <- sample(1:MM_nclus, 1)  # Randomly choose which cluster will be largest
    
    MM_nclusk <- 1:MM_nclus
    MM_GperK <- c(rep(x = largest_cluster, times = largest), rep(x = MM_nclusk[-largest_cluster], each = smallest))
    #GperK <- c(rep(x = 1, times = largest), rep(x = 2:nclus, each = smallest))
    randomized_MM_GperK <- sample(MM_GperK)
  }
  
  # Generate Theta per group depending on reliability
  # Generate sample covariance matrix (sigma) per groups
  Sigma <- array(data = 0, dim = c(p, p, ngroups))

  if (MM_uniqvar == "cluster"){
  Theta <- array(data = 0, dim = c(p, p, MM_nclus))
  for (k in 1:MM_nclus){
    Theta[, , k] <- diag(runif(n = p, min = ((1 - load) - .1), max = ((1 - load) + .1))) #0.5-0.7
  }
  for (g in 1:ngroups){
    Sigma[, , g] <- Lambda[,,randomized_MM_GperK[g]] %*% cov_eta[, , g] %*% t(Lambda[,,randomized_MM_GperK[g]]) + Theta[, , randomized_MM_GperK[g]]
  }
  } else if (MM_uniqvar == "group"){
    Theta <- array(data = 0, dim = c(p, p, ngroups))
    for (g in 1:ngroups){
    Theta[, , g] <- diag(runif(n = p, min = ((1 - load) - .1), max = ((1 - load) + .1))) #0.5-0.7
    Sigma[, , g] <- Lambda[,,randomized_MM_GperK[g]] %*% cov_eta[, , g] %*% t(Lambda[,,randomized_MM_GperK[g]]) + Theta[, , g]
    }
  }
  ####################################combination of small and large groups#######################################################
  #small sample size (per cluster) 12.1.2023
  # Create a list to store selected groups
  selected_groups <- list()
  
  # Iterate over unique clusters in GperK
  unique_clusters <- unique(randomized_MM_GperK)
  
  for (cluster in unique_clusters) {
    # Get indices of groups in the current cluster
    cluster_indices <- which(randomized_MM_GperK == cluster)
    
    # Calculate the number of groups to select in this cluster
    num_groups_to_select <- ceiling(smallG_proportion * length(cluster_indices))
    
    # Sample without replacement from the cluster indices
    selected_indices <- sample(cluster_indices, size = num_groups_to_select, replace = FALSE)
    
    # Add the selected indices to the list
    selected_groups[[as.character(cluster)]] <- selected_indices
  }
  
  # Flatten the list to get the selected group indices
  selected_group_indices <- as.numeric(unlist(selected_groups))
  
  #smallGIdx <- sample(x = 1:ngroups, size = smallG_proportion*ngroups, replace = F)
  N_g_ngroups <- rep(x = N_g, each = ngroups) # Group variable
  N_g_ngroups[selected_group_indices] <- smallG_samplesize
  
  
  # Data Generation final
  # For now, mu would be 0 as we are only interested in centered variables
  SimData <- c()
  for(g in 1:ngroups){
    tmp <- mvrnorm(n = N_g_ngroups[g], mu = rep(0, p), Sigma = Sigma[, , g], empirical = emp)
    SimData <- rbind(SimData, tmp)
  }
  
  # Add the final labels
  group <- rep(x = c(1:ngroups), times = N_g_ngroups) # Group variable
  SimData <- cbind(SimData, group)
  colnames(SimData) <- c(obs_var, "group")
  
  # Return data
  return(list(SimData = SimData, 
                  Sigma = Sigma,
              # exog_vars = list(exog_var1, exog_var2), 
              # exog_cov = exog_cov, 
              # endo_vars = list(endo_var1, endo_var2),
              cov_eta = cov_eta,
              psi_g = psi_g,
              Theta = Theta,
              Lambda = Lambda, #cluster-specific
              beta = beta,
              randomized_MM_GperK = randomized_MM_GperK,
              randomized_SM_GperK = randomized_SM_GperK,
              N_g_ngroups=N_g_ngroups))
}

for (i in 1:400) {
  
k=0
SimData <- NULL # Initialize SimData to NULL
while (is.null(SimData)) { # Repeat until SimData is not NULL
  
  set.seed(as.numeric(paste0(i, k))) # Set seed to i + k    
  SimData <- tryCatch({
    DataGeneration(model = Design[i, "model"],
                   ngroups = Design[i, "ngroups"], #24
                   N_g = Design[i, "N_g"], #100
                   
                   smallG_samplesize = Design[i, "smallG_samplesize"], 
                   smallG_proportion = Design[i, "smallG_proportion"], 
                   MM_nclus = Design[i, "MM_nclus"], #2
                   MM_balance = Design[i, "MM_balance"], #balanced
                   SM_nclus = Design[i, "SM_nclus"], #2
                   SM_balance = Design[i, "SM_balance"], #balanced
                   
                   SM_coeff = Design[i, "SM_coeff"], #0.4
                   MM_loading = Design[i, "MM_loading"], #-0.2
                   MM_uniqvar = Design[i, "MM_uniqvar"], #group
                   reliability = Design[i, "reliability"],
                   emp = Design[i, "emp"])
    
  }, error = function(e) {
    # Code to handle the error, such as printing a message or logging the error
    cat("Error in data ", i, "k=", k, ":", e$message, "\n")
    
  })
  k <- k+1
}

# Save data for replication purposes
save(SimData, file = paste0("data/Data", i, ".RData"))

}
