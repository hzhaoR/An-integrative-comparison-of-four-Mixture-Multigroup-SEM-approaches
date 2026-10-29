group <- "CNTRYID"
nfactors <- M <- 3
post_est_ordered <- vector(mode = "list", M)
matching_columns <- matching_columns_sd <- list()
latvar <- c("att", "mot", "adapt")
latvar_names <- toupper(latvar)

res_temp <- MGBCFA_output_4non 

for (i in 1:3) {
  
# res_temp[[i]] <- mplusModeler(mod, modelout = paste0("cfa", i, ".inp"), run = 1L)

# col_names <- colnames(res_temp[[i]][["results"]][["savedata"]])
# matched_col <- col_names[toupper(group) == col_names][1] #Mplus group name

# post_est_ordered[[i]] <- res_temp[[i]][["results"]][["savedata"]][order(res_temp[[i]][["results"]][["savedata"]][[matched_col]]), ]

post_est_ordered[[i]] <- res_temp[[i]][["results"]][["savedata"]][
  order(res_temp[[i]][["results"]][["savedata"]][[ncol(res_temp[[i]][["results"]][["savedata"]])]]),
]

N_gs <- as.numeric(table(post_est_ordered[[i]][[ncol(post_est_ordered[[i]])]])) #within-group sample size
ngroups <- length(unique(post_est_ordered[[i]][[ncol(post_est_ordered[[i]])]])) #the number of groups

# Regular expression pattern to match either '%W' or nothing, and allow optional spaces before 'Mean'
pattern <- paste0("^(", paste(latvar_names, collapse = "|"), ")\\s*Mean$")
# Locate columns matching the pattern
matching_columns[[i]] <- grep(pattern, names(post_est_ordered[[i]]), value = TRUE)

# Regular expression pattern to match either '%W' or nothing, and allow optional spaces before 'Mean'
pattern_sd <- paste0("^(", paste(latvar_names, collapse = "|"), ")(%W)?\\s*Standard D*")
# Locate columns matching the pattern
matching_columns_sd[[i]] <- grep(pattern_sd, names(post_est_ordered[[i]]), value = TRUE)
}

# Intermediate Step: the single-indicator approach (factor scores) ---------------------------------
# posterior mean estimates of factor scores
postmeanfs <- do.call(cbind, lapply(seq_along(post_est_ordered), function(i) {
  post_est_ordered[[i]][, matching_columns[[i]], drop = FALSE]
}))
# factor scores
fscores <- data.frame(cbind(postmeanfs, post_est_ordered[[1]][[ncol(post_est_ordered[[1]])]])) 
vars <- paste0("single", c(1:nfactors)) #observed var
colnames(fscores) <- c(vars, "group")
fcovdat <- lapply(split(fscores[, 1:nfactors], fscores$group), cov) # factor score covariances
fcovdat <- array(unlist(fcovdat), dim = c(nfactors, nfactors, ngroups))
postmeanfs_list <- split(fscores[, 1:nfactors], fscores$group)

# posterior standard deviation/variance estimates of factor scores
postvarfs <- cbind((do.call(cbind, lapply(seq_along(post_est_ordered), function(i) {
  post_est_ordered[[i]][, matching_columns_sd[[i]], drop = FALSE]
})))^2, post_est_ordered[[1]][[ncol(post_est_ordered[[1]])]])
colnames(postvarfs) <- c(vars, "group")
postvarfs_list <- split(postvarfs[, 1:nfactors], postvarfs$group)

# factor variances computed based on factor scores and their variances
fvar_compute <- replicate(nfactors, numeric(ngroups), simplify = FALSE)
for (i in 1:nfactors) {
  for (j in 1:ngroups) {
    fvar_compute[[i]][j] <- var(postmeanfs_list[[j]][,i]) + mean(postvarfs_list[[j]][,i])
  }
}

# single-indicator loadings (rhos)
rho <- lapply(1:nfactors, function(i) sapply(1:ngroups, function(j) var(postmeanfs_list[[j]][,i]) / fvar_compute[[i]][j]))
lambda <- array(0, dim = c(nfactors, nfactors, ngroups))
for (j in 1:ngroups) lambda[,,j] <- diag(sapply(rho, `[`, j))
dimnames(lambda)[[1]] <- vars
dimnames(lambda)[[2]] <- latvar  

# single-indicator unique variances
resvar <- lapply(1:ngroups, function(j) sapply(1:nfactors, function(i) fvar_compute[[i]][j] * rho[[i]][j] * (1 - rho[[i]][j])))
theta <- array(unlist(lapply(1:ngroups, function(j) {
  diag(resvar[[j]])  # Create diagonal matrix for each group
})), dim = c(nfactors, nfactors, ngroups))
dimnames(theta)[[1]] <- dimnames(theta)[[2]] <-  vars  

# single-indicator fcov: Croon's formula
cov_eta <-  array(data = 0, dim = c(nfactors, nfactors, ngroups))
for (j in 1:ngroups) {
  cov_eta[,,j] <- solve(lambda[,,j]) %*% (fcovdat[,,j] - theta[,,j]) %*% solve(t(lambda[,,j]))
}
dimnames(cov_eta)[[1]] <- dimnames(cov_eta)[[2]] <-  latvar

# single-indicator data (center the data per group, so that the mean for all variables in each group is 0)
centered <- fscores
g_name <- as.character(unique(fscores[, "group"]))
group.idx <- match(fscores[,"group"], g_name)
group.sizes <- tabulate(group.idx)
group.means <- rowsum.default(as.matrix(fscores[,vars]),
                              group = group.idx, reorder = FALSE,
                              na.rm = FALSE)/group.sizes
centered[,vars] <- fscores[,vars] - group.means[group.idx, ,drop = FALSE]

# Get sample covariance matrix per group (used later)
S_unbiased <- lapply(X = unique(centered[, "group"]), FUN = function(x) {cov(centered[centered[, "group"] == x, vars])})

# save output from mplus
# res_temp
mplus_parameter <- lapply(res_temp, function(r) r[["results"]][["parameters"]][["unstandardized"]])

s1output_mbcfa_4non <- list(
  single_data  = centered,
  lambda_gs    = lambda,
  theta_gs     = theta,
  cov_eta      = cov_eta,
  ngroups      = ngroups,
  N_gs         = N_gs,
  S_unbiased   = S_unbiased,                     
  vars         = vars, #single-indicator names
  lat_var      = latvar, #latent variable names
  
  mplus_output = list(mplus_parameter = mplus_parameter,
                      mplus_fs        = post_est_ordered)
  
)
save(s1output_mbcfa_4non, file = "p2_s1output_mbcfa_4non0.05_3F.RData")

sm_model <- '
  # regressions
  mot ~ att 
  adapt ~ att + mot

'
#Step 2 with nclus=k
library(lavaan)
result <- list()
runtime_min <- c()

for (i in 1:9) {
  start.time <- Sys.time()
result[[i]] <- Mix_Step2(s1out = s1output_mbcfa_4non, step2model = sm_model,
                         nclus = i, seed = 100, userStart = NULL,
                         max_it = 10000L, nstarts = 50, printing = T,
                         partition = "hard", Endo2Cov = TRUE, allG = TRUE, fit = "factors")
end.time <- Sys.time()
runtime_min[i] <- difftime(end.time, start.time, units='mins')
}

MMGBSEM <- list(result = result,
               runtime_min = runtime_min)
# [1] 0.06014053 0.56153653 0.76318233 0.83751747 0.91373097
# [6] 1.58906175 1.99002382 2.63461930 2.97107232

save(MMGBSEM, file ="MMGBSEM_9k_3F.RData")

##MGB
s2output_4k_post <- cbind(MMGBSEM$result[[4]][["posteriors"]], s1output_mbcfa_4non$N_gs)

BIC_G <- sapply(MMGBSEM$result, function(x) x$model_sel$BIC$observed$BIC_G)

AIC <- sapply(MMGBSEM$result, function(x) x$model_sel$AIC$observed)

nr.par_obs <- sapply(MMGBSEM$result, function(x) x$NrPar$Obs.nrpar)

LL_S <- sapply(MMGBSEM$result, function(x) x$logLik$obs_loglik)


#BIC/AIC
clus <- c(1:9)
par(mfrow = c(1, 2))
plot(x = clus, y = BIC_G); lines(x = clus, y = BIC_G)#3
plot(x = clus, y = AIC); lines(x = clus, y = AIC)#3

#CHull
CHull(loglik = LL_S, nrpar = nr.par_obs, nsclust = c(1, 9)) #3
hpts_S <- chull(x = nr.par_obs, y = LL_S)
plot(x = nr.par_obs, y = LL_S); lines(x = nr.par_obs[hpts_S], y = LL_S[hpts_S]) #3
[,1]      [,2] [,3]     [,4]
[1,]    1 -239026.4  183       NA
[2,]    2 -238527.7  187 2.226046
[3,]    3 -238303.8  191 1.206819
[4,]    4 -238118.1  195 2.892185
[5,]    5 -238054.0  199 1.308050
[6,]    6 -238004.9  203 2.686726
[7,]    7 -237988.9  207       NA
[8,]    8 -237968.4  211 1.283839
[9,]    9 -237954.2  215       NA

