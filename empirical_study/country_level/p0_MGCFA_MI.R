load("C:/Users/u0159126/wp4_empirical_01052026/country_level_0116_20random/data_subset.RData")
data <- data_subset
#MG-SEM
library(lavaan)
mgcfa_model <- '
  # measurement
  att =~ ST181Q02HA + ST181Q03HA +ST181Q04HA
   mot  =~ ST182Q03HA + ST182Q04HA + ST182Q05HA + ST182Q06HA
   adapt =~ ST216Q01HA + ST216Q02HA + ST216Q03HA + ST216Q04HA + ST216Q05HA + ST216Q06HA
   # compenv =~ ST205Q01HA + ST205Q02HA + ST205Q03HA + ST205Q04HA
   # collenv =~ ST206Q01HA + ST206Q02HA + ST206Q03HA + ST206Q04HA
'

# Fit multigroup SEM treating items as continuous, robust ML
fit1_mg <- sem(mgcfa_model,
               data = data,
               group = "CNTRYID",
               estimator = "MLR",    # robust ML for continuous ordinal-ish data
               std.lv = TRUE)        # identify by fixing latent variances = 1

# Summary: fit measures + standardized estimates
# summary(fit1_mg, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
round(fitmeasures(object = fit1_mg)[c("cfi", "rmsea")], 3)
# cfi rmsea 
# 0.959 0.056 

# metric invariance
fit2_mg <- sem(mgcfa_model,
               data = data,
               group = "CNTRYID",
               estimator = "MLR",    # robust ML for continuous ordinal-ish data
               std.lv = TRUE,
               group.equal = "loadings")        # identify by fixing latent variances = 1
round(fitmeasures(object = fit2_mg)[c("cfi", "rmsea")], 3)
# cfi rmsea 
# 0.953 0.056

##Step1: MM fit with markers##
library(lavaan)
mgcfa_model2 <- '
  # measurement
  att =~ ST181Q04HA + ST181Q02HA + ST181Q03HA 
   mot  =~ ST182Q04HA + ST182Q03HA + ST182Q05HA + ST182Q06HA
   adapt =~ ST216Q03HA + ST216Q01HA + ST216Q02HA + ST216Q04HA + ST216Q05HA + ST216Q06HA
   # compenv =~ ST205Q02HA + ST205Q01HA + ST205Q03HA + ST205Q04HA
   # collenv =~ ST206Q03HA + ST206Q01HA + ST206Q02HA + ST206Q04HA
'
fit1_mgcfa <- sem(mgcfa_model2,
               data = data,
               group = "CNTRYID",
               estimator = "MLR"    # robust ML for continuous ordinal-ish data
               )       

est <- lavInspect(fit1_mgcfa, "est")
load <- lapply(est, "[[", "lambda")

ref <- load[[1]]
idx <- which(ref != 0, arr.ind = TRUE)
length(idx[,1])
# should return 21
loading_vectors <- lapply(seq_len(nrow(idx)), function(i) {
  item   <- rownames(ref)[idx[i, "row"]]
  factor <- colnames(ref)[idx[i, "col"]]
  
  sapply(load, function(x) x[item, factor])
})
names(loading_vectors) <- apply(idx, 1, function(i) {
  paste(rownames(ref)[i["row"]])
})

loading_matrix <- do.call(cbind, loading_vectors)
boxplot(
  loading_matrix,
  outline = TRUE,
  las = 2,
  ylab = "Loading estimate",
  cex.axis = 0.6
)

# partial metric inv
start.time <- Sys.time()
fit2_mgcfa <- cfa(mgcfa_model2,
               data = data,
               group = "CNTRYID",
               estimator = "MLR",
               group.equal = "loadings",
               group.partial = c("att =~ ST181Q02HA",
                                 "mot =~ ST182Q05HA", 
                                 "mot =~ ST182Q06HA",
                                 "adapt =~ ST216Q05HA"))
round(fitmeasures(object = fit2_mgcfa)[c("cfi", "rmsea")], 3)  
end.time <- Sys.time()
diff_time2 <- difftime(end.time, start.time, units='mins')

# cfi rmsea 
# 0.957 0.055 
# Time difference of 0.6152963 mins

##############################MMGSEM_STEP2###################################

sm_model <- '
  # regressions
  mot ~ att 
  adapt ~ att + mot

'
data <- as.data.frame(data)

result <- list()
runtime_min <- c()
for (i in 1:9) {
  start.time <- Sys.time()
  result[[i]] <-  MMGSEM(data, 
                         S1 = mgcfa_model2, S2 = sm_model,
                         group = "CNTRYID", nclus = i, seed = 100, 
                         userStart = NULL, s1_fit = fit2_mgcfa,
                         max_it = 10000L, nstarts = 50, printing = T,
                         partition = "hard")
  end.time <- Sys.time()
  runtime_min[i] <- difftime(end.time, start.time, units='mins')
}

MMGSEM <- list(result = result,
               runtime_min = runtime_min)

save(MMGSEM, file ="MMGSEM_9k_3F.RData")

####model selection####
BIC_G <- sapply(MMGSEM$result, function(x) x$model_sel$BIC$observed$BIC_G)
AIC <- sapply(MMGSEM$result, function(x) x$model_sel$AIC$observed)
nr.par_obs <- sapply(MMGSEM$result, function(x) x$NrPar$Obs.nrpar)
LL_S <- sapply(MMGSEM$result, function(x) x$logLik$obs_loglik)

#BIC/AIC
clus <- c(1:9)
par(mfrow = c(1, 2))
plot(x = clus, y = BIC_G); lines(x = clus, y = BIC_G)#3
plot(x = clus, y = AIC); lines(x = clus, y = AIC)#3

#CHull
CHull(loglik = LL_S, nrpar = nr.par_obs, nsclust = c(1, 9)) #3
hpts_S <- chull(x = nr.par_obs, y = LL_S)
plot(x = nr.par_obs, y = LL_S); lines(x = nr.par_obs[hpts_S], y = LL_S[hpts_S]) #3
[,1]     [,2] [,3]     [,4]
[1,]    1 -1409907  409       NA
[2,]    2 -1409406  413 3.025801
[3,]    3 -1409241  417 1.333696
[4,]    4 -1409116  421 2.382930
[5,]    5 -1409076  425       NA
[6,]    6 -1409020  429       NA
[7,]    7 -1408960  433 2.820774
[8,]    8 -1408948  437       NA
[9,]    9 -1408923  441       NA


clustering <- t(apply(MMGSEM$result[[4]][["posteriors"]], 1, function(x) as.numeric(x == max(x))))
table(data_subset$CNTRYID)
table(data_subset$CNT)
# 7    12    18    21    24    25    28    32    33    34    40    42    43    48    50    55 
# 4475  4402  4434  4808  3886 10634  4954  6338  3461  3599  3343  2467  4502  5542  4454  5620 
# 57    62    63    65 
# 6093  7741  5974  2304 

cbind(clustering, as.numeric(names(table(data_subset$CNTRYID))),
      names(table(data_subset$CNT)))

MMGSEM$result[[4]]$param$beta_ks
uncertainty <- 1 - apply(MMGBSEM$result[[4]][["posteriors"]], 1, max)
meanuncertain <- mean(uncertainty)
# MG: 9.40909e-07; ML: 6.031522e-07; MGB: 6.451639e-05; 2mix: 3.711898e-07
# > round(runtime_min,3)
# [1] 0.248 0.562 0.625 1.044 1.440 1.595 1.886
# [8] 1.818 2.267