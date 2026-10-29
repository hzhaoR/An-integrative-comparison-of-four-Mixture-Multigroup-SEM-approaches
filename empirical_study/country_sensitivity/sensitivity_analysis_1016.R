##2026.10.16 sensitivity analysis##
# empirical data 9775*21 var
load("C:/Users/u0159126/wp4_empirical_01052026/country_level_0116_20random/data_subset.RData")
# data <- data_subset
data <- as.data.frame(data_subset[,c(1,4:24)])
library(MplusAutomation)

start_time <- Sys.time()

# Define variables and Mplus models for each factor 
var_names <- list("ST181Q02HA ST181Q03HA ST181Q04HA",
                  "ST182Q03HA ST182Q04HA ST182Q05HA ST182Q06HA",
                  "ST216Q01HA ST216Q02HA ST216Q03HA ST216Q04HA ST216Q05HA ST216Q06HA")
# "ST205Q01HA ST205Q02HA ST205Q03HA ST205Q04HA",
# "ST209Q01HA ST209Q02HA ST209Q03HA ST209Q04HA")

# latent variable names
latvar_names <- c("att", "mot", "adapt") #, "compenv", "collenv")

# Mplus models  
mplus_models_4non_ml <- list("%within%
att BY ST181Q04HA;
s1 | att BY ST181Q02HA;
s2 | att BY ST181Q03HA;
logatt |att;
logST181Q02HA |ST181Q02HA;
logST181Q03HA |ST181Q03HA;
logST181Q04HA |ST181Q04HA;
[ST181Q02HA-ST181Q04HA@0]",
                          
                          "%within%
mot BY ST182Q04HA ST182Q03HA;
s6 | mot BY ST182Q05HA;
s7 | mot BY ST182Q06HA;
logmot |mot;
logST182Q03HA |ST182Q03HA;
logST182Q04HA |ST182Q04HA;
logST182Q05HA |ST182Q05HA;
logST182Q06HA |ST182Q06HA;
[ST182Q03HA-ST182Q06HA@0]",
                          
                          "%within%
adapt BY ST216Q03HA ST216Q02HA ST216Q06HA ST216Q01HA ST216Q04HA ST216Q05HA;

logadapt |adapt;
logST216Q01HA | ST216Q01HA;
logST216Q02HA | ST216Q02HA;
logST216Q03HA | ST216Q03HA;
logST216Q04HA | ST216Q04HA;
logST216Q05HA | ST216Q05HA;
logST216Q06HA | ST216Q06HA;
[ST216Q01HA-ST216Q06HA@0]"
)


set.seed(100)
#Step 1
s1output_4non_ml <- try(MixML_Step1(data = data, var_names = var_names, latvar_names = latvar_names,
                                 step1model = mplus_models_4non_ml, group = "CNTRYID"))
save(s1output_4non_ml, file = "ML-CFA_4non_3F_ml.RData")

end_time <- Sys.time()
diff_time <- difftime(end_time, start_time, units='mins') 
# Time difference of 10.78228 mins

########################################3MMLSEM_STEP2#########################
sm_model <- '
  # regressions
  mot ~ att 
  adapt ~ att + mot

'
#Step 2 with nclus=k
MMLSEM_s2output <- list()
diff_time <- c()

for (i in 1:9) {
  start.time <- Sys.time()
  MMLSEM_s2output[[i]] <- Mix_Step2(s1out = s1output_4non_ml, step2model = sm_model,
                                    nclus = i, seed = 100, userStart = NULL,
                                    max_it = 10000L, nstarts = 50, printing = T,
                                    partition = "hard", Endo2Cov = TRUE, allG = TRUE, fit = "factors")
  end.time <- Sys.time()
  diff_time[i] <- difftime(end.time, start.time, units='mins')
  
}
diff_time
# 0.06811926 0.36192840 0.52071987 0.77099245 0.80927785 1.03524808 1.29601536
# [8] 1.61392076 1.86504147
save(diff_time, file = "MMLSEM_step2_runtime_ml.RData")
save(MMLSEM_s2output, file = "MixMLSEM_s2output_4non_ml.RData")

###model selection############
BIC_G <- sapply(MMLSEM_s2output, function(x) x$model_sel$BIC$observed$BIC_G)

AIC <- sapply(MMLSEM_s2output, function(x) x$model_sel$AIC$observed)

nr.par_obs <- sapply(MMLSEM_s2output, function(x) x$NrPar$Obs.nrpar)

LL_S <- sapply(MMLSEM_s2output, function(x) x$logLik$obs_loglik)


#BIC/AIC
clus <- c(1:9)
par(mfrow = c(1, 2))
plot(x = clus, y = BIC_G[1:9]); lines(x = clus, y = BIC_G[1:9])#3
plot(x = clus, y = AIC[1:9]); lines(x = clus, y = AIC[1:9])#3

#CHull
CHull(loglik = LL_S[1:9], nrpar = nr.par_obs[1:9], nsclust = c(1, 9)) #3
hpts_S <- chull(x = nr.par_obs[1:9], y = LL_S[1:9])
plot(x = nr.par_obs[1:9], y = LL_S[1:9]); lines(x = nr.par_obs[hpts_S], y = LL_S[hpts_S]) #3
[,1]      [,2] [,3]     [,4]
[1,]    1 -238120.3  183       NA
[2,]    2 -237608.6  187 2.482528
[3,]    3 -237402.5  191 1.140563
[4,]    4 -237221.8  195 2.577623
[5,]    5 -237151.7  199 1.306175
[6,]    6 -237098.1  203 2.170482
[7,]    7 -237073.3  207 2.561204
[8,]    8 -237063.9  211       NA
[9,]    9 -237054.0  215       NA

clustering <- t(apply(MMLSEM_s2output[[4]][["posteriors"]], 1, function(x) as.numeric(x == max(x))))
table(data_subset$CNTRYID)
# 7    12    18    21    24    25    28    32    33    34    40    42    43    48    50    55 
# 4475  4402  4434  4808  3889 10934  4954  9338  3491  3599  3343  2497  4502  5542  4454  5920 
# 57    92    93    95 
# 9093  7741  5974  2304 
[,1] [,2] [,3] [,4]
[1,]    0    0    1    0
[2,]    0    1    0    0
[3,]    1    0    0    0
[4,]    0    1    0    0
[5,]    1    0    0    0
[6,]    0    0    0    1
[7,]    1    0    0    0
[8,]    1    0    0    0
[9,]    0    0    1    0
[10,]    0    0    0    1
[11,]    0    0    1    0
[12,]    1    0    0    0
[13,]    0    0    1    0
[14,]    0    0    0    1
[15,]    0    1    0    0
[16,]    0    0    1    0
[17,]    1    0    0    0
[18,]    1    0    0    0
[19,]    0    0    1    0
[20,]    0    1    0    0
cbind(clustering, as.numeric(names(table(data_subset$CNTRYID))))

##MIxMix
output2 <-mixmgfa(data = centered, N_gs = N_gs, nfactors = nfactors, 
                  cluster.spec = c("loadings"), nsclust = c(1:9), 
                  maxiter = 100000, nruns = 100, design = loading_matrix, 
                  invar_loadings = c(0,1,1,
                                     1,1,0,1,
                                     1,1,1,1,1,1))
plot(output2$overview, AIC = 1)

save(output2, file = "2Mix_MM_100starts_2non.RData")
# output2$overview
# nr of clusters   loglik nrpars   BIC_N   BIC_G     AIC screeratios convergence nr.activated.constraints
# [1,]              1 -1409714    650 2826905 2821375 2820728          NA           1                        0
# [2,]              2 -1409505    653 2826521 2820966 2820316          NA           1                        0
# [3,]              3 -1409264    656 2826075 2820494 2819841    5.001454           1                        0
# [4,]              4 -1409219    659 2826019 2820413 2819757    2.307078           1                        0
# [5,]              5 -1409200    662 2826015 2820383 2819724    2.010647           1                        0
# [6,]              6 -1409193    665 2826036 2820378 2819716          NA           1                        0
# [7,]              7 -1409191    668 2826066 2820383 2819718          NA           1                        0
# [8,]              8 -1409187    671 2826092 2820383 2819715          NA           1                        0
# [9,]              9 -1409161    674 2826075 2820341 2819670          NA           1                        0
start_time <- Sys.time()
s1output <- list()
for (i in 1:9) {

  s1output[[i]] <- MixMix_Step1(data, step1model = mm_model, group = "CNTRYID", MM.cluster.spec=c("loadings"),
                                MM.nclus = i, MM.maxiter = 100000, MM.nruns = 100, MM.design = loading_matrix,
                                # invar_loadings = c(0,0,1,
                                #                    0,1,0,0,
                                #                    0,0,1,0,0,0),
                                invar_loadings = c(0,1,1,
                                                   1,1,0,1,
                                                   1,1,1,1,1,1),
                                markers = markers, seed = 100)
}
end_time <- Sys.time()
diff_time <- difftime(end_time, start_time, units='mins')
# 1.306947 mins
save(s1output, file = "s1output_2non.RData")




##MixMG-BSEM


