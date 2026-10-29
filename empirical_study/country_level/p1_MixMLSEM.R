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
mplus_models_4non <- list("%within%
att BY ST181Q04HA ST181Q03HA;
s1 | att BY ST181Q02HA;
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
adapt BY ST216Q03HA ST216Q02HA ST216Q06HA ST216Q01HA ST216Q04HA;
s12 | adapt BY ST216Q05HA;

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
s1output_4non <- try(MixML_Step1(data = data, var_names = var_names, latvar_names = latvar_names,
                        step1model = mplus_models_4non, group = "CNTRYID"))
save(s1output_4non, file = "ML-CFA_4non_3F.RData")

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
MMLSEM_s2output[[i]] <- Mix_Step2(s1out = s1output_4non, step2model = sm_model,
                      nclus = i, seed = 100, userStart = NULL,
                      max_it = 10000L, nstarts = 50, printing = T,
                      partition = "hard", Endo2Cov = TRUE, allG = TRUE, fit = "factors")
end.time <- Sys.time()
diff_time[i] <- difftime(end.time, start.time, units='mins')

}
diff_time
# [1] 0.05820976 0.29325252 0.48902503
# [4] 0.73080130 0.80515002 1.07252565
# [7] 1.41290322 1.29394528 1.75940452
save(diff_time, file = "MMLSEM_step2_runtime.RData")
save(MMLSEM_s2output, file = "MixMLSEM_s2output_4non.RData")

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
# [,1]      [,2] [,3]     [,4]
# [1,]    1 -239158.6  183       NA
# [2,]    2 -238654.6  187 2.956831
# [3,]    3 -238484.2  191 1.452951
# [4,]    4 -238366.8  195 1.948478
# [5,]    5 -238306.6  199 1.343079
# [6,]    6 -238268.0  203       NA
# [7,]    7 -238217.0  207 1.581892
# [8,]    8 -238188.6  211 1.019407
# [9,]    9 -238160.8  215       NA

clustering <- t(apply(MMLSEM_s2output[[4]][["posteriors"]], 1, function(x) as.numeric(x == max(x))))
table(data_subset$CNTRYID)
# 7    12    18    21    24    25    28    32    33    34    40    42    43    48    50    55 
# 4475  4402  4434  4808  3889 10934  4954  9338  3491  3599  3343  2497  4502  5542  4454  5920 
# 57    92    93    95 
# 9093  7741  5974  2304 

cbind(clustering, as.numeric(names(table(data_subset$CNTRYID))))

> abs(bnew[,,4])-abs(bold[,,1])
att         mot adapt
att    0.000000000 0.000000000     0
mot    0.014770824 0.000000000     0
adapt -0.002143149 0.002993038     0
> abs(bnew[,,1])-abs(bold[,,3])
att         mot adapt
att    0.000000000 0.000000000     0
mot   -0.004231479 0.000000000     0
adapt -0.002449408 0.004322617     0
> abs(bnew[,,2])-abs(bold[,,2])
att         mot adapt
att    0.00000000  0.00000000     0
mot   -0.05220121  0.00000000     0
adapt  0.01059978 -0.01621269     0
> abs(bnew[,,3])-abs(bold[,,4])
att        mot adapt
att    0.000000000 0.00000000     0
mot   -0.023136338 0.00000000     0
adapt  0.003412966 0.01599169     0