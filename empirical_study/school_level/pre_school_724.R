library(haven)
pisa2018 <- read_sav("~/WP4_2025/some prep_July2025/Baysu 2022_mediation models/pisa2018data/STU/CY07_MSU_STU_QQQ.sav")
labels <- sapply(pisa2018, function(x) attr(x, "label"))

data <- pisa2018[,c("CNTRYID", "CNT", "CNTSCHID",
                    "ST181Q02HA","ST181Q03HA", "ST181Q04HA", #Competitive Attitude (Attitude) 4 point scale
                    "ST182Q03HA", "ST182Q04HA", "ST182Q05HA", "ST182Q06HA", #Task Motivation (Motivate) 4 point scale
                    "ST216Q01HA", "ST216Q02HA", "ST216Q03HA", "ST216Q04HA", "ST216Q05HA", "ST216Q06HA", #Cognitive Adaptability (Adapt) 5 point scale
                    "ST205Q01HA", "ST205Q02HA", "ST205Q03HA", "ST205Q04HA", #Competitive Environment (CompEnv) 4 point scale
                    "ST206Q01HA", "ST206Q02HA", "ST206Q03HA", "ST206Q04HA")] #Collaborative Environment (CollEnv) 4 point scale
#612004 obs in total
data <- na.omit(data) #336513 obs
save(data, file = "data_org_noNA.RData")

#country level
max(table(data$CNTRYID)) #66 countries
dat0 <- data[data$CNTRYID==724,]
#Belgium: CNTRYID==56; Q216 adaptability items NA;
#Austria: 40; 2221 obs/79 groups;
#France: 250; 1334/49 groups;
#Germany: 276; >25:0
#Netherlands: 528; Q216 adaptability items NA;
#Switzerland: 756; 344/13 groups;

# dat0 <- data[data$CNTRYID==724,] #23452 obs
max(table(data$CNTSCHID))

library(dplyr)
data_filtered <- dat0 %>%
  group_by(CNTSCHID) %>%
  filter(n() >= 25) %>%
  ungroup()
length(table(data_filtered$CNTSCHID))
#school=40/250/756: 3899 obs left; 141 groups

#school=8/31/32: 5774 obs left; 188 groups
#school=724: 14082 obs left; 462 groups
data <- data_filtered
save(data, file = "data_724_centered_3F.Rdata")

#pre-centered data
g_name <- as.character(unique(data$CNTSCHID))
data$CNTSCHID <- as.character(data$CNTSCHID)
group.idx <- match(data$CNTSCHID, g_name)
group.sizes <- tabulate(group.idx)
vars <- colnames(data)
vars <- vars[4:24]
group.means <- rowsum.default(as.matrix(data[,vars]),
                              group = group.idx, reorder = FALSE,
                              na.rm = FALSE)/group.sizes
data[,vars] <- data[,vars] - group.means[group.idx, , drop = FALSE] 
data$CNTSCHID <- as.numeric(group.idx)
# data <- data_filtered

mm_model <- 
  # measurement
  'att =~ ST181Q02HA + ST181Q03HA + ST181Q04HA 
   mot  =~ ST182Q03HA + ST182Q04HA + ST182Q05HA + ST182Q06HA
   adapt =~ ST216Q01HA + ST216Q02HA + ST216Q03HA + ST216Q04HA + ST216Q05HA + ST216Q06HA'
   # compenv =~ ST205Q01HA + ST205Q02HA + ST205Q03HA + ST205Q04HA
   # collenv =~ ST206Q01HA + ST206Q02HA + ST206Q03HA + ST206Q04HA'

# # structural 
sm_model <- 
'mot ~ att
 adapt ~ att + mot'

# sm_model <- '
#   # regressions
#   att ~ compenv + collenv
#   mot ~ att + compenv + collenv
#   adapt    ~ att + mot
# 
#   # covariance between exogenous latent variables
#   compenv ~~ collenv
# '

#2mixmgsem
##check MM clusters 1-6
##MMGFA

#CFA loading matrix
# MM.design <- matrix(0, ncol = nfactors, nrow = nobsvar)
# rownames(MM.design) <- obs_var
# colnames(MM.design) <- lat_var
# MM.design[1:5,1] <- MM.design[6:10,2] <- MM.design[11:15,3] <- MM.design[16:20,4] <- 1
create_loading_matrix <- function(model_string) {
  # Split the model string into lines and filter out the factor loading lines
  factor_lines <- grep("=~", strsplit(model_string, "\n")[[1]], value = TRUE)
  
  # Create a list where each element contains the indicators for a factor
  loadings <- lapply(factor_lines, function(line) {
    split_line <- strsplit(line, "=~")[[1]]
    factor <- trimws(split_line[1])  # The factor name (left side of =~)
    indicators <- trimws(unlist(strsplit(split_line[2], "\\+")))  # Indicators (right side of =~)
    return(list(factor = factor, indicators = indicators))
  })
  
  # Get unique indicators and factors
  all_indicators <- unique(unlist(lapply(loadings, function(x) x$indicators)))
  all_factors <- sapply(loadings, function(x) x$factor)
  
  # Create a loading matrix of zeros
  design <- matrix(0, nrow = length(all_indicators), ncol = length(all_factors),
                   dimnames = list(all_indicators, all_factors))
  
  # Assign 1s to the appropriate factor-indicator relationships
  for (i in seq_along(loadings)) {
    indicators <- loadings[[i]]$indicators
    factor <- loadings[[i]]$factor
    design[indicators, factor] <- 1
  }
  
  return(design)
}
loading_matrix <- create_loading_matrix(mm_model)

#for the simulated data
data <- data[,c(3:24)]
data <- as.data.frame(data)
which(table(data$CNTSCHID)<20)

library(lavaan)
step1model = mm_model
# # centered data
# g_name <- as.character(unique(data[, group]))
# vars <- lavNames(lavaanify(step1model, auto = TRUE)) #observed var
lat_var <- lavNames(lavaanify(step1model, auto = TRUE), "lv") #latent var

# centered <- data
# group.idx <- match(data[,group], g_name)
# group.sizes <- tabulate(group.idx)
# group.means <- rowsum.default(as.matrix(data[,vars]),
#                               group = group.idx, reorder = FALSE,
#                               na.rm = FALSE)/group.sizes
# centered[,vars] <- data[,vars] - group.means[group.idx, ,drop = FALSE]
# 
N_gs <- group.sizes #as.numeric(table(data[,group]))
nfactors <- length(lat_var) #ncol(MM.design) #4 lv
# ngroups <- length(N_gs)
# 
# s1output <- list()
# for (i in 1:6) {
# 
#   s1output[[i]] <- MixMix_Step1(data, step1model = mm_model, group = "CNTSCHID", MM.cluster.spec=c("loadings"),
#                                 MM.nclus = i, MM.maxiter = 100000, MM.nruns = 50, MM.design = loading_matrix,
#                                 invar_loadings = c(1,rep(0,2), 1,rep(0,3), 1,rep(0,5), 1,rep(0,3), 1,rep(0,3)),
#                                 markers = markers, seed = 100)
# }
# 
# save(s1output, file = "s1output_6k_3schools_converged.RData")


library(mixmgfa)
set.seed(100)
# data <- data[,1:14]
output <-mixmgfa(data = data[,1:14], N_gs = N_gs, nfactors = nfactors,
                  cluster.spec = c("loadings"), nsclust = c(4),
                  maxiter = 100000, nruns = 100, design = loading_matrix,
                 invar_loadings = c(0,0,1,
                                    0,1,0,0,
                                    0,0,1,0,0,0))
                  # invar_loadings = c(0,rep(0,2), 0,rep(0,3), 0,rep(0,5)))
plot(output$overview, AIC=1)
> output[["MMGFAsolutions"]][["4.clusters"]][["clusterspecific.loadings"]][["Cluster_1"]]
Factor_1  Factor_2  Factor_3
ST181Q02HA 0.6015870 0.0000000 0.0000000
ST181Q03HA 0.5816163 0.0000000 0.0000000
ST181Q04HA 0.6631564 0.0000000 0.0000000
ST182Q03HA 0.0000000 0.4140858 0.0000000
ST182Q04HA 0.0000000 0.4843620 0.0000000
ST182Q05HA 0.0000000 0.4538482 0.0000000
ST182Q06HA 0.0000000 0.4524252 0.0000000
ST216Q01HA 0.0000000 0.0000000 0.6215795
ST216Q02HA 0.0000000 0.0000000 0.6227783
ST216Q03HA 0.0000000 0.0000000 0.7245926
ST216Q04HA 0.0000000 0.0000000 0.6548707
ST216Q05HA 0.0000000 0.0000000 0.6309455
ST216Q06HA 0.0000000 0.0000000 0.6426489
> output[["MMGFAsolutions"]][["4.clusters"]][["clusterspecific.loadings"]][["Cluster_2"]]
Factor_1  Factor_2  Factor_3
ST181Q02HA 0.4989102 0.0000000 0.0000000
ST181Q03HA 0.5593684 0.0000000 0.0000000
ST181Q04HA 0.6682672 0.0000000 0.0000000
ST182Q03HA 0.0000000 0.2013049 0.0000000
ST182Q04HA 0.0000000 0.4772611 0.0000000
ST182Q05HA 0.0000000 0.4035556 0.0000000
ST182Q06HA 0.0000000 0.5520416 0.0000000
ST216Q01HA 0.0000000 0.0000000 0.6742474
ST216Q02HA 0.0000000 0.0000000 0.6136325
ST216Q03HA 0.0000000 0.0000000 0.7333597
ST216Q04HA 0.0000000 0.0000000 0.3834277
ST216Q05HA 0.0000000 0.0000000 0.6198588
ST216Q06HA 0.0000000 0.0000000 0.4748668
> output[["MMGFAsolutions"]][["4.clusters"]][["clusterspecific.loadings"]][["Cluster_3"]]
Factor_1  Factor_2  Factor_3
ST181Q02HA 0.6328321 0.0000000 0.0000000
ST181Q03HA 0.5215936 0.0000000 0.0000000
ST181Q04HA 0.6678024 0.0000000 0.0000000
ST182Q03HA 0.0000000 0.5197698 0.0000000
ST182Q04HA 0.0000000 0.3816109 0.0000000
ST182Q05HA 0.0000000 0.5399326 0.0000000
ST182Q06HA 0.0000000 0.3588396 0.0000000
ST216Q01HA 0.0000000 0.0000000 0.6738172
ST216Q02HA 0.0000000 0.0000000 0.6939233
ST216Q03HA 0.0000000 0.0000000 0.8080415
ST216Q04HA 0.0000000 0.0000000 0.4400984
ST216Q05HA 0.0000000 0.0000000 0.5345495
ST216Q06HA 0.0000000 0.0000000 0.4687560
> output[["MMGFAsolutions"]][["4.clusters"]][["clusterspecific.loadings"]][["Cluster_4"]]
Factor_1  Factor_2  Factor_3
ST181Q02HA 0.5371490 0.0000000 0.0000000
ST181Q03HA 0.6563688 0.0000000 0.0000000
ST181Q04HA 0.6879834 0.0000000 0.0000000
ST182Q03HA 0.0000000 0.5008378 0.0000000
ST182Q04HA 0.0000000 0.4081221 0.0000000
ST182Q05HA 0.0000000 0.5375571 0.0000000
ST182Q06HA 0.0000000 0.4408261 0.0000000
ST216Q01HA 0.0000000 0.0000000 0.3361139
ST216Q02HA 0.0000000 0.0000000 0.4826206
ST216Q03HA 0.0000000 0.0000000 0.5200258
ST216Q04HA 0.0000000 0.0000000 0.7773904
ST216Q05HA 0.0000000 0.0000000 0.7578751
ST216Q06HA 0.0000000 0.0000000 0.8575530

# output1$overview
nr of clusters    loglik nrpars
[1,]              1 -196661.4  14792
[2,]              2 -196582.4  14803
[3,]              3 -196513.9  14814
[4,]              4 -196474.7  14826
[5,]              5 -196451.7  14837
BIC_N    BIC_G      AIC
[1,] 534625.7 484080.1 422906.8
[2,] 534572.7 483989.5 422770.8
[3,] 534540.8 483920.0 422655.8
[4,] 534577.0 483915.3 422601.4
[5,] 534636.2 483936.8 422577.5
screeratios convergence
[1,]          NA           1
[2,]    1.153927           1
[3,]    1.906682           1
[4,]    1.564448           1
[5,]          NA           1

output2 <- ScaleRotateMixmgfa_pinvar(output, N_gs = N_gs,
                                     nsclust = 4, design = loading_matrix, rescale=1, markers = markers,
                                     rotation=0,targetT=0,targetW=0)


# #MG-SEM
# library(lavaan)
# mgcfa_model <- '
#   # measurement
#   att =~ ST181Q04HA + ST181Q02HA + ST181Q03HA  
#    mot  =~ ST182Q03HA + ST182Q04HA + ST182Q05HA + ST182Q06HA
#    adapt =~ ST216Q03HA + ST216Q01HA + ST216Q02HA + ST216Q04HA + ST216Q05HA + ST216Q06HA
#    compenv =~ ST205Q01HA + ST205Q02HA + ST205Q03HA + ST205Q04HA
#    collenv =~ ST206Q01HA + ST206Q02HA + ST206Q03HA + ST206Q04HA
# '
# 
# # Fit multigroup SEM treating items as continuous, robust ML
# fit1_mg <- sem(mgcfa_model,
#                data = data,
#                group = "CNTSCHID",
#                estimator = "MLR",    # robust ML for continuous ordinal-ish data
#                std.lv = TRUE)        # identify by fixing latent variances = 1
# Error: cannot allocate vector of size 808.3 Mb
# In addition: Warning messages:
#   1: In doTryCatch(return(expr), name, parentenv, handler) :
#   display list redraw incomplete
# 2: In doTryCatch(return(expr), name, parentenv, handler) :
#   invalid graphics state
# 3: In doTryCatch(return(expr), name, parentenv, handler) :
#   invalid graphics state

# # Summary: fit measures + standardized estimates
# summary(fit1_mg, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
# round(fitmeasures(object = fit1_mg)[c("cfi", "rmsea")], 3)  
# # cfi rmsea 
# # 0.960 0.046 
# 
# # Get parameter table and standardized solution
# # pe <- parameterEstimates(fit1_mg) # unstandardized estimates, group column included
# # est <- lavInspect(fit1_mg, "est")
# # load <- lapply(est, "[[", "lambda") 
# 
# # metric invariance
# fit2_mg <- sem(mgcfa_model,
#                data = data[1:117862,],
#                group = "CNTRYID",
#                estimator = "MLR",    # robust ML for continuous ordinal-ish data
#                std.lv = TRUE, 
#                group.equal = "loadings")        # identify by fixing latent variances = 1
# round(fitmeasures(object = fit2_mg)[c("cfi", "rmsea")], 3)  
# # cfi rmsea 
# # 0.956 0.047 
# est2 <- lavInspect(fit2_mg, "est")
# save(est2, file = "est_fit2.RData")
# 
