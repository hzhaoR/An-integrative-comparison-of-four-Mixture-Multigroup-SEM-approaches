##0315##
library(MplusAutomation)
data <- as.data.frame(data_subset[,c(1,4:24)])

# data$CNTRYID_new <- as.integer(factor(data$CNTRYID))
####att####
start_time <- Sys.time()
mod1_ami <- mplusObject(
  TITLE =  "BSEM APPROXIMATE MI;",
  VARIABLE = "USEVARIABLES = CNTRYID ST181Q02HA-ST181Q04HA;
CLASSES = c(20); 
KNOWNCLASS = c(CNTRYID= 7 62 43 65 32 48 12 25 18 21 50 34 55 33 40 28 42 57 63 24);
",
  
  ANALYSIS = "TYPE = MIXTURE;
ESTIMATOR = BAYES;
PROCESSORS = 2;
BITERATIONS = 10000;
MODEL=ALLFREE;
",
  
  MODEL ="
%OVERALL%
att BY ST181Q04HA ST181Q02HA ST181Q03HA (lam#_1-lam#_3);

ST181Q02HA-ST181Q04HA;

[ST181Q02HA-ST181Q04HA@0];

",
  MODELPRIORS = "DIFF (lam2_3-lam20_3)~N(0,0.1);",

  
  SAVEDATA = paste0("\nFILE IS ", "mod1_ami", "Fscores1.dat", ";\nSAVE = FSCORES(100 10); \nFORMAT IS FREE;"),
  
  OUTPUT = "TECH1; TECH8; TECH9;",
  rdata = as.data.frame(data[,c(1,2:4)]))

res1_ami <- mplusModeler(mod1_ami, modelout = paste0("ml_mod1_ami.inp"), run = 1L)
end_time <- Sys.time()
diff_time <- difftime(end_time, start_time, units='mins') 
# Time difference of 3.43448 mins

# *** ERROR in MODEL PRIORS command
# The DO function requires 2 or 3 arguments for the range and index string.
# And invalid number of arguments found for the DO function.
# Problem with:  DO(3)  DIFF(LAM2_#-LAM20_#) ~N(0,0.1) ;
#                            
                           

####mot####
start_time <- Sys.time()
mod2_ami <- mplusObject(
  TITLE =  "BSEM APPROXIMATE MI;",
  VARIABLE = "USEVARIABLES = CNTRYID ST182Q03HA-ST182Q06HA;
CLASSES = c(20); 
KNOWNCLASS = c(CNTRYID = 7 62 43 65 32 48 12 25 18 21 50 34 55 33 40 28 42 57 63 24);
",
  
  ANALYSIS = "TYPE = MIXTURE;
ESTIMATOR = BAYES;
PROCESSORS = 2;
BITERATIONS = 10000;
MODEL=ALLFREE;
",
  
  MODEL ="
%OVERALL%
mot BY ST182Q04HA ST182Q03HA ST182Q05HA ST182Q06HA (lam#_1-lam#_4);

ST182Q03HA-ST182Q06HA;

[ST182Q03HA-ST182Q06HA@0];

",
  MODELPRIORS = "DIFF (lam2_2-lam20_2)~N(0,0.1);",
  # DO(2,6) DIFF (res2_#-res20_#)~N(0, 0.1);
  
  SAVEDATA = paste0("\nFILE IS ", "mod2_ami", "Fscores1.dat", ";\nSAVE = FSCORES(100 10); \nFORMAT IS FREE;"),
  
  OUTPUT = "TECH1; TECH8; TECH9;",
  rdata = as.data.frame(data[,c(1,5:8)]))

res2_ami <- mplusModeler(mod2_ami, modelout = paste0("ml_mod2_ami.inp"), run = 1L)
end_time <- Sys.time()
diff_time <- difftime(end_time, start_time, units='mins') 
# Time difference of 3.111779 mins

####adapt####
start_time <- Sys.time()
mod3_ami <- mplusObject(
  TITLE =  "BSEM APPROXIMATE MI;",
  VARIABLE = "USEVARIABLES =  CNTRYID ST216Q01HA-ST216Q06HA;
CLASSES = c(20); 
KNOWNCLASS = c(CNTRYID= 7 62 43 65 32 48 12 25 18 21 50 34 55 33 40 28 42 57 63 24);
",
  
  ANALYSIS = "TYPE = MIXTURE;
ESTIMATOR = BAYES;
PROCESSORS = 2;
BITERATIONS = 10000;
MODEL=ALLFREE;
",
  
  MODEL ="
%OVERALL%
adapt BY ST216Q03HA ST216Q01HA ST216Q02HA ST216Q04HA ST216Q06HA ST216Q05HA(lam#_1-lam#_6);

ST216Q01HA-ST216Q06HA;

[ST216Q01HA-ST216Q06HA@0];

",
  MODELPRIORS = "DO (2,5) DIFF (lam2_#-lam20_#)~N(0,0.1);",
  # DO(2,6) DIFF (res2_#-res20_#)~N(0, 0.1);
  
  SAVEDATA = paste0("\nFILE IS ", "mod3_ami", "Fscores1.dat", ";\nSAVE = FSCORES(100 10); \nFORMAT IS FREE;"),
  
  OUTPUT = "TECH1; TECH8; TECH9;",
  rdata = as.data.frame(data[,c(1,9:14)]))

res3_ami <- mplusModeler(mod3_ami, modelout = paste0("ml_mod3_ami.inp"), run = 1L)
end_time <- Sys.time()
diff_time <- difftime(end_time, start_time, units='mins') 
# Time difference of 4.066734 mins

# ####compenv####
# start_time <- Sys.time()
# mod4_ami <- mplusObject(
#   TITLE =  "BSEM APPROXIMATE MI;",
#   VARIABLE = "USEVARIABLES =  CNTRYID ST205Q01HA-ST205Q04HA;
# CLASSES = c(20); 
# KNOWNCLASS = c(CNTRYID= 7 62 43 65 32 48 12 25 18 21 50 34 55 33 40 28 42 57 63 24);
# ",
#   
#   ANALYSIS = "TYPE = MIXTURE;
# ESTIMATOR = BAYES;
# PROCESSORS = 2;
# BITERATIONS = 10000;
# MODEL=ALLFREE;
# ",
#   
#   MODEL ="
# %OVERALL%
# compenv BY ST205Q02HA ST205Q01HA ST205Q03HA ST205Q04HA (lam#_1-lam#_4);
# 
# ST205Q01HA-ST205Q04HA;
# 
# [ST205Q01HA-ST205Q04HA@0];
# 
# ",
#   MODELPRIORS = "DIFF (lam2_3-lam20_3)~N(0,0.1);",
#   # DO(2,6) DIFF (res2_#-res20_#)~N(0, 0.1);
#   
#   SAVEDATA = paste0("\nFILE IS ", "mod4_ami", "Fscores1.dat", ";\nSAVE = FSCORES(100 10); \nFORMAT IS FREE;"),
#   
#   OUTPUT = "TECH1; TECH8; TECH9;",
#   rdata = as.data.frame(data[,c(1,15:18)]))
# 
# res4_ami <- mplusModeler(mod4_ami, modelout = paste0("ml_mod4_ami.inp"), run = 1L)
# end_time <- Sys.time()
# 
# ####collenv####
# start_time <- Sys.time()
# mod5_ami <- mplusObject(
#   TITLE =  "BSEM APPROXIMATE MI;",
#   VARIABLE = "USEVARIABLES =  CNTRYID ST206Q01HA-ST206Q04HA;
# CLASSES = c(20); 
# KNOWNCLASS = c(CNTRYID= 7 62 43 65 32 48 12 25 18 21 50 34 55 33 40 28 42 57 63 24);
# ",
#   
#   ANALYSIS = "TYPE = MIXTURE;
# ESTIMATOR = BAYES;
# PROCESSORS = 2;
# BITERATIONS = 10000;
# MODEL=ALLFREE;
# ",
#   
#   MODEL ="
# %OVERALL%
# collenv BY ST206Q03HA ST206Q01HA ST206Q02HA ST206Q04HA (lam#_1-lam#_4);
# 
# ST206Q01HA-ST206Q04HA;
# 
# [ST206Q01HA-ST206Q04HA@0];
# 
# ",
#   MODELPRIORS = "DO(2,4) DIFF (lam2_#-lam20_#)~N(0,0.1);",
#   # DO(2,6) DIFF (res2_#-res20_#)~N(0, 0.1);
#   
#   SAVEDATA = paste0("\nFILE IS ", "mod5_ami", "Fscores1.dat", ";\nSAVE = FSCORES(100 10); \nFORMAT IS FREE;"),
#   
#   OUTPUT = "TECH1; TECH8; TECH9;",
#   rdata = as.data.frame(data[,c(1,19:22)]))
# 
# res5_ami <- mplusModeler(mod5_ami, modelout = paste0("ml_mod5_ami.inp"), run = 1L)
# end_time <- Sys.time()

MGBCFA_output_4non_0.1 <- list(res1_ami = res1_ami,
                res2_ami = res2_ami,
                res3_ami = res3_ami)
save(MGBCFA_output_4non_0.1, file = "MGBCFA_output_4non_ami_0.1_3F.RData")

# 
# library(dplyr)
# #factor1_SEHS
# parameter1 <- res1[["results"]][["parameters"]][["unstandardized"]]
# #factor scores
# fscores1 <- res3_ami[["results"]][["savedata"]]
# #arrange by country
# fscores1 <- fscores1 %>%
#   arrange(COUNTRY)
