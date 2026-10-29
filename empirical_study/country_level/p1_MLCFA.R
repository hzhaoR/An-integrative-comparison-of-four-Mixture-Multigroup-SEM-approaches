library(MplusAutomation)
library(dplyr)
# load("C:/Users/u0159126/wp4_empirical_01052026/data_66countries_centered.Rdata")
# set.seed(115)
# ids <- sample(c(1:66),20)
# data_subset <- data[data$CNTRYID %in% ids, ]
data <- as.data.frame(data_subset[,c(1,4:24)])

############################with markers#######################################
start_time <- Sys.time()
mod22 <- mplusObject(
  TITLE =  "Multilevel test;",
  VARIABLE = "USEVARIABLES = CNTRYID ST182Q03HA-ST182Q06HA;
WITHIN = ST182Q03HA-ST182Q06HA; 
CLUSTER = CNTRYID;",
  
  DEFINE = "CENTER ST182Q03HA-ST182Q06HA (GROUPMEAN);",
  
  ANALYSIS = "TYPE = RANDOM TWOLEVEL;
ESTIMATOR = BAYES;
ITERATIONS = 1000000;
CONVERGENCE = 0.000001;
",
  
  MODEL ="
%within%
mot BY ST182Q04HA;
s4 | mot BY ST182Q03HA;
s6 | mot BY ST182Q05HA;
s7 | mot BY ST182Q06HA;
logmot |mot;
logST182Q03HA |ST182Q03HA;
logST182Q04HA |ST182Q04HA;
logST182Q05HA |ST182Q05HA;
logST182Q06HA |ST182Q06HA;
[ST182Q03HA-ST182Q06HA@0]
",
  
  SAVEDATA = paste0("\nFILE IS ", "mod2_marker_", "Fscores.dat", ";\nSAVE = FSCORES(100 10); \nFORMAT IS FREE;"),
  
  OUTPUT = "FSCOMPARISON; TECH1; TECH4; TECH8;",
  rdata = as.data.frame(data[,c(1,5:8)]))

res22 <- mplusModeler(mod22, modelout = paste0("mod2_marker_", "factors.inp"), run = 1L)
end_time <- Sys.time()
diff_time <- difftime(end_time, start_time, units='mins') 
# Time difference of 3.574178 mins
####att####
start_time <- Sys.time()
mod11 <- mplusObject(
  TITLE =  "Multilevel test;",
  VARIABLE = "USEVARIABLES = CNTRYID ST181Q02HA-ST181Q04HA;
WITHIN = ST181Q02HA-ST181Q04HA; 
CLUSTER = CNTRYID;",
  
  DEFINE = "CENTER ST181Q02HA-ST181Q04HA (GROUPMEAN);",
  
  ANALYSIS = "TYPE = RANDOM TWOLEVEL;
ESTIMATOR = BAYES;
ITERATIONS = 1000000;
CONVERGENCE = 0.000001;
",
  
  MODEL ="
%within%
att BY ST181Q04HA;
s1 | att BY ST181Q02HA;
s2 | att BY ST181Q03HA;
logatt |att;
logST181Q02HA |ST181Q02HA;
logST181Q03HA |ST181Q03HA;
logST181Q04HA |ST181Q04HA;
[ST181Q02HA-ST181Q04HA@0]

",
  
  SAVEDATA = paste0("\nFILE IS ", "mod1_marker_", "Fscores.dat", ";\nSAVE = FSCORES(100 10); \nFORMAT IS FREE;"),
  
  OUTPUT = "FSCOMPARISON; TECH1; TECH4; TECH8;",
  rdata = as.data.frame(data[,c(1,2:4)]))

res11 <- mplusModeler(mod11, modelout = paste0("mod1_marker_", "factors.inp"), run = 1L)
end_time <- Sys.time()
diff_time <- difftime(end_time, start_time, units='mins') 
# Time difference of 4.019263 mins

####adapt####
start_time <- Sys.time()
mod33 <- mplusObject(
  TITLE =  "Multilevel test;",
  VARIABLE = "USEVARIABLES = CNTRYID ST216Q01HA-ST216Q06HA;
WITHIN = ST216Q01HA-ST216Q06HA; 
CLUSTER = CNTRYID;",
  
  DEFINE = "CENTER ST216Q01HA-ST216Q06HA (GROUPMEAN);",
  
  ANALYSIS = "TYPE = RANDOM TWOLEVEL;
ESTIMATOR = BAYES;
ITERATIONS = 1000000;
CONVERGENCE = 0.000001;
",
  
  MODEL ="
%within%
adapt BY ST216Q03HA;
s8 | adapt BY ST216Q01HA;
s9 | adapt BY ST216Q02HA;
s11 | adapt BY ST216Q04HA;
s12 | adapt BY ST216Q05HA;
s13 | adapt BY ST216Q06HA;
logadapt |adapt;
logST216Q01HA | ST216Q01HA;
logST216Q02HA | ST216Q02HA;
logST216Q03HA | ST216Q03HA;
logST216Q04HA | ST216Q04HA;
logST216Q05HA | ST216Q05HA;
logST216Q06HA | ST216Q06HA;
[ST216Q01HA-ST216Q06HA@0]
",
  
  SAVEDATA = paste0("\nFILE IS ", "mod33_marker_", "Fscores.dat", ";\nSAVE = FSCORES(100 10); \nFORMAT IS FREE;"),
  
  OUTPUT = "FSCOMPARISON; TECH1; TECH4; TECH8;",
  rdata = as.data.frame(data[,c(1,9:14)]))

res33 <- mplusModeler(mod33, modelout = paste0("mod3_marker_", "factors.inp"), run = 1L)
end_time <- Sys.time()
diff_time <- difftime(end_time, start_time, units='mins')


##########################no markers all random##############################
####mot####
start_time <- Sys.time()
mod2 <- mplusObject(
  TITLE =  "Multilevel test;",
  VARIABLE = "USEVARIABLES = CNTRYID ST182Q03HA-ST182Q06HA;
WITHIN = ST182Q03HA-ST182Q06HA; 
CLUSTER = CNTRYID;",
  
  DEFINE = "CENTER ST182Q03HA-ST182Q06HA (GROUPMEAN);",
  
  ANALYSIS = "TYPE = RANDOM TWOLEVEL;
ESTIMATOR = BAYES;
ITERATIONS = 1000000;
CONVERGENCE = 0.000001;
",
  
  MODEL ="
%within%
s4 | mot BY ST182Q03HA;
s5 | mot BY ST182Q04HA;
s6 | mot BY ST182Q05HA;
s7 | mot BY ST182Q06HA;

logmot |mot;

logST182Q03HA |ST182Q03HA;
logST182Q04HA |ST182Q04HA;
logST182Q05HA |ST182Q05HA;
logST182Q06HA |ST182Q06HA;

[ST182Q03HA-ST182Q06HA@0]
",
  
  SAVEDATA = paste0("\nFILE IS ", "mod2", "Fscores.dat", ";\nSAVE = FSCORES(100 10); \nFORMAT IS FREE;"),
  
  OUTPUT = "FSCOMPARISON; TECH1; TECH4; TECH8;",
  rdata = as.data.frame(data[,c(1,5:8)]))

res2 <- mplusModeler(mod2, modelout = paste0("mod2", "factors.inp"), run = 1L)
end_time <- Sys.time()
diff_time_mlcfa1 <- difftime(end.time, start.time, units='mins') 

####att####
start_time <- Sys.time()
mod1 <- mplusObject(
  TITLE =  "Multilevel test;",
  VARIABLE = "USEVARIABLES = CNTRYID ST181Q02HA-ST181Q04HA;
WITHIN = ST181Q02HA-ST181Q04HA; 
CLUSTER = CNTRYID;",
  
  DEFINE = "CENTER ST181Q02HA-ST181Q04HA (GROUPMEAN);",
  
  ANALYSIS = "TYPE = RANDOM TWOLEVEL;
ESTIMATOR = BAYES;
ITERATIONS = 1000000;
CONVERGENCE = 0.000001;
",
  
  MODEL ="
%within%
s1 | att BY ST181Q02HA;
s2 | att BY ST181Q03HA;
s3 | att BY ST181Q04HA;

logatt |att;

logST181Q02HA |ST181Q02HA;
logST181Q03HA |ST181Q03HA;
logST181Q04HA |ST181Q04HA;
[ST181Q02HA-ST181Q04HA@0]
",
  
  SAVEDATA = paste0("\nFILE IS ", "mod1", "Fscores.dat", ";\nSAVE = FSCORES(100 10); \nFORMAT IS FREE;"),
  
  OUTPUT = "FSCOMPARISON; TECH1; TECH4; TECH8;",
  rdata = as.data.frame(data[,c(1,2:4)]))

res1 <- mplusModeler(mod1, modelout = paste0("mod1", "factors.inp"), run = 1L)
end_time <- Sys.time()
diff_time_mlcfa2 <- difftime(end.time, start.time, units='mins') 

####adapt####
start_time <- Sys.time()
mod3 <- mplusObject(
  TITLE =  "Multilevel test;",
  VARIABLE = "USEVARIABLES = CNTRYID ST216Q01HA-ST216Q06HA;
WITHIN = ST216Q01HA-ST216Q06HA; 
CLUSTER = CNTRYID;",
  
  DEFINE = "CENTER ST216Q01HA-ST216Q06HA (GROUPMEAN);",
  
  ANALYSIS = "TYPE = RANDOM TWOLEVEL;
ESTIMATOR = BAYES;
ITERATIONS = 1000000;
CONVERGENCE = 0.000001;
",
  
  MODEL ="
%within%
s8 | adapt BY ST216Q01HA;
s9 | adapt BY ST216Q02HA;
s10 | adapt BY ST216Q03HA;
s11 | adapt BY ST216Q04HA;
s12 | adapt BY ST216Q05HA;
s13 | adapt BY ST216Q06HA;

logadapt |adapt;

logST216Q01HA | ST216Q01HA;
logST216Q02HA | ST216Q02HA;
logST216Q03HA | ST216Q03HA;
logST216Q04HA | ST216Q04HA;
logST216Q05HA | ST216Q05HA;
logST216Q06HA | ST216Q06HA;

[ST216Q01HA-ST216Q06HA@0]

",
  
  SAVEDATA = paste0("\nFILE IS ", "mod3", "Fscores.dat", ";\nSAVE = FSCORES(100 10); \nFORMAT IS FREE;"),
  
  OUTPUT = "FSCOMPARISON; TECH1; TECH4; TECH8;",
  rdata = as.data.frame(data[,c(1,9:14)]))

res3 <- mplusModeler(mod3, modelout = paste0("mod3", "factors.inp"), run = 1L)
end_time <- Sys.time()
diff_time_mlcfa3 <- difftime(end.time, start.time, units='mins')


