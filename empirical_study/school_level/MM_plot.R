##check MM clusters 1-5
# centered data
library(lavaan)
group = "CNTSCHID"
g_name <- as.character(unique(data[, group]))
vars <- lavNames(lavaanify(mm_model, auto = TRUE)) #observed var
lat_var <- lavNames(lavaanify(mm_model, auto = TRUE), "lv") #latent var

centered <- data
group.idx <- match(data[,group], g_name)
group.sizes <- tabulate(group.idx)
group.means <- rowsum.default(as.matrix(data[,vars]),
                              group = group.idx, reorder = FALSE,
                              na.rm = FALSE)/group.sizes
centered[,vars] <- data[,vars] - group.means[group.idx, ,drop = FALSE]

N_gs <- group.sizes #as.numeric(table(data[,group]))
nfactors <- length(lat_var) #ncol(MM.design) #4 lv
ngroups <- length(N_gs)

set.seed(100)
output1 <-mixmgfa(data = data[,3:16], N_gs = N_gs, nfactors = nfactors, 
                  cluster.spec = "loadings", nsclust = c(1:5), 
                  maxiter = 100000, nruns = 100, design = loading_matrix, invar_loadings = c(1,1,1,
                                                                                             0,1,0,1,
                                                                                             1,1,1,0,0,0))
save(output1,file = "mmgfa_output_MMplot.RData")
plot(output1$overview, AIC=1)
