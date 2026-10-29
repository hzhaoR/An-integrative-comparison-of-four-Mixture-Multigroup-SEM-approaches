library(mixmgfa)
library(lavaan)

ScaleRotateMixmgfa_pinvar <- function(OutputObject,N_gs,nsclust,design=0,rescale=0,markers=0,rotation=0,targetT=0,targetW=0){
  
  if(length(OutputObject)==2){
    Output=OutputObject$MMGFAsolutions
  } else {
    Output=OutputObject
  }# Output should be a list of one in case it includes only one number of clusters
  nsclustOutput=lapply(Output,function(x) ncol(x$clustermemberships)) 
  
  # check whether 'cluster.spec' is correctly specified (whether 'residuals' is cluster-specific or not doesn't matter for this function)
  # if(is.element("loadings",cluster.spec) && is.null(Output[[1]]$clusterspecific.loadings)){
  #   stop("The specification of cluster.spec does not match the OutputObject. Please specify cluster.spec correctly and try again.")
  # }
  # if(is.element("intercepts",cluster.spec) && is.null(Output[[1]]$clusterspecific.intercepts)){
  #   stop("The specification of cluster.spec does not match the OutputObject. Please specify cluster.spec correctly and try again.")
  # }
  # derive 'cluster.spec' from input (whether 'residuals' is cluster-specific or not doesn't matter for this function)
  cluster.spec=c()
  if(is.null(Output[[1]]$clusterspecific.loadings)==FALSE){
    cluster.spec=c(cluster.spec,"loadings")
  }
  if(is.null(Output[[1]]$clusterspecific.intercepts)==FALSE){
    cluster.spec=c(cluster.spec,"intercepts")
  }
  
  if(is.element("loadings",cluster.spec)){
    nfactors=ncol(Output[[1]]$clusterspecific.loadings[[1]])
    nvar=nrow(Output[[1]]$clusterspecific.loadings[[1]])
  } else {
    nfactors=ncol(Output[[1]]$invariant.loadings)
    nvar=nrow(Output[[1]]$invariant.loadings)
  }
  # browser( )
  if(is.element("loadings",cluster.spec)){ # check for invariant loadings across clusters
    ind_nclust=which(nsclustOutput==nsclust[1])
    loadingsmatrix<-matrix(unlist(Output[[ind_nclust]]$clusterspecific.loadings), ncol = nfactors*nsclust[1], nrow=nrow(Output[[ind_nclust]]$clusterspecific.loadings[[ind_nclust]]), byrow = FALSE)
    equalloadings=matrix(0,nvar,nfactors)
    for(q in 1:nfactors){
      indcols=seq(from=1, to=ncol(loadingsmatrix),by=nfactors)
      uniqueloads=apply(matrix(loadingsmatrix[,indcols], nrow = nvar),1,unique)
      equalloadings[which(lapply(uniqueloads, function(x) length(x)==1)==TRUE),q]=1
    }
    invar.loadings=matrix(0,1,nvar)
    invar.loadings[which(apply(equalloadings,1,sum)==nfactors)]=1
    # in case of some invariant loadings over clusters, rescaling and rotation needs to be done over all groups rather than per cluster, and markers should take this into account
  } else {
    invar.loadings=matrix(1,1,nvar)
  }
  
  if(max(design)!=0){
    if(ncol(design)!=nfactors){
      stop("The specified design matrix does not match the number of factors.")
    }
  }
  if(sum(design)==0){
    EFA=1
  } else {
    EFA=0
  }
  if(max(design)!=0){
    zl=1-design
    if(is.element("loadings",cluster.spec)){
      loadings1=(Output[[1]]$clusterspecific.loadings[[1]])
    } else {
      loadings1=(Output[[1]]$invariant.loadings)
    }
    zl=loadings1*zl # loadings that should be zero according to design matrix
    if(max(abs(zl))>0){
      stop("The specified design matrix does not match the inputted loadings (the zero loadings are not zero). Check whether the columns of design should be in a different order.")
    }
  }
  
  if(rotation==1 && (nfactors==1 || EFA==0)){
    cat("You requested rotation for a one-factor or CFA solution. Rotation will not be performed.")
    cat("\n")
    rotation=0
  }
  
  if(rotation!=0){
    if(rotation=="varimax" || rotation=="Varimax" || rotation=="VARIMAX"){
      rotation="varimax"
    } else if(rotation=="oblimin" || rotation=="Oblimin" || rotation=="OBLIMIN"){
      rotation="oblimin"
    } else if(rotation=="geomin" || rotation=="Geomin" || rotation=="GEOMIN"){
      rotation="geomin"
    } else if(rotation=="target" || rotation=="Target" || rotation=="TARGET"){
      rotation="target"
    } else {
      cat("You seem to have misspelled your rotation criterion. Oblimin rotation will be performed.")
      cat("\n")
      rotation="oblimin"
    }
  }
  if(rotation!=0){
    if(rotation=="target"){
      if(max(targetT)==0){
        stop("You failed to specify a targetT matrix or it contains only zeros. Please specify a valid targetT matrix")
      } else if(max(targetT)!=0){
        if(ncol(targetT)!=nfactors){
          stop("The specified targetT matrix does not match the specified number of factors.")
        }
      }
      if(max(targetW)==0){
        targetW=1-ceiling(targetT)
        if(sum(targetW)<(nfactors*(nfactors-1))){
          targetW=matrix(1,ncol(targetT),nfactors)
          cat("Note: targetW is unspecified. A valid weight matrix could not be derived from targetT, so fully specified target rotation will be performed.")
          cat("\n")
        } else {
          cat("Note: targetW is unspecified. A weight matrix is derived from targetT so that the zeros in targetT are approximated.")
          cat("\n")
        }
      }
    }
    if(!require(GPArotation)){ # test if package is loaded
      loadtest<-try(library(GPArotation),silent=TRUE) # test if package can be loaded (if it is installed)
      if(class(loadtest)=="try-error"){
        install.packages("GPArotation")
      }
    }
  }
  
  if(rescale==1 && min(colSums(as.matrix(markers)))==0){ 
    if(EFA==0){
      cat("You requested rescaling for a CFA solution, but did not specify marker variables (or not for all factors). Default marker variables will be derived from the design matrix.")
      cat("\n")
      markers=matrix(0,nrow(design),ncol(design))
      ssitems=rowSums(design)==1 # items with simple structure in the design
      if(max(invar.loadings)==1){
        ssitems=which(ssitems==1 & invar.loadings==1) # default markers take partial invariance of loadings into account
      } else {
        ssitems=which(ssitems==1)
      }
      markerind=ssitems[apply(design[ssitems,],2,which.max)]
      markers[markerind,seq_len(nfactors)]=1
      markers=markers*design
    } else {
      cat("You requested rescaling for an EFA solution, but did not specify marker variables (or not for all factors). Rescaling will not be performed.")
      cat("\n")
      rescale=0
    }
  } else if(rescale==1 && (EFA==0 || max(targetT)!=0 || (max(invar.loadings)==1 && is.element("loadings",cluster.spec)))){ # checking whether markers match design or targetT or partial invariance of loadings across clusters
    markeritems=which(rowSums(markers)>0)
    if(EFA==0){
      indhl=apply(design[markeritems,],1,which.max) # index of factor each marker item has non-zero loading on according to 'index'
    } else if (max(targetT)!=0){
      indhl=apply(targetT[markeritems,],1,which.max) # index of factor each marker item has non-zero loading on according to 'index'
    } else {
      indhl=seq_len(nfactors)
    }
    if(length(unique(indhl))<nfactors){
      cat("It could be that you specified items with crossloadings as markers. This is not recommended.")
      cat("\n")
    } else {
      ind1=apply(markers[markeritems,],1,which.max) # index of factor each item has a 1 for in 'markers'
      permm=match(ind1,indhl) # match the two to check whether permutation of 'markers' is needed
      markers=markers[,permm]
    }
    if((is.element("loadings",cluster.spec) && max(invar.loadings)==1 && min(invar.loadings[markeritems])==0) || sum(colSums(markers)>0)<nfactors){
      if(EFA==0){
        cat("Invariant loadings across clusters are detected and you specified markers that do not correspond to one of these invariant loadings, or another misspecification occurred. Default markers will be derived.")
        cat("\n")
        markers=matrix(0,nrow(design),ncol(design))
        ssitems=rowSums(design)==1 # items with simple structure in the design
        if(max(invar.loadings)==1){
          ssitems=which(ssitems==1 & invar.loadings==1)
        } else {
          ssitems=which(ssitems==1)
        }
        markerind=ssitems[apply(design[ssitems,],2,which.max)]
        markers[markerind,seq_len(nfactors)]=1
        markers=markers*design
      } else {
        cat("You requested rescaling for an EFA solution, but specified markers that do not correspond to loadings that are detected to be invariant across clusters, or another misspecification occurred. Rescaling will not be performed.")
        cat("\n")
        rescale=0
      }
    }
  }
  
  if(is.element("loadings",cluster.spec)){
    ngroup=nrow(Output[[1]]$group.and.clusterspecific.factorcovariances)
  } else {
    ngroup=nrow(Output[[1]]$groupspecific.factorcovariances)
  }
  
  if(length(nsclust)>2){ # nsclust needs to be of length two
    cat("You specified more than two numbers of clusters. The minimum and maximum of these numbers will be used.")
    cat("\n")
    nsclust=c(min(nsclust),max(nsclust))
  }
  if(length(nsclust)==2){
    if(nsclust[2]>ngroup){
      cat("Note: the number of clusters cannot exceed the number of groups.")
      cat("\n")
      nsclust[2]=ngroup
    }
  } else if(length(nsclust)==1){
    if(nsclust[1]>ngroup){
      cat("Note: the number of clusters cannot exceed the number of groups.")
      cat("\n")
      nsclust[1]=ngroup
    }
    nsclust=c(nsclust,nsclust)
  }
  # Check whether the specified numbers of clusters match the ones in the Output object
  nclust_max=length(Output[[length(Output)]]$clusterproportions)
  nclust_min=length(Output[[1]]$clusterproportions)
  correction=0
  if(nsclust[1]<nclust_min){
    nsclust[1]=nclust_min
    correction=1
  }
  if(nsclust[2]>nclust_max){
    nsclust[2]=nclust_max
    correction=1
  }
  if(correction==1){
    cat("Note: the specified number(s) of clusters exceed the number(s) of clusters in the Output object. This is corrected.")
    cat("\n")
  }
  
  N=sum(N_gs)
  # start of loop over numbers of clusters
  for(nclust in nsclust[1]:nsclust[2]){
    ind_nclust=which(nsclustOutput==nclust)
    # if rotation is requested in case of EFA, make factors orthogonal
    # per cluster or overall
    if(EFA==1 && rotation!=0){
      if(is.element("loadings",cluster.spec) && max(invar.loadings)==0){ # per cluster
        z_gks=Output[[ind_nclust]]$clustermemberships
        pi_ks=Output[[ind_nclust]]$clusterproportions
        N_gks=sweep(z_gks,N_gs,MARGIN=1,'*',check.margin = FALSE) #N_gks=diag(N_gs[,1])%*%z_gks
        N_ks=colSums(N_gks)
        for(k in 1:nclust){
          if(N_ks[k]>0){
            phi_k=matrix(0,nfactors,nfactors)
            for(g in 1:ngroup){
              phi_gk=Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]]
              phi_k=phi_k+(N_gks[g,k]/N_ks[k])*phi_gk;
            }
            # find matrix square root via eigenvalue decomposition
            ed=eigen(phi_k)
            sqrtFscale=ed$vectors%*%diag(ed$values)^(1/2)%*%solve(ed$vectors)
            invsqrtFscale=solve(sqrtFscale) # This normally has a diagonal of ones, so that it only orthogonalizes, without rescaling
            for(g in 1:ngroup){
              phi_gk=Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]]
              phi_gk=invsqrtFscale%*%phi_gk%*%invsqrtFscale;
              Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]]=((phi_gk+t(phi_gk))/2); # enforce perfect symmetry
            }
            Output[[ind_nclust]]$clusterspecific.loadings[[k]]=Output[[ind_nclust]]$clusterspecific.loadings[[k]]%*%sqrtFscale # compensate for (re)scaling of factors in the loadings
            if(is.element("intercepts",cluster.spec)){
              for(g in 1:ngroup){ # transform the factor means accordingly
                Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]]=Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]]%*%invsqrtFscale
              }
            }
          }
        }
      } else { # overall (in case of (some) invariant loadings)
        phi=matrix(0,nfactors,nfactors)
        if(is.element("loadings",cluster.spec) && max(invar.loadings)==1){ # cluster-specific loadings with some invariances
          z_gks=Output[[ind_nclust]]$clustermemberships
          pi_ks=Output[[ind_nclust]]$clusterproportions
          N_gks=sweep(z_gks,N_gs,MARGIN=1,'*',check.margin = FALSE) #N_gks=diag(N_gs[,1])%*%z_gks
          N_ks=colSums(N_gks)
          for(k in 1:nclust){
            if(N_ks[k]>0){
              for(g in 1:ngroup){
                phi_gk=Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]]
                phi=phi+(N_gks[g,k]/N)*phi_gk
              }
            }
          }
        } else { # invariant loadings
          for(g in 1:ngroup){
            phi_g=Output[[ind_nclust]]$groupspecific.factorcovariances[[g]]
            phi=phi+(N_gs[g]/N)*phi_g
          }
        }
        # find matrix square root via eigenvalue decomposition
        ed=eigen(phi)
        sqrtFscale=ed$vectors%*%diag(ed$values)^(1/2)%*%solve(ed$vectors)
        invsqrtFscale=solve(sqrtFscale)
        if(is.element("loadings",cluster.spec) && max(invar.loadings)==1){ # cluster-specific loadings with some invariances
          for(k in 1:nclust){
            if(N_ks[k]>0){
              for(g in 1:ngroup){
                phi_gk=Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]]
                phi_gk=invsqrtFscale%*%phi_gk%*%invsqrtFscale
                Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]]=((phi_gk+t(phi_gk))/2)
              }
              Output[[ind_nclust]]$clusterspecific.loadings[[k]]=Output[[ind_nclust]]$clusterspecific.loadings[[k]]%*%sqrtFscale # compensate for (re)scaling of factors in the loadings
            }
          }
        } else { # invariant loadings
          for(g in 1:ngroup){
            phi_g=Output[[ind_nclust]]$groupspecific.factorcovariances[[g]]
            phi_g=invsqrtFscale%*%phi_g%*%invsqrtFscale
            Output[[ind_nclust]]$groupspecific.factorcovariances[[g]]=((phi_g+t(phi_g))/2) # enforce perfect symmetry
          }
          Output[[ind_nclust]]$invariant.loadings=Output[[ind_nclust]]$invariant.loadings%*%sqrtFscale # compensate for (re)scaling/orthogonalizing factors in the loadings
        }
        if(is.element("intercepts",cluster.spec)){ # for "(loadings)intercepts" and "(loadings)interceptsresiduals" variants, group-and-cluster-specific factor means need to be adjusted
          for(g in 1:ngroup){ # transform the factor means accordingly
            for(k in 1:nclust){
              Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]]=Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]]%*%invsqrtFscale
            }
          }
        } else {
          if(is.element("loadings",cluster.spec)==0){ # for "residuals" variant, group-specific factor means need to be adjusted
            for(g in 1:ngroup){ # transform the factor means accordingly
              Output[[ind_nclust]]$groupspecific.factormeans[[g]]=Output[[ind_nclust]]$groupspecific.factormeans[[g]]%*%invsqrtFscale
            }
          }
        }
      }
    }
    
    if(is.element("loadings",cluster.spec)){ # rotation/rescaling needs to be performed per cluster
      if(sum(design)==0){ # EFA
        if(max(invar.loadings)>0 && rotation!=0){
          # determine and rotate mean loading matrix
          meanLambda=matrix(0,nvar,nfactors)
          for(k in 1:nclust){
            meanLambda=meanLambda+Output[[ind_nclust]]$clusterspecific.loadings[[k]]
          }
          meanLambda=meanLambda/nclust
          if(rotation=="varimax"){ # rotation based on mean loading matrix
            rot<-GPForth(meanLambda, method="varimax")
          } else if(rotation=="oblimin"){
            rot<-GPFoblq(meanLambda, method="oblimin")
          } else if(rotation=="geomin"){
            rot<-GPFoblq(meanLambda, method="geomin")
          } else if(rotation=="target"){
            rot<-GPFoblq(meanLambda, method="pst",methodArgs =list(W=targetW,Target=targetT))
          }
          # browser()
          meanrotatedloadings=rot$loadings
          strongloadcutoff=mean(apply(abs(meanrotatedloadings),2,max))/2
          strongloadind=which(abs(meanrotatedloadings)>=strongloadcutoff)
          strongload=matrix(0,nrow=nvar,ncol=nfactors)
          strongload[strongloadind]=meanrotatedloadings[strongloadind]
          toreflect=colSums(strongload>0)<colSums(strongload<0)
          reflcommon=matrix(1,nrow=1,ncol=nfactors)
          reflcommon[toreflect]=-1
        }
        for(k in 1:nclust){
          if(rotation!=0){
            if(max(invar.loadings)==0){
              if(rotation=="varimax"){
                rot<-GPForth(Output[[ind_nclust]]$clusterspecific.loadings[[k]], method="varimax")
              } else if(rotation=="oblimin"){
                rot<-GPFoblq(Output[[ind_nclust]]$clusterspecific.loadings[[k]], method="oblimin")
              } else if(rotation=="geomin"){
                rot<-GPFoblq(Output[[ind_nclust]]$clusterspecific.loadings[[k]], method="geomin")
              } else if(rotation=="target"){
                rot<-GPFoblq(Output[[ind_nclust]]$clusterspecific.loadings[[k]], method="pst",methodArgs =list(W=targetW,Target=targetT))
              }
              rotatedloadings=rot$loadings
            }
            if(rotation=="varimax"){
              T_matrix=rot$Th
            } else {
              T_matrix=t(solve(rot$Th))
            }
            invT=solve(T_matrix)
            if(max(invar.loadings)>0){ # apply general rotation matrix to each cluster
              rotatedloadings=Output[[ind_nclust]]$clusterspecific.loadings[[k]]%*%T_matrix
            }
            if(rotation!="target"){ # reflect and permute in case of rotated EFA (but not target rotation)
              if(max(invar.loadings)>0){ # apply general reflection to each cluster
                refl=reflcommon
                rotatedloadings=matrix(refl,nrow=nvar,ncol=nfactors,byrow=TRUE)*rotatedloadings
              } else {
                if(k==1){
                  strongloadcutoff=mean(apply(abs(rotatedloadings),2,max))/2
                  strongloadind=which(abs(rotatedloadings)>=strongloadcutoff)
                  strongload=matrix(0,nrow=nvar,ncol=nfactors)
                  strongload[strongloadind]=rotatedloadings[strongloadind]
                  toreflect=colSums(strongload>0)<colSums(strongload<0)
                  refl=matrix(1,nrow=1,ncol=nfactors)
                  refl[toreflect]=-1
                  rotatedloadings=matrix(refl,nrow=nvar,ncol=nfactors,byrow=TRUE)*rotatedloadings
                  rotatedloadings_cl1=rotatedloadings
                } else {
                  agreem_cl1=matrix(0,nfactors,nfactors)
                  agreem_cl1_refl=agreem_cl1
                  for(q in 1:nfactors){
                    load=rotatedloadings[,q,drop=FALSE]
                    for(q1 in 1:nfactors){
                      agreem_cl1[q,q1]=sum(((load)-rotatedloadings_cl1[,q1])^2)
                      agreem_cl1_refl[q,q1]=sum(((-1*load)-rotatedloadings_cl1[,q1])^2)
                    }
                  }
                  mi=apply(agreem_cl1,2,sort)
                  mi_ind=apply(agreem_cl1,2,order)
                  mi_refl=apply(agreem_cl1_refl,2,sort)
                  mi_refl_ind=apply(agreem_cl1_refl,2,order)
                  refl=matrix(1,nrow=1,ncol=nfactors)
                  refl[mi[1,]>mi_refl[1,]]=-1
                  # browser()
                  perm=mi_ind[1,]
                  perm[refl==-1]=mi_refl_ind[1,refl==-1]
                  if (length(unique(perm))<length(perm)){
                    nonuniq=which(tabulate(perm)>1)
                    missing=which(is.element(1:nfactors,perm)==FALSE)
                    #perm[is.element(perm,nonuniq)]
                    for(nu in 1:length(nonuniq)){
                      fnu=nonuniq[nu]
                      pos=which(perm==fnu)
                      pos=pos[order(mi[1,perm==fnu])]
                      perm[pos[2:length(pos)]]=0
                    }
                    perm[perm==0]=missing # replacement is not optimized at this time
                  }
                  rotatedloadings=rotatedloadings[,perm]
                  rotatedloadings=matrix(refl,nrow=nvar,ncol=nfactors,byrow=TRUE)*rotatedloadings
                }
              }
              Rm=matrix(refl,nrow=nfactors,ncol=nfactors,byrow=TRUE)
            }
          } else {
            rotatedloadings=Output[[ind_nclust]]$clusterspecific.loadings[[k]]
            refl=matrix(1,nrow=1,ncol=nfactors)
            Rm=matrix(refl,nrow=nfactors,ncol=nfactors,byrow=TRUE)
            invT=diag(nfactors)
          }
          
          if(rescale==1){ # determine rescaling constant per factor, taking into account reflection and permutation
            if(k==1){
              if(rotation!="target"){ # also when rotation=0 (object can/should be previously rotated)
                markeritems=which(rowSums(markers)>0)
                indhl=apply(rotatedloadings[markeritems,],1,which.max) # index of factor each marker item has highest loading on
                ind1=apply(markers[markeritems,],1,which.max) # index of factor each item has a 1 for in 'markers'
                permm=match(ind1,indhl) # match the two to check whether permutation of 'markers' is needed
                if (length(unique(permm))<length(permm)){ # check whether this includes no duplicate elements and thus pertains to a proper permutation
                  nonuniq=which(tabulate(permm)>1)
                  missing=which(is.element(1:nfactors,permm)==FALSE)
                  for(nu in 1:length(nonuniq)){
                    fnu=nonuniq[nu] # non-unique element
                    pos=which(permm==fnu) # positions of non-unique element
                    load=rotatedloadings*markers[,permm] # corresponding loading values
                    load=load[markeritems[fnu],]
                    pos=pos[order(load[permm==fnu],decreasing=TRUE)]
                    permm[pos[2:length(pos)]]=0 # only keep element at position corresponding to highest loading
                  }
                  permm[permm==0]=missing # replacement is not optimized at this time
                }
                markers2=markers[,permm] # Note that markers2 is also used for k>1
              } else { # for target rotation, the permutation of 'markers' is already checked at the beginning, in correspondence to targetT
                markers2=markers
              }
            }
            markload=rotatedloadings*markers2 # rescaling also happens for fully specified target rotation
            markload=markload[markload!=0]
            sqrtFscale=diag(1/markload)
            invsqrtFscale=diag(markload)
            rotatedloadings=rotatedloadings%*%sqrtFscale
          }
          Output[[ind_nclust]]$clusterspecific.loadings[[k]]=rotatedloadings
          
          if(is.element("intercepts",cluster.spec)){ # factor means also need to be counterrotated/rescaled
            for(g in 1:ngroup){ # counter-rotate all sets of factor (co)variances and factor means
              crotatedFcov=invT%*%Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]]%*%t(invT)
              crotatedFmeans=Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]]%*%t(invT)
              if(k>1 && rotation!="target" && rotation!=0){
                crotatedFcov=crotatedFcov[perm,perm]
                crotatedFmeans=crotatedFmeans[,perm]
              }
              if(rotation!="target"){
                crotatedFcov=Rm*crotatedFcov*t(Rm)
                crotatedFmeans=refl*crotatedFmeans
              }
              if(rescale==1){ # correct for rescaling as well
                crotatedFcov=invsqrtFscale%*%crotatedFcov%*%invsqrtFscale
                crotatedFmeans=crotatedFmeans%*%invsqrtFscale
              }
              Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]]=crotatedFcov
              Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]]=crotatedFmeans
            }
          } else {
            for(g in 1:ngroup){ # counter-rotate all corresponding sets of factor (co)variances
              crotatedFcov=invT%*%Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]]%*%t(invT)
              if(k>1 && rotation!="target" && rotation!=0 && max(invar.loadings)==0){
                crotatedFcov=crotatedFcov[perm,perm]
              }
              if(rotation!="target"){
                crotatedFcov=Rm*crotatedFcov*t(Rm)
              }
              if(rescale==1){
                crotatedFcov=invsqrtFscale%*%crotatedFcov%*%invsqrtFscale
              }
              Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]]=crotatedFcov
            }
          }
        }
      }
      if(EFA==0 || nfactors==1 || (rotation=="target" & sum(targetT==0 & targetW==1)==sum(targetW))){ # in case of CFA or zero-approximating target rotation, only reflect
        for(k in 1:nclust){
          loadings=Output[[ind_nclust]]$clusterspecific.loadings[[k]]
          nvar=nrow(loadings)
          if(k==1){ # reflect to have positive strong loadings
            if(rotation=="target"){
              strongloadcutoff=mean(apply(abs(loadings),2,max))/2
              strongloadind=which(abs(loadings)>=strongloadcutoff)
              strongload=matrix(0,nrow=nvar,ncol=nfactors)
              strongload[strongloadind]=loadings[strongloadind]
            } else {
              strongload=loadings # for CFA, all non-zero loadings are considered strong loadings
            }
            toreflect=colSums(strongload>0)<colSums(strongload<0)
            refl=matrix(1,nrow=1,ncol=nfactors)
            refl[toreflect]=-1
            reflloadings=matrix(refl,nrow=nvar,ncol=nfactors,byrow=TRUE)*loadings
            reflloadings_cl1=reflloadings
          } else { # reflect to agreement with cluster 1
            agreem_cl1=matrix(0,nfactors,nfactors)
            agreem_cl1_refl=agreem_cl1
            for(q in 1:nfactors){
              load=loadings[,q,drop=FALSE]
              agreem_cl1[q,q]=sum(((load)-reflloadings_cl1[,q])^2)
              agreem_cl1_refl[q,q]=sum(((-1*load)-reflloadings_cl1[,q])^2)
            }
            refl=matrix(1,nrow=1,ncol=nfactors)
            refl[diag(agreem_cl1)>diag(agreem_cl1_refl)]=-1
            reflloadings=matrix(refl,nrow=nvar,ncol=nfactors,byrow=TRUE)*loadings
          }
          if(rescale==1 && EFA==0){ # determine rescaling constant per factor, taking into account reflection
            # for EFA with target rotation, rescaling already happened above
            markload=reflloadings*markers
            markload=markload[markload!=0]
            sqrtFscale=diag(1/markload)
            invsqrtFscale=diag(markload)
            reflloadings=reflloadings%*%sqrtFscale
          }
          Output[[ind_nclust]]$clusterspecific.loadings[[k]]=reflloadings
          Rm=matrix(refl,nrow=nfactors,ncol=nfactors,byrow=TRUE)
          for(g in 1:ngroup){ # reflect all corresponding sets of factor (co)variances
            Fcov=Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]]
            Fcov=Rm*Fcov*t(Rm)
            if(rescale==1 && EFA==0){
              Fcov=invsqrtFscale%*%Fcov%*%invsqrtFscale
            }
            Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]]=Fcov
            if(is.element("intercepts",cluster.spec)){ # factor means also need to be counterrotated/rescaled
              Fmeans=Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]]
              Fmeans=refl*Fmeans
              if(rescale==1 && EFA==0){
                Fmeans=Fmeans%*%invsqrtFscale
              }
              Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]]=Fmeans
            }
          }
        }
      }
      
      prefix="Factor"
      suffix=seq(1,nfactors)
      factorlabels=noquote(paste(prefix,suffix,sep="_"))
      for(k in 1:nclust){
        colnames(Output[[ind_nclust]]$clusterspecific.loadings[[k]])<-factorlabels
        for(g in 1:ngroup){
          colnames(Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]])<-factorlabels
          rownames(Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]])<-factorlabels
          if(is.element("intercepts",cluster.spec)){
            colnames(Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]])<-factorlabels
          }
        }
      }
      
    } else { # rotation/rescaling needs to be performed once, across all groups
      if(sum(design)==0){
        if(rotation!=0){
          if(rotation=="varimax"){
            rot<-GPForth(Output[[ind_nclust]]$invariant.loadings, method="varimax")
          } else if(rotation=="oblimin"){
            rot<-GPFoblq(Output[[ind_nclust]]$invariant.loadings, method="oblimin")
          } else if(rotation=="geomin"){
            rot<-GPFoblq(Output[[ind_nclust]]$invariant.loadings, method="geomin")
          } else if(rotation=="target"){
            rot<-GPFoblq(Output[[ind_nclust]]$invariant.loadings, method="pst",methodArgs =list(W=targetW,Target=targetT))
          }
          rotatedloadings=rot$loadings
          nvar=nrow(rotatedloadings)
          if(rotation=="varimax"){
            T_matrix=rot$Th
          } else {
            T_matrix=t(solve(rot$Th))
          }
          invT=solve(T_matrix)
          if(rotation!="target" || (rotation=="target" & sum(targetT==0 & targetW==1)==sum(targetW))){ # reflect in case of rotated EFA (not for target rotation, unless it only approximates zeros)
            strongloadcutoff=mean(apply(abs(rotatedloadings),2,max))/2
            strongloadind=which(abs(rotatedloadings)>=strongloadcutoff)
            strongload=matrix(0,nrow=nvar,ncol=nfactors)
            strongload[strongloadind]=rotatedloadings[strongloadind]
            toreflect=colSums(strongload>0)<colSums(strongload<0)
            refl=matrix(1,nrow=1,ncol=nfactors)
            refl[toreflect]=-1
            rotatedloadings=matrix(refl,nrow=nvar,ncol=nfactors,byrow=TRUE)*rotatedloadings
          } else {
            refl=matrix(1,nrow=1,ncol=nfactors)
          }
        } else {
          rotatedloadings=Output[[ind_nclust]]$invariant.loadings
          refl=matrix(1,nrow=1,ncol=nfactors)
          invT=diag(nfactors)
        }
        Rm=matrix(refl,nrow=nfactors,ncol=nfactors,byrow=TRUE)
        if(rescale==1){
          if(rotation!="target"){
            markeritems=which(rowSums(markers)>0)
            indhl=apply(rotatedloadings[markeritems,],1,which.max) # index of factor each marker item has highest loading on
            ind1=apply(markers[markeritems,],1,which.max) # index of factor each item has a 1 for in 'markers'
            permm=match(ind1,indhl) # match the two to check whether permutation of 'markers' is needed
            if (length(unique(permm))<length(permm)){ # check whether this includes no duplicate elements and thus pertains to a proper permutation
              nonuniq=which(tabulate(permm)>1)
              missing=which(is.element(1:nfactors,permm)==FALSE)
              for(nu in 1:length(nonuniq)){
                fnu=nonuniq[nu] # non-unique element
                pos=which(permm==fnu) # positions of non-unique element
                load=rotatedloadings*markers[,permm] # corresponding loading values
                load=load[markeritems[fnu],]
                pos=pos[order(load[permm==fnu],decreasing=TRUE)]
                permm[pos[2:length(pos)]]=0 # only keep element at position corresponding to highest loading
              }
              permm[permm==0]=missing # replacement is not optimized at this time
            }
            markers2=markers[,permm]
          } else {
            markers2=markers
          }
          markload=rotatedloadings*markers2
          markload=markload[markload!=0]
          sqrtFscale=diag(1/markload)
          invsqrtFscale=diag(markload)
          rotatedloadings=rotatedloadings%*%sqrtFscale
        }
        Output[[ind_nclust]]$invariant.loadings=rotatedloadings
        
        if(is.element("intercepts",cluster.spec)){ # intercepts are cluster-specific, so factor means are group-and-cluster-specific
          for(g in 1:ngroup){ # counter-rotate all sets of factor (co)variances and factor means
            crotatedFcov=invT%*%Output[[ind_nclust]]$groupspecific.factorcovariances[[g]]%*%t(invT)
            crotatedFcov=Rm*crotatedFcov*t(Rm)
            if(rescale==1){
              crotatedFcov=invsqrtFscale%*%crotatedFcov%*%invsqrtFscale
            }
            Output[[ind_nclust]]$groupspecific.factorcovariances[[g]]=crotatedFcov
            for(k in 1:nclust){
              crotatedFmeans=Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]]%*%t(invT)
              crotatedFmeans=refl*crotatedFmeans
              if(rescale==1){
                crotatedFmeans=crotatedFmeans%*%invsqrtFscale
              }
              Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]]=crotatedFmeans
            }
          }
        } else { # intercepts are invariant, so factor means are group-specific
          for(g in 1:ngroup){ # counter-rotate all sets of factor (co)variances and factor means
            crotatedFcov=invT%*%Output[[ind_nclust]]$groupspecific.factorcovariances[[g]]%*%t(invT)
            crotatedFmeans=Output[[ind_nclust]]$groupspecific.factormeans[[g]]%*%t(invT)
            crotatedFcov=Rm*crotatedFcov*t(Rm)
            crotatedFmeans=refl*crotatedFmeans
            if(rescale==1){
              crotatedFcov=invsqrtFscale%*%crotatedFcov%*%invsqrtFscale
              crotatedFmeans=crotatedFmeans%*%invsqrtFscale
            }
            Output[[ind_nclust]]$groupspecific.factorcovariances[[g]]=crotatedFcov
            Output[[ind_nclust]]$groupspecific.factormeans[[g]]=crotatedFmeans
          }
        }
      }
      if(sum(design)>0 || nfactors==1){ # in case of CFA, only reflect
        loadings=Output[[ind_nclust]]$invariant.loadings
        nvar=nrow(loadings)
        strongload=loadings # for CFA, all non-zero loadings are considered strong loadings
        toreflect=colSums(strongload>0)<colSums(strongload<0)
        refl=matrix(1,nrow=1,ncol=nfactors)
        refl[toreflect]=-1
        reflloadings=matrix(refl,nrow=nvar,ncol=nfactors,byrow=TRUE)*loadings
        if(rescale==1){
          markload=reflloadings*markers
          markload=markload[markload!=0]
          sqrtFscale=diag(1/markload)
          invsqrtFscale=diag(markload)
          reflloadings=reflloadings%*%sqrtFscale
        }
        Output[[ind_nclust]]$invariant.loadings=reflloadings
        
        Rm=matrix(refl,nrow=nfactors,ncol=nfactors,byrow=TRUE)
        for(g in 1:ngroup){ # reflect all corresponding sets of factor (co)variances and factor means
          Fcov=Output[[ind_nclust]]$groupspecific.factorcovariances[[g]]
          Fcov=Rm*Fcov*t(Rm)
          if(rescale==1){
            Fcov=invsqrtFscale%*%Fcov%*%invsqrtFscale
          }
          Output[[ind_nclust]]$groupspecific.factorcovariances[[g]]=Fcov
          if(is.element("intercepts",cluster.spec)){
            for(k in 1:nclust){
              Fmeans=Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]]
              Fmeans=refl*Fmeans
              if(rescale==1){
                Fmeans=Fmeans%*%invsqrtFscale
              }
              Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]]=Fmeans
            }
          } else {
            Fmeans=Output[[ind_nclust]]$groupspecific.factormeans[[g]]
            Fmeans=refl*Fmeans
            Output[[ind_nclust]]$groupspecific.factormeans[[g]]=Fmeans
          }
        }
      }
      
      prefix="Factor"
      suffix=seq(1,nfactors)
      factorlabels=noquote(paste(prefix,suffix,sep="_"))
      colnames(Output[[ind_nclust]]$invariant.loadings)<-factorlabels
      for(g in 1:ngroup){
        colnames(Output[[ind_nclust]]$groupspecific.factorcovariances[[g]])<-factorlabels
        rownames(Output[[ind_nclust]]$groupspecific.factorcovariances[[g]])<-factorlabels
        if(is.element("intercepts",cluster.spec)){
          for(k in 1:nclust){
            colnames(Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]])<-factorlabels
          }
        } else {
          colnames(Output[[ind_nclust]]$groupspecific.factormeans[[g]])<-factorlabels
        }
      }
    }
  }
  
  if(length(OutputObject)==2){
    OutputObject$MMGFAsolutions=Output
  } else {
    OutputObject=Output
  }
  return(OutputObject)
}


# Step 1 function --------------------------------------------------------------------------------
MixMix_Step1 <- function(data, step1model, group = "group", MM.cluster.spec=c("loadings"),
                         MM.nclus = 1, MM.maxiter = 10000, MM.nruns = 50, MM.design = loading_matrix, invar_loadings = invar_loadings,
                         markers, seed = 100){
  start_time_step1 <- Sys.time()  # Start time for Step 1
  
  # centered data
  g_name <- as.character(unique(data[, group]))
  vars <- lavNames(lavaanify(step1model, auto = TRUE)) #observed var
  lat_var <- lavNames(lavaanify(step1model, auto = TRUE), "lv") #latent var
  
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
  
  # Get sample covariance matrix per group (used later)
  S_unbiased <- lapply(X = unique(centered[, group]), FUN = function(x) {cov(centered[centered[, group] == x, vars])})
  
  # browser()
  set.seed(seed)
  #MMGFA: loadings cluster-specific; 1-6 clusters;
  output1 <-mixmgfa(data = centered, N_gs = N_gs, nfactors = nfactors, 
                    cluster.spec = MM.cluster.spec, nsclust = MM.nclus, 
                    maxiter = MM.maxiter, nruns = MM.nruns, design = MM.design, invar_loadings = invar_loadings)
  # browser()
  # if(length(MM.nclus)>1){
  #   MM.nclus.final <- which.max(output1[["overview"]][, "screeratios"])
  # } else {
  #   MM.nclus.final <- MM.nclus
  # }
  
  #MMGFA: rescaling factors using marker variables
  output2 <- ScaleRotateMixmgfa_pinvar(output1, N_gs = N_gs,
                                       nsclust = MM.nclus, design = MM.design, rescale=1, markers = markers,
                                       rotation=0,targetT=0,targetW=0)
  # browser()
  # Get the list of cluster names dynamically from MMGFAsolutions
  cluster_names <- grep("^\\d+\\.clusters$", names(output2[["MMGFAsolutions"]]), value = TRUE)
  
  # Initialize lists to store extracted information and results for each cluster
  # factor_cov_list <- list()
  # cluster_memb_list <- list()
  # lambda_list <- list()
  # theta_list <- list()
  # cov_eta <- list()  # cov_eta as a list for each cluster
  # 
  
  # Loop over each cluster by name and extract data
  # for (nk in cluster_names) {
  uniquevar_key <- grep("uniquevariances$", names(output2[["MMGFAsolutions"]][[cluster_names]]), value = TRUE)
  
  # Extract cluster-specific information
  factor_cov_list<- output2[["MMGFAsolutions"]][[cluster_names]][["group.and.clusterspecific.factorcovariances"]]
  cluster_memb_list <- output2[["MMGFAsolutions"]][[cluster_names]][["clustermemberships"]]
  lambda_list <- output2[["MMGFAsolutions"]][[cluster_names]][["clusterspecific.loadings"]]
  theta_list <- output2[["MMGFAsolutions"]][[cluster_names]][[uniquevar_key]] #group or cluster
  
  # Initialize an array for the weighted sum specific to this cluster
  weighted_sum_mg <- array(data = NA, dim = c(nfactors, nfactors, ngroups))
  
  # Compute weighted sum for this cluster
  for (j in 1:ngroups) {
    # Multiply each cluster-specific matrix by its corresponding cluster membership and sum the results
    weighted_sum_mg[,,j] <- Reduce(`+`, Map(`*`, factor_cov_list[j, ], cluster_memb_list[j, ]))
  }
  
  # Store the result in cov_eta list for this cluster
  cov_eta <- weighted_sum_mg 
  dimnames(cov_eta)[[1]] <- dimnames(cov_eta)[[2]] <- lat_var
  
  # ######################################################################################################
  # ## single-indicator approach with MM_nclus
  # data2 <- split(centered[,2:21], centered[,1])
  # nfactors <- 4
  # # ngroups <- length(data2)
  # # browser()
  # 
  # # Extract cluster-specific covariance matrices (column-based extraction)
  # cov_matrix <- factor_cov_list #s1output[[paste0(MM_nclus, " clusters")]][["mmgfa_output"]][["factor_cov"]]
  # cov_list <- lapply(1:MM.nclus, function(k) as.list(cov_matrix[, k]))
  # 
  # # Function to extract list-based parameters (lambda and theta)
  # extract_param <- function(param_name, cluster_num) {
  #   # s1output[[paste0(MM_nclus, " clusters")]][["mmgfa_output"]][[param_name]]
  #   param_name[[paste0("Cluster_", cluster_num)]]
  # }
  # 
  # # Extract lambda and theta for all clusters
  # lambda_list_fs <- lapply(1:MM.nclus, function(k) extract_param(lambda_list, k))
  # theta_list_fs <- theta_list[,1] #lapply(1:MM.nclus, function(k) extract_param(theta_list, k))
  # 
  # # Extract cluster membership matrix
  # z_gk <- cluster_memb_list#s1output[[paste0(MM_nclus, " clusters")]][["mmgfa_output"]][["cluster_memb"]]
  # 
  # # Function to compute factor scores (FS)
  # compute_fs <- function(cov_k, lambda_k, theta_k, data2) {
  #   lapply(1:ngroups, function(x) {
  #     fsmatrix <- diag(diag(cov_k[[x]])) %*% t(lambda_k) %*% 
  #       solve(lambda_k %*% diag(diag(cov_k[[x]])) %*% t(lambda_k) + theta_k[[x]])
  #     as.matrix(data2[[x]]) %*% t(fsmatrix)
  #   })
  # }
  # 
  # # Compute FS for all clusters
  # fs_list <- mapply(compute_fs, cov_list, lambda_list_fs, MoreArgs = list(theta_k = theta_list_fs, data2 = data2), SIMPLIFY = FALSE)
  # 
  # # Compute weighted FS
  # fs_weighted <- lapply(1:ngroups, function(x) {
  #   Reduce(`+`, lapply(1:MM.nclus, function(k) fs_list[[k]][[x]] * z_gk[x, k]))
  # })
  # 
  # # Function to compute FS variance
  # compute_fsvar <- function(cov_k, lambda_k, theta_k) {
  #   lapply(1:ngroups, function(x) {
  #     diag(diag(cov_k[[x]])) - 
  #       t(diag(diag(cov_k[[x]]))) %*% t(lambda_k) %*% 
  #       solve(lambda_k %*% diag(diag(cov_k[[x]])) %*% t(lambda_k) + theta_k[[x]]) %*% 
  #       lambda_k %*% diag(diag(cov_k[[x]]))
  #   })
  # }
  # 
  # # Compute FS variance for all clusters
  # fsvar_list <- mapply(compute_fsvar, cov_list, lambda_list_fs, MoreArgs = list(theta_k = theta_list_fs), SIMPLIFY = FALSE)
  # 
  # # Compute weighted FS variance
  # fsvar_weighted <- lapply(1:ngroups, function(x) {
  #   sqrt(Reduce(`+`, lapply(1:MM.nclus, function(k) fsvar_list[[k]][[x]] * z_gk[x, k])))
  # })  
  # 
  # 
  # ####SINGLE-INDICATOR
  # #single indicator all factors together
  # postmeanfs_list <- fs_weighted #split(fs_ind[, c(26, 28,30,32)], fs_ind[,1])
  # # postvarfs_list <- split(fs_ind[, c(27, 29,31,33)], fs_ind[,1])
  # 
  # # factor variances computed based on factor scores and their variances
  # fvar_compute <- replicate(nfactors, numeric(ngroups), simplify = FALSE)
  # 
  # for (f in 1:nfactors) {
  #   for (g in 1:ngroups) {
  #     fvar_compute[[f]][g] <- var(postmeanfs_list[[g]][,f]) + fsvar_weighted[[g]][f,f]^2
  #   }
  # }
  # 
  # # single-indicator loadings (rhos)
  # rho <- lapply(1:nfactors, function(i) sapply(1:ngroups, function(j) var(postmeanfs_list[[j]][,i]) / fvar_compute[[i]][j]))
  # lambda <- array(0, dim = c(nfactors, nfactors, ngroups))
  # for (j in 1:ngroups) lambda[,,j] <- diag(sapply(rho, `[`, j))
  # dimnames(lambda)[[1]] <- c("x", "z", "m", "y")
  # dimnames(lambda)[[2]] <- c("F1", "F2", "F3", "F4")  
  # 
  # # single-indicator unique variances
  # resvar <- lapply(1:ngroups, function(j) sapply(1:nfactors, function(i) fvar_compute[[i]][j] * rho[[i]][j] * (1 - rho[[i]][j])))
  # theta <- array(unlist(lapply(1:ngroups, function(j) {
  #   diag(resvar[[j]])  # Create diagonal matrix for each group
  # })), dim = c(nfactors, nfactors, ngroups))
  # dimnames(theta)[[1]] <- dimnames(theta)[[2]] <-  c("x", "z", "m", "y") 
  # 
  # # single-indicator data cov
  # fcovdat_cov <- lapply(postmeanfs_list, function(x) cov(x[,1:4]))
  # fcovdat <- array(unlist(fcovdat_cov), dim = c(nfactors, nfactors, ngroups))
  # 
  # # single-indicator fcov: Croon's formula
  # cov_eta_fs <-  array(data = 0, dim = c(nfactors, nfactors, ngroups))
  # for (j in 1:ngroups) {
  #   cov_eta_fs[,,j] <- solve(lambda[,,j]) %*% (fcovdat[,,j] - theta[,,j]) %*% solve(t(lambda[,,j]))
  # }
  # dimnames(cov_eta_fs)[[1]] <- dimnames(cov_eta_fs)[[2]] <- c("F1", "F2", "F3", "F4") 
  
  
  end_time_step1 <- Sys.time()  # End time for Step 1
  step1_time <- difftime(end_time_step1, start_time_step1, units = "mins")
  # cat("Step 1 completed in", step1_duration, "minutes\n")
  
  
  return(list(# data = centered,
    cov_eta = cov_eta,
    ngroups = ngroups,
    N_gs = N_gs,
    S_unbiased   = S_unbiased,
    vars         = vars, #var names
    lat_var      = lat_var, #latent variable names
    mmgfa_output = list(output1 = output1,
                        output2 = output2,
                        # MM.nclus.final = MM.nclus.final,
                        cluster_memb = cluster_memb_list,
                        factor_cov = factor_cov_list, #gk
                        lambda_gs = lambda_list, #gk
                        theta_gs = theta_list), #g
    # ##single-indicator FS
    # fs_weighted = fs_weighted,
    # fsvar_weighted = fsvar_weighted,
    # cov_eta_fs = cov_eta_fs,
    step1_time = step1_time))
  
}

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

#marker variable 
markers <- loading_matrix
markers[1:2,1] <- markers[c(4,6,7),2] <- markers[c(8,9,11:13),3]<- 0

#for the simulated data
# data <- data[,c(3:24)]
# write.csv(data, file = "school_level_data.csv")
# data <- as.data.frame(data)
# which(table(data$CNTSCHID)<25)

# s1output_3k <- list()
# for (i in 3) {
start.time <- Sys.time()
  s1output_4MMk <- MixMix_Step1(data[,1:14], step1model = mm_model, group = "CNTSCHID", MM.cluster.spec=c("loadings"),
                                MM.nclus = 4, MM.maxiter = 100000, MM.nruns = 100, MM.design = loading_matrix,
                                invar_loadings = c(1,1,1,
                                                   0,1,0,1,
                                                   1,1,1,0,0,0),
                                markers = markers, seed = 100)
# }
end.time <- Sys.time()
diff_time <- difftime(end.time, start.time, units='mins')
# Time difference of 17.34394 mins
save(s1output_4MMk, file = "s1output_4MMk_724.RData")

uncertainty <- 1 - apply(s1output_4MMk[["mmgfa_output"]][["cluster_memb"]], 1, max)
meanuncertain <- mean(uncertainty)
# 0.1193098
clustering <- t(apply(s1output[[2]][["mmgfa_output"]][["cluster_memb"]], 1, function(x) as.numeric(x == max(x))))
colMeans(clustering)


# ##check MM clusters 1-5
# # centered data
# g_name <- as.character(unique(data[, group]))
# vars <- lavNames(lavaanify(step1model, auto = TRUE)) #observed var
# lat_var <- lavNames(lavaanify(step1model, auto = TRUE), "lv") #latent var
# 
# centered <- data
# group.idx <- match(data[,group], g_name)
# group.sizes <- tabulate(group.idx)
# group.means <- rowsum.default(as.matrix(data[,vars]),
#                               group = group.idx, reorder = FALSE,
#                               na.rm = FALSE)/group.sizes
# centered[,vars] <- data[,vars] - group.means[group.idx, ,drop = FALSE]
# 
# N_gs <- group.sizes #as.numeric(table(data[,group]))
# nfactors <- length(lat_var) #ncol(MM.design) #4 lv
# ngroups <- length(N_gs)
# 
# set.seed(100)
# output1 <-mixmgfa(data = data, N_gs = N_gs, nfactors = nfactors, 
#                   cluster.spec = MM.cluster.spec, nsclust = c(1:6), 
#                   maxiter = MM.maxiter, nruns = MM.nruns, design = MM.design, invar_loadings = invar_loadings)
# plot(output1$overview)
# 
# result <- as.matrix(output1[["overview"]])
# 
# plot(x = c(1:6), y = result[,2]); lines(x = c(1:6), y = result[,2]) #LL ?
# plot(x = c(1:6), y = result[,4]); lines(x = c(1:6), y = result[,4]) #BIC_N 2
# plot(x = c(1:6), y = result[,5]); lines(x = c(1:6), y = result[,5]) #BIC_G 2/3
# plot(x = c(1:6), y = result[,6]); lines(x = c(1:6), y = result[,6]) #AIC ?


# nr of clusters    loglik nrpars    BIC_N    BIC_G      AIC screeratios
# [1,]              1 -295094.3  26339 841796.0 751793.3 642866.7          NA
# [2,]              2 -294960.8  26350 841633.9 751593.7 642621.5    3.193629
# [3,]              3 -294880.9  26371 841674.9 751562.8 642503.8    1.178097
# [4,]              4 -294832.5  26386 841721.3 751558.0 642437.0    1.058836
# [5,]              5 -294792.9  26399 841766.3 751558.6 642383.8          NA
# convergence nr.activated.constraints
# [1,]           1                       11
# [2,]           1                       17
# [3,]           1                       13
# [4,]           1                       15
# [5,]           1                       19
