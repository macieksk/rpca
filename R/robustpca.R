######################
### TRPCA FUNCTION ###
######################

#library(compiler)  #Include directive in DESCRIPTION
#library(irlba)     #Include directive in DESCRIPTION
#library(Matrix)    #Include directive in DESCRIPTION

F2norm<-cmpfun(function(M)sqrt(sum(M^2)))
l2.norm <- F2norm #cmpfun(function(x){ sqrt(sum(x^2)) })
               
thr.op<-cmpfun(function(x,thr){sign(x)*pmax(abs(x)-thr,0)})
thr.op.sparse<-cmpfun(function(x,thr){
    x[(x<=thr)&(x>=(-thr))]<-0
    x=x-(x>thr)*thr
    x=x+(x<(-thr))*thr
})                          
              
trpca.thresh.nuclear<-cmpfun(function(M, curr.imu, k, prev.res = NULL,
                                      maxit=200,
                                      thr.op.fun=thr.op,
                                      matrix.mul.fun=function(a,b){a%*%b}){
    s<-irlba(M,k+1,maxit=maxit,work=k+20,v=prev.res)
    s$thr <- tail(s$d,1)
    if (is.na(curr.imu)) curr.imu<-mean(tail(s$d,2))
    if (s$thr>curr.imu) {
        s$k.too.small<-TRUE        
        return(s)
    }
    s$k.too.small<-FALSE 
    s$thr <- curr.imu                                          
    s$vt<-t(s$v)
    dd<-thr.op.fun(s$d,s$thr)
    id<-which(dd!=0)
    s$d<-dd[id]
    s$u<-s$u[,id,drop=FALSE]
    s$vt<-s$vt[id,,drop=FALSE]
    s$v<-s$v[,id,drop=FALSE]    
    s$L<-matrix.mul.fun(s$u, s$d*s$vt)  # s$u%*%(s$d*s$vt)
    s
})      

#Sparse version in fact requires sparse svd decomposition ? #TOTHINK #TODO later
m2sparseM<-cmpfun(function(m)sparseMatrix(i=row(m),j=col(m),x=as.vector(m)))

trpca.thresh.nuclear.sparse1<-cmpfun(function(...)trpca.thresh.nuclear(...,
                                                         thr.op.fun=thr.op.sparse,
                                                         matrix.mul.fun=function(a,b){m2sparseM(a)%*%b}))

# The one below is probably faster than the above for large matrices due to Matrix sparse matrices not being parallel
trpca.thresh.nuclear.sparse2<-cmpfun(function(...)trpca.thresh.nuclear(...,
                                                         thr.op.fun=thr.op.sparse,
                                                         matrix.mul.fun=function(a,b)m2sparseM(a%*%b)))
                                                                      
thresh.l1<-cmpfun(function(x,thr){
    thr.op(x,thr)
})
                      
thresh.l1.sparse<-cmpfun(function(x,thr){
    thr.op.sparse(x,thr)
})                      

zero.matrix<-function(n,m)matrix(0,nrow=n,ncol=m)                                                                       
zero.matrix.sparse<-function(n,m)sparseMatrix(i=numeric(0),j=numeric(0),x=numeric(0),dims=c(n,m))                                                                       

trpca<-function(M, #k,
             k.start=1,
             lambda = 1/sqrt(max(dim(M))), #This is ok only for dense matrices
             lambda2 = 100*lambda, #TODO needs proper L1 sparse vs L2 noise weight setting
             L2noise = TRUE, #Do decompose into M=L+S+E, or just M=L+S if FALSE
             mu = prod( dim(M)) / (4*sum(abs(M)) ), #This is ok only for dense matrices
             mu.max = mu*100,  #Stops mu from getting to large too fast 
                               #(i.e. from caring too much for constraint than objective.function)
             mu.min = mu/200,   #If smallest computed SV is larger than 1/mu.min we increase k.current
                               #and compute one more SV in next iteration
             mu.growth.ratio=1.1,
             term.delta=10^(-7),
             max.iter=5000,
             trace=FALSE,
             message.iter=100,
             n.iter.without.L2noise=5, #Number of start iterations without decomposing L2noise
             #thresh.nuclear.fun=trpca.thresh.nuclear.sparse2,
             #thresh.l1.fun=thresh.l1.sparse,
             #zero.matrix.fun=zero.matrix.sparse,             
             thresh.nuclear.fun=trpca.thresh.nuclear,
             thresh.l1.fun=thresh.l1,
             zero.matrix.fun=zero.matrix,                
             F2norm.fun=F2norm){
    tolerance <- term.delta # used for E, (1-dd)>0+tolerance calculation
    dm<-dim(M)
    term.norm<-term.delta*F2norm.fun(M)
        
    k.current<-k.start
    #q<-(max.iter/2)^(1/((k-k.start+1))) #No longer used, TO DO FIX to match approx. max.iter
    #local.max.iter.phase<-2*q
    #local.max.iter<-0+local.max.iter.phase
             
    stats<-c()
    stats.objective.f<-c()
    jumps<-c()
    imu.vec<-c()
    L.rank.vec<-c()
        
    #S<-Yimu<-zero.matrix.fun(dm[1],dm[2])
    #MS<-M
    #resid.norm<-NA
    #objective.f.value<-NA    
    
    imu<-1/mu.max     
    
    j<-0
    converged<-FALSE
    #while (k.current<=k && !converged) {
    while (j<max.iter && !converged) {
    message(paste("***** \n LOOP WITH k=",k.current,sep="")) #", q=",q
    
    if( j==0 || (imu>=(1/mu.min))) {
      L<-S<-Yimu<-zero.matrix.fun(dm[1],dm[2])
      MS<-M
      MLS<-M
      #MSE <- M # M-S-E
      imu<-NA
      imu.old<-NA
      resid.norm<-NA
      objective.f.value<-NA
      dd <- 1
    }
        
    jumps<-c(jumps,length(stats)+1)    
       
    i<-(-1)    
    converged<-FALSE
    while(TRUE){        
        i<-i+1
        j<-j+1
        #L.svd<-thresh.nuclear.fun(MS+Yimu, imu, k.current)
        #L.svd<-thresh.nuclear.fun(MSE+Yimu, imu, k.current)
        
        YYimu = Yimu + MLS  # <--> YYimu = Yimu + (M-L-S)  
        # thus m*YYimu can be used in E calculation

        #E=(1-l2/|Y + m(M-L-S)|)/m * (Y + m(M-L-S))
        #E=(1-l2/(m|YYimu|))/m * m * YYimu = (1-l2/(m*|YYimu|))*YYimu
        # L.svd = thr.nuclear.l1( Yimu + MLS - E ) = thr.nuclear.l1(YYimu-(1-l2/(m*|YYimu|))*YYimu)

        # 1. L opt step
        # L.svd = thr.nuclear.l1( (lamda2/((1/imu) * l2.norm(YYimu)))*YYimu + L )
        L.svd <- thresh.nuclear.fun( dd*YYimu + L, imu, k.current)
        
        if(trace & (i%%message.iter==0)){
          print("----- SUMMARY AFTER L opt step -----")
          print(rbind(L.svd.d=round(L.svd$d,4)))
          cat(paste("\nmu.min=", round(mu.min,4), ", mu=", round(1/imu,4), 
                    ", mu.max=", round(mu.max,4), ", imu=", round(imu,4), sep=""))
          cat(paste("\nk.too.small = ", L.svd$k.too.small,", L.svd.thr = ",round(L.svd$thr,4)))
          #print(data.frame(dd=dd,normYYimu=l2.norm(YYimu),normL=F2norm(L)))
          cat("\n----- ----- ----- ----- ----- -----")
        }
        #stopifnot(i < 10)
        
        if (is.na(imu)) { 
            #First iteration-like case            
            
            stopifnot(sum(Yimu)==0)
            imu<-L.svd$thr                        
            imu.old<-imu
            limu<-lambda*imu
            
            #TOTHINK what if imu is still very small: imu<1/mu.max ? 
            #--> it will be grown to at.least 1/mu.max below in "else" clause
            if (trace & (i%%message.iter==0)) { #DEBUG
            print("----- First iteration like case is.na(imu) -----")
            print(data.frame(k.too.small=L.svd$k.too.small,L.svd.thr=L.svd$thr))
            ##TODO optimize below l2.norm computation with the crossprod trick
            print(data.frame(dd=dd,normYYimu=l2.norm(YYimu),normL=F2norm(L),
                mu=1/imu,imu=imu))
            print("----- ----- ----- -----  ----- ----- ----- -----")
            }
        }
        if (L.svd$k.too.small) { #k.current too small <-> imu too large (mu to small) 
                                 #though is it possible that imu is too small? #TODO #TOTHINK
          #Strategy: imu grows to current (k+1)th sv until either it reaches 1/mu.min, and then k->k+1,
          #          or k.too.small==FALSE
          imu.old<-imu
          #imu<-min(L.svd$thr,1/mu.min)          
          imu<-min(tail(L.svd$d,1),1/mu.min)
          k.current = k.current + 2  
          stopifnot(L.svd$thr==tail(L.svd$d,1))
        } else {
            k.current <- min(min(dim(M))-2,length(L.svd$d)+1)
            ##########################################
            # Thus general loop structure looks like this

            # Moved to the begining of while: YYimu = Yimu + MLS  # <--> YYimu = Yimu + (M-L-S)  

            # 1. L opt step
            # L.svd = thr.nuclear.l1( (lamda2/((1/imu) * l2.norm(YYimu)))*YYimu + L )
            # Moved to the begining of while: L.svd = thresh.nuclear.fun((lambda2/((1/imu) * l2.norm(YYimu)))*YYimu + L, imu, k.current)
            L <- L.svd$L

            # 2. E opt step
            # Update E by updating YYimu 
            YYimu <- Yimu + MS - L        # newL from above
            #normYYimu <- l2.norm(YYimu) # sqrt(sum(YYimu^2))
            ## The trick below is to change dimension to long vector
            ## and compute crossprod which uses matrix operations
            ## (BTW if this is done inside a code block, like a function
            ## it will cause long rewrite of the whole matrix, however this
            ## code below shouldn't)
            .dYYimu<-dim(YYimu)
            dim(YYimu)<-c(prod(.dYYimu),1)
            normYYimu<-sqrt(crossprod(YYimu))
            dim(YYimu)<-.dYYimu
            ##
            stopifnot(normYYimu != 0)
            # IGNORE l2 NOISE 
            if (j<n.iter.without.L2noise | ! L2noise ) {
                dd<-1                
            } else {
                dd <- lambda2*imu / normYYimu            
                dd<-min(1,dd)
                #print(c(dd=dd))
                #stopifnot(dd<(1-tolerance))
            }

            # 3. S opt step
            ### argmin S =   takie samo wyliczenie
            ### = thr.op.l1( (l2/(m*|YYimu|))*YYimu + S )
            ### S = thr.op.l1( dd*YYimu + S )
            S <- thresh.l1.fun( dd*YYimu + S, limu)

            # Update help variables 
            MS <- M-S
            MLS<- MS-L

            # 4. E opt step
            # Update E by updating YYimu 
            YYimu <- Yimu + MLS
            #normYYimu <- l2.norm(YYimu)
            ##For the description of the trick read some lines above
            .dYYimu<-dim(YYimu)
            dim(YYimu)<-c(prod(.dYYimu),1)
            normYYimu<-sqrt(crossprod(YYimu))
            dim(YYimu)<-.dYYimu
            ##
            stopifnot(normYYimu != 0)            
            # IGNORE L2 NOISE 
            if (j<n.iter.without.L2noise | ! L2noise ) {
                dd<-1                
            } else {
                dd <- lambda2*imu/normYYimu
                if (dd > 1 ) cat(paste("\n[DD INFO] After optimization dd > 1 (dd = ", round(dd,4), ")", sep=""))
                dd<-min(1,dd)
                
                #stopifnot(dd<(1-tolerance))
            }


            # Update Yimu and MLSE

            # MLSE = MLS-E # This can be computed more efficiently
            # E=(1-dd)*YYimu = (1-dd) * (Yimu+M-L-S)
            # M-L-S-E = M-L-S- (1-dd)* (Yimu+M-L-S) = M-L-S -(1-dd)*(M-L-S) - (1-dd)*Yimu
            # MLSE    = dd *(M-L-S) - (1-dd)*Yimu

            #MLSE = dd * (MLS) - (1-dd)*Yimu
            #Yimu = Yimu + MLSE
            
            #Wait! The above two lines can be improved!
            #Yimu = Yimu + dd * MLS - (1-dd)*Yimu
            #Yimu =  dd*MLS + dd*Yimu

            Yimuold <- Yimu
            Yimu <- dd*YYimu #= dd*(Yimu+MLS) 
            #if (i==0) Yimu<-Yimu/2
            #This above one line is enough, though MLSE needs to be computed for convergence check
            MLSE = Yimu - Yimuold
                    
          #Old code TODELETE? 
          # L<-L.svd$L
          # #S<-thresh.l1.fun(M-L+Yimu,limu)
          # C <- M - L + Yimu  
          # S <- thresh.l1.fun(C, limu)# - E,limu)
          
          # #E <- E# lambda2/( (1-1/imu)* sqrt(((C - S)/imu)^2)) * (C-S/imu)
            
          # MS = M-S
          # MLS = MS-L
          # #MLSE = MLS - E
            
          ## compute resid.norm every 3rd iteration for speedup
          if(i%%3 == 0 || i==max.iter){              
            #resid.norm<-F2norm.fun(MLSE)
            ##For the description of the trick read some lines above
            .dMLSE<-dim(MLSE)
            dim(MLSE)<-c(prod(.dMLSE),1)
            resid.norm<-sqrt(crossprod(MLSE))
            dim(MLSE)<-.dMLSE
            ##
          } else { 
            resid.norm<-NA
          }
            
          if(trace){
            ## compute every 6th iteration for speedup
            l2.normE<-(1-dd)*normYYimu
            #objective.f.value<-sum(L.svd$d)+lambda*sum(abs(S))+lambda2*sqrt(sum(E^2))
            if(i%%6==0){
                objective.f.value<-c(sum(L.svd$d),lambda*sum(abs(S)),lambda2*l2.normE)
                cat(paste("\nScaled norm values: Lnorm = ", round(objective.f.value[1], 4), 
                          ", Snorm = ", round(objective.f.value[2],4), 
                          ", Enorm = ", round(objective.f.value[3],4)))
                #objective.f.value<-c(sum(L.svd$d),lambda*sum(abs(S)),lambda2*l2.norm(E^2))
            }
          }
          # Yimu = Yimu + MLS#E      
          
          # #Strategy: grow mu <--> shrink imu
          imu.old<-imu
          imu<-max(1/mu.max,imu/mu.growth.ratio)  
          # #TODO consider not growing k.current if everything goes well <-> if reached here  
          # #QUESTION: It is worth to consider if k.current should be set to length(L.svd$L$d)+C (here and further in the code)
        }
        if (imu!=imu.old) {
          limu<-lambda*imu
          Yimu<-Yimu*(imu/imu.old)
        }
                                          
        imu.vec<-c(imu.vec,imu)
        L.rank.vec<-c(L.rank.vec,if(!is.null(L.svd)) ifelse(L.svd$k.too.small,NA,length(L.svd$d)) else NA)
        #TODO consider setting k basing on L.rank.vec (and imu shrinkage/growth?)
            
        stats<-c(stats,resid.norm)
        if(trace){ 
            #stats.objective.f<-c(stats.objective.f,objective.f.value)
            stats.objective.f <- rbind(stats.objective.f, objective.f.value)
        }
        pfun <- if (i%%message.iter==0) function(...)message(paste(paste(names(...),collapse=", "),"\n",
                                                                   paste(...,collapse=", "))) else print
        if (trace || (i%%message.iter==0)) {
            pfun(c(iter=as.integer(j),k=k.current,L.rank=tail(L.rank.vec,1),
                   resid.norm=resid.norm*(term.delta/term.norm),mu=1/imu,imu=imu))
            pfun("----------------------------------------------------")
        }
        converged<-ifelse(is.na(resid.norm), FALSE, resid.norm<term.norm)
        #if ((i>=local.max.iter)||converged||imu>=(1/mu.min)){ #mu<=mu.min 
        if (converged||imu>=(1/mu.min)){ #mu<=mu.min 
            if (trace) { #DEBUG
                print("[!]Breaking the inner loop:")
                print(c(converged=converged,"mu<=mu.min"=imu>=(1/mu.min)))
            }
            break;        
        }
    } ## end while(TRUE)
    E<-(1-dd)*YYimu
    
    #local.max.iter.phase<-local.max.iter.phase*q #No longer used
    #local.max.iter<-(local.max.iter-i)+local.max.iter.phase    
    k.current=k.current+1        ##TO DO DONE? is this still needed? TOTEST when commented out 
                                 ##Looks like this come into play when (L.svd$k.too.small) above 
                                 ##and  imu<-min(tail(L.svd$d,1),1/mu.min) is set to the 1/mu.min
    } ## end while (j<max.iter && !converged)
    final.delta<-resid.norm*term.delta/term.norm
    if (!converged) {
        warning(paste("rpca did not converge after",i,"iterations, final.delta=",final.delta))        
        if (imu>=(1/mu.min))
            warning(paste("rpca did not converge beacuse given k was apparently too small!"))
    }
    list(L=L,S=S,E=E,L.svd=L.svd,
         convergence=list(converged=converged,
                          iterations=length(stats),
                          final.delta=final.delta,
                          all.delta=stats*(term.delta/term.norm),
                          stats=stats,
                          stats.objective.f=stats.objective.f,
                          mu.vec=1/imu.vec,
                          L.rank.vec=L.rank.vec,
                          mu=mu,
                          mu.max=mu.max,
                          mu.min=mu.min,
                          jumps=jumps))
}
               
trpca.nonsparse<-function(...)trpca(...,
             thresh.nuclear.fun=trpca.thresh.nuclear,
             thresh.l1.fun=thresh.l1,
             zero.matrix.fun=zero.matrix)
                      
                      

#### Previous version 0.2.3 of rpca package
#### Repeating code from above removed

#La.svd.cmp<-cmpfun(La.svd)
thresh.nuclear<-cmpfun(function(M,thr){
    s<-La.svd(M)
    dd<-thr.op(s$d,thr)
    id<-which(dd!=0)
    s$d<-dd[id]
    s$u<-s$u[,id,drop=FALSE]
    s$vt<-s$vt[id,,drop=FALSE]
    s$L<-s$u%*%(s$d*s$vt)
    s
})

rpca<-function(M,
             lambda = 1/sqrt(max(dim(M))),
             mu = prod(dim(M))/(4*sum(abs(M))),
             term.delta=10^(-7),
             max.iter=5000,
             trace=FALSE,
             thresh.nuclear.fun=thresh.nuclear,
             thresh.l1.fun=thresh.l1,
             F2norm.fun=F2norm){
    dm<-dim(M)
    term.norm<-term.delta*F2norm.fun(M)
    S<-Yimu<-matrix(0,nrow=dm[1],ncol=dm[2])

    imu<-1/mu
    limu<-lambda/mu

    i<-0
    stats<-c()
    converged<-FALSE
    while(TRUE){
        i<-i+1
        L.svd<-thresh.nuclear.fun(M-S+Yimu,imu)
        L<-L.svd$L
        S<-thresh.l1.fun(M-L+Yimu,limu)

        MLS = M-L-S

        resid.norm<-F2norm.fun(MLS)
        stats<-c(stats,resid.norm)
        if (trace) 
            print(c(iter=i,resid.norm=resid.norm))
        converged<-resid.norm<term.norm
        if ((i>max.iter)||converged)
            break;

        Yimu = Yimu + MLS
    }
    final.delta<-resid.norm*term.delta/term.norm
    if (!converged)
        warning(paste("rpca did not converge after",i,"iterations, final.delta=",final.delta))
    list(L=L,S=S,L.svd=L.svd,
         convergence=list(converged=converged,
                          iterations=i,
                          final.delta=final.delta,
                          all.delta=stats*(term.delta/term.norm)))
}

# #require(gputools)
# 
# F2norm.gpu<-cmpfun(function(M)sqrt(gputools::gpuCrossprod(as.vector(M))))
# 
# #Below function could be in principle be gpu improved, 
# #but gputools do not expose required simple multiplication functions
# thr.op.gpu<-cmpfun(function(x,thr){sign(x)*pmax(abs(x)-thr,0)})
# 
# thresh.nuclear.gpu<-cmpfun(function(M,thr){
#     s<-gputools::gpuSvd(M)
#     dd<-thr.op.gpu(s$d,thr)
#     id<-which(dd!=0)
#     s$d<-dd[id]
#     s$u<-s$u[,id,drop=FALSE]
#     s$v<-s$v[,id,drop=FALSE]
#     s$vt<-t(s$v)
#     s$L<-gputools::gpuMatMult(s$u,s$d*s$vt)
#     s
# })
# 
# 
# thresh.l1.gpu<-cmpfun(function(M,thr){
#     thr.op.gpu(M,thr)
# })
# 
# rpca.gpu<-function(M,
#              lambda = 1/sqrt(max(dim(M))),
#              mu = prod(dim(M))/(4*sum(abs(M))),
#              term.delta=10^(-5),
#              max.iter=5000,
#              trace=FALSE,
#              gpu.to.choose=NULL){
#     if (!requireNamespace("gputools", quietly = TRUE)) {
#         stop("package 'gputools' version 0.26 is needed for this function to work. Please install it with CULA library: install.packages(\"gputools_0.26.tar.gz\", configure.args=\"--with-cuda-home=/opt/cuda --with-cula-home=/opt/cula\")", call. = FALSE)
#     }
#     try(attachNamespace("gputools"),silent=TRUE)
#     if (!exists("gpuSvd")) 
#         stop("library(gputools) needs to be compiled with CULA library, and its version (0.26) needs to implement gpuSvd: install.packages(\"gputools_0.26.tar.gz\",configure.args=\"--with-cula-home=/opt/cula\")")
#     if (!is.null(gpu.to.choose))
#         gputools::chooseGpu(gpu.to.choose)    
#     rpca(M,lambda=lambda,mu=mu,term.delta=term.delta,max.iter=max.iter,trace=trace,
#          thresh.nuclear.fun=thresh.nuclear.gpu,thresh.l1.fun=thresh.l1.gpu,F2norm.fun=F2norm.gpu)
# }


