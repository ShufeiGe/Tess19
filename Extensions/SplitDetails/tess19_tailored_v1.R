# remove(list = ls())
  if(l.max==Inf){
    l.max <- dim(xyz)[1]
  }
  split.seed=123 #set the seed you want to use
   
  cut.type = 1
 
   ALPHA = 1/10000
  
   # l.max=4  # set to 10 for test (max number of cut allowed)

if (split.seed == -1) {
  split.seed = as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31)
  set.seed(split.seed)
}
tau=Inf
sample.type=1
accept.rate.min=0.05
# N=200 #number of particles 
rep.max=1 #number of trees

 
w.normal = NULL
 

 
suppressMessages(library(purrr)    )

cat(sprintf('Loading data.\n'))
 
V.all <-  as.matrix(xyz) #*****
if(l.max==Inf){
  l.max=dim(data)[1]
}

min <- matrix(apply(V.all,2,min),nrow=dim(V.all)[1],ncol=dim(V.all)[2],byrow=TRUE)
max <- matrix(apply(V.all,2,max),nrow=dim(V.all)[1],ncol=dim(V.all)[2],byrow=TRUE)
V.all <- (V.all-min)/(max-min)
rm("min","max")

 
V.label <- group #*****
#id=which(V.label==-1)
#V.label[id]<-2
id.test <- which(is.na(V.label)) #****
V.label[id.test]<- NA #set the unknow labels as NA; otherwise, they will be grouped as a new group.
group <- as.factor(V.label)

group.level <- levels(group)
group.len <- length(group.level)

Cut_plane_standard <- function(V,V.ID,V.label,sample.type,cut.type,accept.rate.min,dist.max=NA,w.normal=NULL){
  d <- dim(V)[2]; #cat(c("dimension:",d,"\n"))
  P.u <- rep(NA,d)
  theta <- rep(NA,(d-1))
  mu <- c(NA)
  len.projection <- c(NA)
  normal.v <- rep(NA,d)
  sample.cut.count <- 1
  V.left.ID <-NULL
  V.right.ID <- NULL
  t.scale.mu <- NULL
  if(is.null(w.normal)){
    w.normal <- rep(1/d,d)
  }
  
  if(sum(!is.na(V.label))>0){
    V.unique <- unique(V[!is.na(V.label),])
    d1 <- dim(V.unique)[1]
    if(is.na(dist.max)){
      if(d1==1){
        dist.max <- 0
      }else{
        dist.max <- max(dist(V.unique))
      }
    }
    #skip the cut if all obs are from the same group(pausing condition)
    # or if labled vertices are same
    skip.index <- (sum(table(V.label)!=0)==1) | dist.max==0
  }else{
    skip.index =1 # skip if all vertices are unlabled
  }
  
  
  if(!skip.index){
    if(sample.type==1){
      # use the largest distance between vertices as upper bound instead of the diameter of the
      # smallest d-sphere, since
      # length(projection.line.segment(Vertices))<= max{dist between Vertices}
      #                                         <= the diameter of the smallest sphere
      # dist.max would be the smallest boundary of the rejection sampling
      run.index <- 1
      if(cut.type==1){
        while(run.index){
          
          
          normal.v <- rnorm(d,sd=0.1)
          normal.v <- w.normal*normal.v
          normal.v <- normal.v/sqrt(sum(normal.v^2))
          
          t.scale <-  V%*%matrix(normal.v,nrow=d,ncol=1)/sum(normal.v^2) #****
          t.scale.train <- t.scale[!is.na(V.label)]
          t.scale.ends <- range(t.scale.train,na.rm = TRUE)
          V.projection.ends <-rbind(t.scale.ends[1]*normal.v,t.scale.ends[2]*normal.v)
          len.projection <- diff(t.scale.ends)*sqrt(sum(normal.v^2))
          
          
          mu <- runif(1,min=0,max=dist.max)
          mu <- mu/len.projection
          if(mu<1){
            t.scale.mu <- (1-mu)*t.scale.ends[1]+mu*t.scale.ends[2]
            #if(any(t.scale<t.scale.mu) && any(t.scale>t.scale.mu)){
            if(any(t.scale.train<=t.scale.mu) && any(t.scale.train>t.scale.mu)){
              run.index <- 0
              P.u <- (1-mu)*V.projection.ends[1,]+mu*V.projection.ends[2,]
              V.left.ID <- which(t.scale<=t.scale.mu)
              V.right.ID <- which(t.scale>t.scale.mu)
            }
          }else{
            sample.cut.count <- sample.cut.count+1
            run.index <- (1/sample.cut.count>accept.rate.min)*1
            skip.index <- 1-run.index
          }
        }
      }else{
        while(run.index){
          
          normal.v<- c(rmultinom(1,size=1,prob=w.normal))
          
          
          t.scale <-  V%*%matrix(normal.v,nrow=d,ncol=1)/sum(normal.v^2)
          t.scale.train <- t.scale[!is.na(V.label)]
          t.scale.ends <- range(t.scale.train,na.rm = TRUE)
          V.projection.ends <-rbind(t.scale.ends[1]*normal.v,t.scale.ends[2]*normal.v)
          len.projection <- diff(t.scale.ends)*sqrt(sum(normal.v^2))
          
          
          mu <- runif(1,min=0,max=dist.max)
          mu <- mu/len.projection
          
          if(mu<1){
            t.scale.mu<-(1-mu)*t.scale.ends[1]+mu*t.scale.ends[2]
            if(any(t.scale.train<=t.scale.mu) && any(t.scale.train>t.scale.mu)){
              run.index <- 0
              P.u <- (1-mu)*V.projection.ends[1,]+mu*V.projection.ends[2,]
              V.left.ID <- which(t.scale<=t.scale.mu)
              V.right.ID <- which(t.scale>t.scale.mu)
            }
            
          }else{
            sample.cut.count <- sample.cut.count+1
            run.index <- (1/sample.cut.count>accept.rate.min)*1
            skip.index <- 1-run.index
          }
        }
      }
    }else{
      run.index <- 1
      while(run.index){
        
        if(cut.type==1){
          normal.v <- rnorm(d,sd=0.1)
          normal.v <- w.normal*normal.v
          normal.v <- normal.v/sqrt(sum(normal.v^2))
          
          
        }else{
          normal.v<- c(rmultinom(1,size=1,prob=w.normal))
          
        }
        
        t.scale <-  V%*%matrix(normal.v,nrow=d,ncol=1)/sum(normal.v^2)
        t.scale.train <- t.scale[!is.na(V.label)]
        t.scale.ends <- range(t.scale.train,na.rm = TRUE)
        V.projection.ends <-rbind(t.scale.ends[1]*normal.v,t.scale.ends[2]*normal.v)
        len.projection <- diff(t.scale.ends)*sqrt(sum(normal.v^2))
        
        mu <- runif(1,0,1)
        t.scale.mu<-(1-mu)*t.scale.ends[1]+mu*t.scale.ends[2]
        
        if(any(t.scale.train<=t.scale.mu) && any(t.scale.train>t.scale.mu)){
          run.index <- 0
          P.u <- (1-mu)*V.projection.ends[1,]+mu*V.projection.ends[2,]
          V.left.ID <- which(t.scale<=t.scale.mu)
          V.right.ID <- which(t.scale>t.scale.mu)
        }else{
          sample.cut.count <- sample.cut.count+1
          run.index <- (sample.cut.count<11)*1
          skip.index <- 1-run.index
        }
        
      }
      
    }
  }
 
  return(list(N_vec= normal.v,Pu=P.u,
              Len.projection=len.projection,dist.max=dist.max,
              sample.cut.count=sample.cut.count,
              skip.index=skip.index,
              V.left.ID=V.left.ID, #***
              V.right.ID=V.right.ID, #***
              t.scale.left = t.scale[V.left.ID], #*** projection
              t.scale.right = t.scale[V.right.ID],# #*** projection
              V.left.trueID=V.ID[V.left.ID], #***
              V.right.trueID=V.ID[V.right.ID], #***
              t.scale.mu = t.scale.mu
              
  ))
}
 

Generative_Process <- function(partition,V.all,V.label,tau,group.level,group.len,sample.type,cut.type,accept.rate.min,w.normal){
  l <- length(partition$Polytopes)
  d <- dim(V.all)[2]
  tau.v <- partition$tau
  result <- c(0,0,0)
  Cut <- list()
  cut.idx <- 0
  skip.all <- 0  #only set it to 1 when  (1) neither cond1 or cond2 is  true OR (2) cost exceeds the budget .
  group.level2 <- as.factor(group.level)
  order <- NULL #****
  if(tau.v[l]>=tau){
    skip.all <- 1
  }else{
    Lambdas <- sapply(partition$Polytopes, "[[", 2)
    
    ##shrink the candidate space, only choose from polytopes whose
    #  (i) Lambda>0
    #       AND
    #  (ii) obs are not from same group (the obs with known labels)
    cond1 <- (Lambdas>0)
    cond2 <- (apply(is.infinite(sapply(partition$Polytopes,"[[",3)/0),2,sum)>1)
    cond12 <- sum(cond1&cond2)
    if(cond12==0){
      skip.all <- 1
    }else{
      sample.space <- which(cond1&cond2)
      if(length(sample.space)==1){
        j <- sample.space
      }else{
        j <- sample(x=c(sample.space),size=1,prob = Lambdas[sample.space])
      }
      order = j #****
      V.temp.ID <- partition$Polytopes[[j]]$V.ID
      V.temp <- V.all[V.temp.ID,]
      V.temp.label <- V.label[V.temp.ID]
      dist.max <- partition$Polytopes[[j]]$dist.max
      Cut <- Cut_plane_standard(V.temp,V.temp.ID,V.temp.label,sample.type,cut.type,accept.rate.min,dist.max,w.normal)
      if(Cut$skip.index!=1){
        index.left <- Cut$V.left.ID
        index.right <- Cut$V.right.ID
        
        if(length(index.left)==1){
          V.temp.left <- t(as.matrix(V.temp[index.left,]))
        }else{
          V.temp.left <- V.temp[index.left,]
        }
        
        if(length(index.right)==1){
          V.temp.right <- t(as.matrix(V.temp[index.right,]))
        }else{
          V.temp.right <- V.temp[index.right,]
        }
        
        V.temp.ID.left <- V.temp.ID[index.left]
        V.temp.ID.right <- V.temp.ID[index.right]
        
        cut.idx <- 1
        l <- l+1
        
        
        dist.max.left <- max(dist.all[V.temp.ID.left,V.temp.ID.left])
        dist.max.right <-  max(dist.all[V.temp.ID.right,V.temp.ID.right])
        Lambda.left <- dist.max.left/2
        Lambda.right <- dist.max.right/2
        
        
        label.est.l <- getmode(V.label[V.temp.ID.left])
        label.est.r <- getmode(V.label[V.temp.ID.right])
        
        countbygroup.left <- as.numeric(table(c(V.label[V.temp.ID.left],group.level2))-1)
        countbygroup.right <- as.numeric(table(c(V.label[V.temp.ID.right],group.level2))-1)
        partition$Polytopes[[j]] <- list(V.ID=V.temp.ID.left,Lambda=Lambda.left,CountByGroup=countbygroup.left,dist.max=dist.max.left,label.est=label.est.l)
        partition$Polytopes[[l]] <- list(V.ID=V.temp.ID.right,Lambda=Lambda.right,CountByGroup=countbygroup.right,dist.max=dist.max.right,label.est=label.est.r)
        
   
        Lambdas <- sapply(partition$Polytopes, "[[", 2)
        tau.v.plus1 <- rexp(n=1,sum(Lambdas))+tau.v[l-1]
        tau.v <- c(tau.v,tau.v.plus1)
        partition$tau <- tau.v
      }
    }
    
  }
  
 
  return(list(partition=partition,
              cut.indx=cut.idx,
              Cut = Cut, 
              skip.all=skip.all))
}


getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


partition_initial <- function(group.level,group.len,V.all,V.label,l.max,tau,sample.type,cut.type,accept.rate.min){
  group.level2 <- as.factor(group.level)
  n <- length(V.label)
  V.temp.ID <- c(1:n)
  V.temp <- V.all[V.temp.ID,]
  d <- dim(V.all)[2]
  Polytopes <- list()
  #countbygroup <- rep(NA,group.len)
  dist.all <- as.matrix(dist(V.all))
  
  
  
  countbygroup <- as.numeric(table(c(V.label,group.level2))-1)

   
  dist.max <- max(dist.all)
  Polytopes[[1]] <- list(V.ID=V.temp.ID,  Lambda=dist.max/2, CountByGroup=countbygroup,dist.max=dist.max,label.est =getmode(V.label) )
  
  
  Partition.Inital <- list(Polytopes=Polytopes,tau=rexp(1,rate=Polytopes[[1]]$Lambda))
  
  output <- list(Partition.Inital=Partition.Inital,dist.all=dist.all)
  return(output)
}


#----inference -----start from here--------

#------------------Partition versus no. cuts (No. of cuts =l.max or tau=Inf)------------------
alpha0 <- matrix((round(table(group)/table(group)[1],2))*ALPHA,ncol=1) #try math paramter
t0 <- proc.time()
table.out <- c()

cat(sprintf("\nEstimated runtime: ?\n"))
pb <- txtProgressBar(min = 0, max = N*rep.max, style = 3)

set.seed(split.seed)

part_init <- partition_initial(group.level,group.len,V.all,V.label,l.max,tau,sample.type,cut.type,accept.rate.min)

Partition.Inital <- part_init$Partition.Inital
dist.all <- part_init$dist.all


Partition.t0 <- list()
Lambda0 <- Partition.Inital$Polytopes[[1]]$Lambda
for(i in 1:N){
  tau.temp <- rexp(1,rate=Lambda0)
  Partition.t0[[i]] <- list(Polytopes=Partition.Inital$Polytopes,tau=tau.temp)
}

cbg0 <- sapply(Partition.Inital$Polytopes, "[[", 3)
multi.beta.0 <- apply(cbg0,2, (function(x){x+alpha0}))
logl.tminus0 <- rep(sum(lgamma(multi.beta.0))-sum(lgamma(apply(multi.beta.0,2,sum)))-dim(cbg0)[2]*(sum(lgamma(alpha0))-lgamma(sum(alpha0))),N)

Wt0 <- rep(1/N,N)   # Weights vector for N particles at time t=0
logWt0 <- log(Wt0)   # log Weights vector for N particles at time t=0

label.predict.all <- NULL
for(rep in 1:rep.max){ 
  print(c("rep is",rep)) 
  setTxtProgressBar(pb, ((rep-1)*N))
  
  Partition.t <- Partition.t0
  Wt <- Wt0        # Weights vector for N particles at time t
  logWt <- logWt0  # log Weights vector for N particles at time t
  logl.tminus1 <- logl.tminus0 #log likelihood for N particles at time t-1
  logl.t <- logl.tminus1 #log likelihood for N particles at time t
  # create progress bar
  tau.li <- Partition.Inital$tau
  
  end.cond <- FALSE
  
  skip.index.N <- rep(0,N)
  
  t <- 1  # here t is the number of cuts; we use vector tau to denote the cost at each cut.
  t.idx <- rep(0,l.max)
  Cut.all <- list()
  for(i in 1:N){
    Cut.all[[i]]<- list()
  }
  
  tt <- 1  #*****
  while(tau.li<tau){
     
    
    #resampling-------------start
    if(t>1){
      
      jC <- sample(x=1:N,size=N,prob =Wt,replace = TRUE )
      
      
      Partition.t.new <- list()
      Cut.all.new <- list()
      for(i in 1:N){
        Partition.t.new[[i]] <- Partition.t[[jC[i]]]
        Cut.all.new[[i]] <- Cut.all[[jC[i]]]
      }
      
      
      Wt[1:N] <- 1/N
      
      
      logWt <- log(Wt)
      logl.t <- logl.t[jC]
      logl.tminus1 <- logl.tminus1[jC]
      skip.index.N <- skip.index.N[jC]
      Partition.t <- Partition.t.new
      Cut.all <- Cut.all.new
      rm("Partition.t.new","Cut.all.new")
    }
    #resampling-------------end
    
    
    for (i in 1:N){
      if(!skip.index.N[i]){
         
        partition.temp <- Partition.t[[i]]
        Out.temp <- Generative_Process(partition.temp,V.all,V.label,tau,group.level,group.len,sample.type,cut.type,accept.rate.min,w.normal)
 
         
        
        if(Out.temp$cut.indx){
          Cut.len.seg <- length(Out.temp$partition$Polytopes)-1
          Cut.all[[i]][[Cut.len.seg]] <- Out.temp$Cut #****
          cbg.new <- sapply(Out.temp$partition$Polytopes, "[[", 3)  
          multi.beta.a <- apply(cbg.new,2, (function(x){x+alpha0}))
          logl.t[i] <- sum(lgamma(multi.beta.a))-sum(lgamma(apply(multi.beta.a,2,sum)))-dim(cbg.new)[2]*(sum(lgamma(alpha0))-lgamma(sum(alpha0)))
          logWt[i] <- logWt[i]+logl.t[i]-logl.tminus1[i]
          logl.tminus1[i] <- logl.t[i]
          Partition.t[[i]] <- Out.temp$partition
        }else{
          skip.index.N[i] <- Out.temp$skip.all
          i.rep <- 0
          if(!skip.index.N[i]){
            while(i.rep<10){
              Out.temp <- Generative_Process(partition.temp,V.all,V.label,tau,group.level,group.len,sample.type,cut.type,accept.rate.min,w.normal)
               
              i.rep <- i.rep+1
              if(Out.temp$cut.indx){
                Cut.len.seg <- length(Out.temp$partition$Polytopes)-1
                Cut.all[[i]][[Cut.len.seg]] <- Out.temp$Cut #****
                cbg.new <- sapply(Out.temp$partition$Polytopes, "[[", 3)  
                multi.beta.a <- apply(cbg.new,2, (function(x){x+alpha0}))
                logl.t[i] <- sum(lgamma(multi.beta.a))-sum(lgamma(apply(multi.beta.a,2,sum)))-dim(cbg.new)[2]*(sum(lgamma(alpha0))-lgamma(sum(alpha0)))
                logWt[i] <- logWt[i]+logl.t[i]-logl.tminus1[i]
                logl.tminus1[i] <- logl.t[i]
                Partition.t[[i]] <- Out.temp$partition
                i.rep <- Inf
              }else{
                skip.index.N[i] <- Out.temp$skip.all
              }
            }
          }
        }
      }
      
    }
    
    
    Wt <- exp((logWt-max(logWt)))/sum(exp((logWt-max(logWt))))
    
    t <- max(length(Partition.t[[which.max(Wt)]]$tau)-1,1)
    
    
    tau.li <- min(map_dbl(Partition.t,~rev(.$tau)[1]))
    end.cond <- min(skip.index.N)
    
    t.idx[t] <- t.idx[t]+1
    
     
    if(end.cond | (t+1)>l.max | tau.li>tau|(t.idx[t]>100) ){
      tau.li <- Inf
    }
    
    if(tau.li==Inf){
      
      count.by.group.temp <- sapply(Partition.t[[which.max(Wt)]]$Polytopes, "[[", 3)
      count.by.group.temp.idx <- apply(count.by.group.temp,2,which.max)
      
      ID.temp <- sapply(Partition.t[[which.max(Wt)]]$Polytopes, "[[", 1)
      label.predic <- rep(NA,length(V.label))
      for(id.j in 1:length(ID.temp)){
        label.predic[ID.temp[[id.j]]] <- as.numeric(group.level[count.by.group.temp.idx[id.j]])
      }
      
      pred.t <- label.predic #**** all prediction
      label.predict.all<- cbind(label.predict.all,pred.t)
      
      Cut.opt <- Cut.all[[which.max(Wt)]]
      Partition.t.opt <- Partition.t[[which.max(Wt)]]
    }
  }
  
  
  if(rep==1){
    zz = round((proc.time()-t0)[3]/60/60*rep.max,4)
    if (zz < 1) {
      zz = round(zz * 60)
      if (zz <= 1) {
        zz = "1 minute"
      } else {
        zz = sprintf("%d minutes", zz)
      }
    } else if (zz > 24) {
      zz = round(zz/24)
      if (zz <= 1) {
        zz = "1 day"
      } else {
        zz = sprintf("%d days", zz)
      }
    } else {
      zz = round(zz)
      if (zz <= 1) {
        zz = "1 hour"
      } else {
        zz = sprintf("%d hours", zz)
      }
    }
    cat(sprintf("\r\033[K\033[1A\r\033[K\r%s\n", paste("Estimated runtime: ", zz, sep="" )))
  }
}
label.predict.all <- as.matrix(t(na.omit(t(label.predict.all))))
label.predict.RF  <- apply(label.predict.all,1,getmode)

cat(sprintf("\r\033[K\033[K\r\nDone.\n"))
