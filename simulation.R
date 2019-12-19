setwd("E:/lingam")
library(R.matlab)
library(igraph)
library(pcalg)
library(bnlearn)
set.seed(100)
## number of variables i.e. # of nodes
##n = 50, 200,500
##p = 100,200 
n <- 100
p <- 50
degree <-4
gpprob <- degree/(p-1)

tprl <- c()
fprl <- c()
fdrl <- c()
SHDl <- c()
RFl <- c()

for(i in 1:10){
  
  rDAG <- randomDAG(p,prob = gpprob)
  rgDAG <- graph_from_graphnel(rDAG)
  
  #rgDAG <- pgp
  
  ##row out
  ##column in
  A <- as.matrix(get.adjacency(rgDAG))
  sum(A)
  
  mycoord <- layout_with_fr(rgDAG)
  plot.igraph(rgDAG,layout = layout_with_fr,edge.color="orange", vertex.label=V(rgDAG)$number)
  
  # sample <- rmvDAG(n,rDAG, errDist = "cauchy")
  # writeMat("sample.mat",ts = sample)
  
  
  ## generate joint related graph ####
  
  rate <- c(0.05,0.05,0.05)
  
  ## group 1 ###
  A1 <- A
  num1 <- sum(A1)
  B1 <- A1
  B1[B1==1] <- 2*(rbinom(num1,1,0.5)-0.5)* runif(num1, min = 0.3, max = 0.8)
  
  graph1 <- graph_from_adjacency_matrix(A1, mode="directed")
  bn1 <- empty.graph(as.character(1:p))
  amat(bn1) <- A1
  
  ## sample data ####
  ## chisq distribution
  ## exponential distribution
  X1 <- matrix(0,n,p)
  for(i in 1:p){
    # if(permu1[i,2] == 0){
    #   X1[,i] <- rchisq(n,df=1)
    # } else{
    w <- B1[,i]
    X1[,i] <-  X1%*%w + rchisq(n,df=1)
    #X1[,i] <-  X1%*%w + rexp(n,rate =1)
  }
  
  writeMat("X1.mat", ts = X1)
  
  library(huge)
  library(equSA)
  gcm1 <- huge.npn(X1)
  rs1 <- equSAR(gcm1, ALPHA1 =0.05, ALPHA2 = 0.2)
  
  pcor1 <- rs1$Adj
  
  prior <- pcor1
  prior[prior==1] <- -1
  
  #mb1 <- huge(gcm1,lambda = 0.3, method = "mb")
  #mbpcor1 <- as.matrix(mb1$path[[1]])
  sum(A1)
  sum(pcor1)/2
  #sum(mbpcor1)
  sum(A1==1&pcor1==1)
  
  writeMat("pcor.mat",Prior = prior)
  #writeMat("X1.mat", ts = X1)
  
  library(matlabr)
  options(matlab.path = "C:/Program Files/MATLAB/R2016a/bin")
  have_matlab()
  
  code = c("load('E:/lingam/X1.mat')",
           "load('E:/lingam/pcor.mat')",
           "K = permu(transpose(ts),'pk', Prior)",
           "save('Kpk.txt','K','-ascii')"
  )
  
  res = run_matlab_code(code,add_clear_all = TRUE)
  file.exists("Kpk.txt")
  
  
  K <- as.numeric(read.table("Kpk.txt"))
  Pr <- pcor1[K,K]
  Pr[Pr!=0] <-1
  Pr[lower.tri(Pr)]<-0
  
  #Porder <- Pr[order(K),order(K)]
  ts <- X1
  X <- ts[,K]
  B <- matrix(0,p,p)
  Pv <- matrix(0,p,p)
  
  for(i in 1:p){
    if(sum(Pr[,i])>0){
      ind <- which(Pr[,i]==1)
      lm.fit <- lm(X[,i] ~ X[,ind])
      B[ind,i] <- summary(lm.fit)$coefficients[-1,1]
      Pv[ind,i] <- summary(lm.fit)$coefficients[-1,4]
    }
  }
  
  #p.adjust(Pv[upper.tri(Pv)])
  
  B <- B[order(K),order(K)]
  Pv <- Pv[order(K),order(K)]
  adj_pdl <- Pv
  
  
  adj_pdl[adj_pdl!=0] <- p.adjust(Pv[Pv!=0])
  
  adj_pdl[adj_pdl>0.05] <- 0
  adj_pdl[adj_pdl!=0] <-1
  
  
  ## TPR FPR ####
  tpr <- sum(A1==1&adj_pdl==1)/sum(A1)
  fpr <- sum(A1==0&adj_pdl==0)/sum(A1==0)
  fdr <- 1 - sum(A1==1&adj_pdl==1)/sum(adj_pdl==1)
  
  library(bnlearn)
  e <- empty.graph(as.character(1:p))
  amat(e) <- adj_pdl
  
  SHD <- bnlearn::shd(e,bn1) 
  
  RF <- sqrt(sum((B - B1)*(B-B1)))/sqrt(sum(B1*B1))
  
  tprl <- c(tprl,tpr)
  fprl <- c(fprl,fpr)
  fdrl <- c(fdrl,fdr)
  SHDl <- c(SHDl, SHD)
  RFl <- c(RFl,RF)
  
  
}

mean(tprl)
mean(fprl)
mean(fdrl)
mean(SHDl)
mean(RFl)

sd(tprl)
sd(fprl)
sd(fdrl)
sd(SHDl)
sd(RFl)
