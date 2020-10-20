# K =3
U <- cbind(rs1$score,rs2$score[,3],rs3$score[,3])

U_ORG <- U

###############################################################
########### Run FBIA ##########################################
###############################################################
library(MCMCpack)
library(pscl)
library(plyr)

#U_initial <- read.table("score.dat")
### pre_function ####
Ind <- function(x){
  k1=k2=0;
  for(i in 1:(length(x)-1))
  {
    if(x[i]!=x[i+1]){
      k1=k1+1
    }else{
      k2=k2+1 
    }
  }
  return(c(k1,k2))
}

#### setting ####

M <- 3
N <- p*(p-1)/2
#r <- 2

#### hyperparameter  ###

a1 <- 1
b1 <- 10
a2 <- 1
b2 <- 1

##### S configurations #### 
permu <- function(n){
  S <- NULL;
  for(m in 0:n)
  {
    S <- rbind(S,t(apply(combn(1:n,m=m),2,function(cm) replace(rep(0,n),cm,1))))
  }
  return(S)
}
S <- permu(M)
#U <- U_initial
psi_final_hat <- NULL;
num <- dim(S)[1]

for(kk in 1:N)
{
  psi <- as.numeric(U[kk,-c(1,2)])
  wi <- NULL;
  psi_com_hat <- NULL;
  for(ite in 1:num)
  {
    n0 <- length(which(S[ite,]==0))
    n1 <- length(which(S[ite,]==1))
    k1 <- Ind(S[ite,])[1]  ## change
    k2 <- M-1-k1  ## not change
    psi0 <- psi[S[ite,]==0]
    psi1 <- psi[S[ite,]==1]
    ## combined psi score ##
    psi_com <- NULL;
    psi0_com <- sum(psi0)/sqrt(length(psi0))
    psi1_com <- sum(psi1)/sqrt(length(psi1))
    psi_com[S[ite,]==0] <- psi0_com
    psi_com[S[ite,]==1] <- psi1_com
    psi_com_hat <- rbind(psi_com_hat,psi_com)
    #### calculate Prob ####
    p0 <- 1/sqrt(n0)*(1/sqrt(2*pi))^n0*gamma(n0/2-1/2+a2)*(1/2*sum(psi0^2)-1/(2*n0)*(sum(psi0))^2+b2)^(-(n0/2-1/2+a2))
    p1 <- 1/sqrt(n1)*(1/sqrt(2*pi))^n1*gamma(n1/2-1/2+a2)*(1/2*sum(psi1^2)-1/(2*n1)*(sum(psi1))^2+b2)^(-(n1/2-1/2+a2))
    if(n0==0)
    {
      wi[ite] <- beta(k1+a1,k2+b1)*p1
    }else if(n1==0){
      wi[ite] <- beta(k1+a1,k2+b1)*p0
    }else{
      wi[ite] <- beta(k1+a1,k2+b1)*p0*p1
    }
  } 
  Prob <- wi/sum(wi) 
  psi_final <- t(t(psi_com_hat)%*%Prob)
  psi_final_hat <- rbind(psi_final_hat,psi_final)
}

score_com <- cbind(U[,c(1,2)],psi_final_hat)
score_com <- data.frame(score_com)
#write.table(score_com,"score_joint.dat",row.name=FALSE,col.name=FALSE) # joint psi scores.

### hypothesis test for all psi scores ###
psi_tran <- function(z){
  q<-pnorm(-abs(z), log.p=TRUE)
  q<-q+log(2.0)
  s<-qnorm(q,log.p=TRUE)
  s<-(-1)*s
  return(s)
}

#score_com <- read.table("score_joint.dat")
tran_psi <- sapply(score_com[,-(1:2)],psi_tran)
tran_psi <- cbind(score_com[,1:2],tran_psi)


psi_all <- as.vector(as.matrix(tran_psi[,-c(1,2)]))
num <- 1:length(psi_all)
U <- cbind(num,num,psi_all)
Uu <- cbind(num,num,psi_all)

#write(round(t(U),6), ncol=3, file="com.pcor.est")
U<-U[order(U[,3]), 1:3]
N<-length(U[,1])
ratio<-ceiling(N/100000)
m<-floor(N/ratio)
m0<-N-m*ratio
s<-sample.int(ratio,m,replace=TRUE)
for(i in 1:length(s)) s[i]<-s[i]+(i-1)*ratio
if(m0>0){
  s0<-sample.int(m0,1)+length(s)*ratio
  s<-c(s,s0)
}
Us<-U[s,]


#write(round(t(Us),6), ncol=3, file="com.pcor.est.sub")
#system("gcc pcorsel2.c -lm -o pcorsel2")
#system("./pcorsel2 com.pcor.est.sub")
#aaa <- read.table("aaa.fdr")[,-1]
#write.table(aaa,paste("aaa_all.fdr",sep=""),col.names=FALSE,row.name=FALSE)

aaa <- pcorselR(Us,ALPHA2 = 0.2)

# score_com[,8] <- rowMeans(score_com[,3:7])
# test1 <- score_com[,3:7] - score_com[,8]

adj1 <- matrix(0,p,p)
adj1[lower.tri(adj1)] <- score_com[,3]
adj1[adj1<aaa] <-0
adj1[adj1>=aaa] <-1
adj1[upper.tri(adj1)] <- t(adj1)[upper.tri(t(adj1))]

#write.table(adj1,"adj1.edge",col.names = FALSE, row.names = FALSE)

adj2 <- matrix(0,p,p)
adj2[lower.tri(adj2)] <- score_com[,4]
adj2[adj2<aaa] <-0
adj2[adj2>=aaa] <-1
adj2[upper.tri(adj2)] <- t(adj2)[upper.tri(t(adj2))]

#write.table(adj2,"adj2.edge",col.names = FALSE, row.names = FALSE)

adj3 <- matrix(0,p,p)
adj3[lower.tri(adj3)] <- score_com[,5]
adj3[adj3<aaa] <-0
adj3[adj3>=aaa] <-1
adj3[upper.tri(adj3)] <- t(adj3)[upper.tri(t(adj3))]

#write.table(adj3,"adj3.edge",col.names = FALSE, row.names = FALSE)

# adj4 <- matrix(0,p,p)
# adj4[upper.tri(adj4)] <- score_com[,6]
# adj4[adj4<aaa] <-0
# adj4[adj4>=aaa] <-1
# adj4[lower.tri(adj4)] <- t(adj4)[lower.tri(t(adj4))]
# 
# write.table(adj4,"adj4.edge",col.names = FALSE, row.names = FALSE)
# 
# adj5 <- matrix(0,p,p)
# adj5[upper.tri(adj5)] <- score_com[,7]
# adj5[adj5<aaa] <-0
# adj5[adj5>=aaa] <-1
# adj5[lower.tri(adj5)] <- t(adj5)[lower.tri(t(adj5))]
# 
# write.table(adj5,"adj5.edge",col.names = FALSE, row.names = FALSE)

### ABOVE ACQUIRE PRIOR MATRICES ########################
### LINGAM NEXT #########################################

pcor1 <- adj1

prior <- pcor1
prior[prior==1] <- -1

#mb1 <- huge(gcm1,lambda = 0.3, method = "mb")
#mbpcor1 <- as.matrix(mb1$path[[1]])
sum(A3)
sum(adj1)/2
sum(rs1$Adj)/2
#sum(mbpcor1)
sum(A3==1&adj3==1)
sum(A3==1&rs3$Adj==1)

cadj1 <- matrix(0,p,p)
cadj1[adj1==1|rs1$Adj==1] <- 1


sum(cadj1)/2
sum(A1==1&cadj1==1)

cadj2 <- matrix(0,p,p)
cadj2[adj2==1|rs2$Adj==1] <- 1

sum(cadj2)/2
sum(A2==1&cadj2==1)

cadj3 <- matrix(0,p,p)
cadj3[adj3==1|rs3$Adj==1] <- 1

sum(cadj3)/2
sum(A3==1&cadj3==1)


prior <- cadj3
prior[prior==1] <- -1

writeMat("pcor.mat",Prior = prior)
#writeMat("X1.mat", ts = X1)

library(matlabr)
options(matlab.path = "E:/MATLAB/bin")
have_matlab()

code = c("load('E:/research/my paper/jointLiNGAM/X3.mat')",
         "load('E:/research/my paper/jointLiNGAM/pcor.mat')",
         "K = permu(transpose(ts),'pk', Prior)",
         "save('Kpk.txt','K','-ascii')"
)

res = run_matlab_code(code,add_clear_all = TRUE)

file.exists("Kpk.txt")


K <- as.numeric(read.table("Kpk.txt"))
Pr <- prior[K,K]
Pr[Pr!=0] <-1
Pr[lower.tri(Pr)]<-0

#Porder <- Pr[order(K),order(K)]
ts <- X3
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


#adj_pdl[adj_pdl!=0] <- p.adjust(Pv[Pv!=0])

adj_pdl[adj_pdl>0.05] <- 0
adj_pdl[adj_pdl!=0] <-1

sum(adj_pdl)

## TPR FPR ####
tpr <- sum(A3==1&adj_pdl==1)/sum(A3)
fpr <- sum(A3==0&adj_pdl==0)/sum(A3==0)
fdr <- 1 - sum(A3==1&adj_pdl==1)/sum(adj_pdl==1)
