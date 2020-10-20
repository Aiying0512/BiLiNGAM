setwd("E:/research/my paper/jointLiNGAM")
library(R.matlab)
library(igraph)
library(pcalg)
library(bnlearn)
set.seed(200)
## number of variables i.e. # of nodes
##p = 200
## Total sample size N = 750
##K =   3, 5
##n = 250, 150 
##degree = 1, 2, 5
n <- 250
p <- 200
degree <- 5
gpprob <- degree/(p-1)
#nei <- 3
# ## Generate random graph ###
# rangp <- erdos.renyi.game(p, totedge, type = "gnm", directed = TRUE,
#                  loops = FALSE)
# plot.igraph(rangp,layout = layout_with_kk, vertex.label=V(rangp)$number)



# powergp <- randDAG(p, nei, "er")
# pgp <- graph_from_graphnel(powergp)
# plot.igraph(pgp,layout = layout_with_fr,edge.color="orange", vertex.label=V(pgp)$number)

# weight <- powergp@edgeData
# w <- unlist(weight@data)
# 
# adj <- as.matrix(get.adjacency(pgp))
# W <- adj
# W[W==1] <- w
# W[2,25]

## random graph #####

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

## group 2 ###
A2 <- A1
A2[lower.tri(A2)] <- 2
diag(A2) <- 2

sec0 <- sample(which(A2==1),size=rate[1]*length(which(A2==1)))
sec1 <- sample(which(A2==0),size=length(sec0))

B2 <- B1
B2[sec0] <- 0
B2[sec1] <- runif(length(sec1),min=0.1,max=0.5)

A2 <- B2
A2[A2!=0] <-1

graph2 <- graph_from_adjacency_matrix(B2, mode="directed", weighted = TRUE)

bn2 <- empty.graph(as.character(1:p))
amat(bn2) <- A2
## group 3,4,5 can generate like this as follow
A3 <- A2
A3[lower.tri(A3)] <- 2
diag(A3) <- 2

sec0 <- sample(which(A3==1),size=rate[1]*length(which(A3==1)))
sec1 <- sample(which(A3==0),size=length(sec0))

B3 <- B2
B3[sec0] <- 0
B3[sec1] <- runif(length(sec1),min=0.1,max=0.5)

A3 <- B3
A3[A3!=0] <-1

graph3 <- graph_from_adjacency_matrix(B3, mode="directed", weighted = TRUE)

bn3 <- empty.graph(as.character(1:p))
amat(bn3) <- A3

#test <- A1[order(permu1$node),order(permu1$node)]
##############################################
##permutation before generate data ###########
##############################################
# permu <- data.frame(1:p,colSums(A1))
# permu1 <- permu[order(permu$colSums.A1.),] 
# colnames(permu1) <- c("node","degree")
# A1p <- A1[permu1$node,permu1$node]
# B1p <- B1[permu1$node,permu1$node]
# 
# permu <- data.frame(1:p,colSums(A2))
# permu2 <- permu[order(permu$colSums.A2.),] 
# colnames(permu2) <- c("node","degree")
# 
# A2p <- A2[permu2$node,permu2$node]
# B2p <- B2[permu2$node,permu2$node]


## sample data ####
## chisq distribution
## exponential distribution
X1 <- matrix(0,n,p)
for(i in 1:p){
  # if(permu1[i,2] == 0){
  #   X1[,i] <- rchisq(n,df=1)
  # } else{
    w <- B1[,i]
    X1[,i] <-  X1%*%w + rchisq(n,df=1) - 1
    # X1[,i] <-  X1%*%w + rexp(n,rate =1) 
  }


#X1 <- X1p[,order(permu1$node)]
# 
# colnames(X1) <- permu1$node
#oX1 <- X1[, order(permu1$node)]

X2 <- matrix(0,n,p)

for(i in 1:p){
  # if(permu2[i,2] == 0){
  #   X2[,i] <- rchisq(n,df=1)
  # } else{
    w <- B2[,i]
    X2[,i] <-  X2%*%w + rchisq(n,df=1) - 1
    # X2[,i] <-  X2%*%w + rexp(n,rate =1) 
  }

X3 <- matrix(0,n,p)
for(i in 1:p){
  # if(permu2[i,2] == 0){
  #   X2[,i] <- rchisq(n,df=1)
  # } else{
    w <- B3[,i]
    X3[,i] <-  X3%*%w + rchisq(n,df=1) - 1
  #  X3[,i] <-  X3%*%w + rexp(n,rate =1) 
}

writeMat("X1.mat", ts = X1)
writeMat("X2.mat", ts = X2)
writeMat("X3.mat", ts = X3)
###
library(huge)
library(equSA)
gcm1 <- huge.npn(X1)
gcm2 <- huge.npn(X2)
gcm3 <- huge.npn(X3)
rs1 <- equSAR(gcm1, ALPHA1 =0.05, ALPHA2 = 0.2)
rs2 <- equSAR(gcm2, ALPHA1 =0.05, ALPHA2 = 0.2)
rs3 <- equSAR(gcm3, ALPHA1 =0.05, ALPHA2 = 0.2)

#### RUN FBIA ####

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
options(matlab.path = "E:/MATLAB/bin")
have_matlab()

code = c("load('E:/research/my paper/gemeng/X1.mat')",
         "load('E:/research/my paper/gemeng/pcor.mat')",
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





 
 tpr <- sum(A1==1&adj_ges==1)/sum(A1)
 fpr <- sum(A1==0&adj_ges==0)/sum(A1==0)
 fdr <- 1 - sum(A1==1&adj_ges==1)/sum(adj_ges==1)
 
 tpr <- sum(A1==1&adj_ica==1)/sum(A1)
 fpr <- sum(A1==0&adj_ica==0)/sum(A1==0)
 fdr <- 1 - sum(A1==1&adj_ica==1)/sum(adj_ica==1)
## structural hamming distance ####


#library(bnlearn)
e <- empty.graph(as.character(1:p))
amat(e) <- adj_pdl
amat(e) <- as.matrix(adj_ica)
amat(e) <- adj_ges
amat(e) <- adj_pc
bnlearn::shd(e,bn1) 

sqrt(sum((B - B1)*(B-B1)))/sqrt(sum(B1*B1))


#library(PRROC)
###################################
#################PC################

library(pcalg)
ts <- X1

### pc algorithm ########
gCPDAG <- pc(suffStat = list(C = cor(ts), n = n),
             indepTest = gaussCItest, ## (partial correlations)
             alpha = 1.25*10^(-6), p=p, verbose = FALSE)


#plot(gCPDAG,layout = mycoord)
#gpc <- graph_from_graphnel(gCPDAG)
adj_pc <- showAmat(gCPDAG)
adj_pc[adj_pc==1] <- 0
adj_pc[adj_pc!=0] <- 1
pcgp <- graph_from_adjacency_matrix(adj_pc)

tpr_PC <- sum(A1==1&adj_pc==1)/sum(A1)
fpr_PC <- sum(A1==0&adj_pc==0)/sum(A1==0)
fdr_PC <- 1 - sum(A1==1&adj_pc==1)/sum(adj_pc==1)

e <- empty.graph(as.character(1:p))
amat(e) <- adj_pc
SHD_PC <- bnlearn::shd(e,bn1)


## lingam ICA ####
lingam.fit <- lingam(ts)
g_lingam <- graph_from_adjacency_matrix(t(lingam.fit$Bpruned),mode = "directed",weighted = TRUE,diag = FALSE)
plot.igraph(g_lingam,layout = mycoord,edge.color="orange")

adj_ica <- get.adjacency(g_lingam)

tpr_ICA <- sum(A1==1&adj_ica==1)/sum(A1)
fpr_ICA <- sum(A1==0&adj_ica==0)/sum(A1==0)
fdr_ICA <- 1 - sum(A1==1&adj_ica==1)/sum(adj_ica==1)

e <- empty.graph(as.character(1:p))
amat(e) <- as.matrix(adj_ica)
SHD_ICA <- bnlearn::shd(e,bn1)
