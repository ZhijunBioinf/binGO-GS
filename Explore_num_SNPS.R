library(rrBLUP)
library(data.table)

chgo <- as.data.frame(fread('aimSNP.raw',header=T))
chgo <- chgo[,7:ncol(chgo)]
#goname <- read.csv('aimSNPname.csv')
#colnames(chgo) <- goname[,1]
#colne <- read.csv('go_colnames.csv')
#colnames(chgo) <- colne[,1]
#load('R2_s100.RData')
fi <- read.csv('fixedboxFN.csv')
load('cvidFN.RData')
val <- read.csv('pheno31.csv')
Y <- na.omit(val[,7])
snpna <- which(is.na(val[,7]))
if(length(snpna)>0){
  chgo <- chgo[-snpna,]
}

idx <- 1
id <- which(fi[,idx]==1)
#load(paste0('searchFN_',idx,'.RData'))

rrb <- function(snp, Y, trcvid){
  rr <- matrix(NA,length(Y),1)
  for(i in 1:5){
  id <- trcvid[[i]]
  impute1 <- as.matrix(snp[-id,])
  ar_1test <- as.matrix(snp[id,])
  pheno1 <- Y[-id]
  y_anawer <- mixed.solve(pheno1,Z=impute1,K=NULL,SE=FALSE,return.Hinv = FALSE)
  gc()
  rr[id,1] <- ar_1test %*% y_anawer$u + c(y_anawer$beta)
  }
 R <- cor(rr,Y)^2;return(R)
}

data <- chgo[-id,]
y <- Y[-id]
rm(chgo)

Simul <- function(data,Y,i,snp_num,trcvid){
  set.seed(i)
  aa <- sample(1:ncol(data),snp_num)
  a2 <- rrb(data[,aa],Y,trcvid) 
  return(a2)
}

area1 <- seq(1000,10000,1000)
area2 <- seq(12000,80000,2000)
area <- c(area1,area2)
rm(area1)
rm(area2)

ab <- c()
for(i in 1:length(area)){
   rlt <- c()
   for(j in 1:40){
     ac <- Simul(data, y, j, area[i], cvid)
     rlt <- c(rlt,ac)
     }
   ab <- c(ab,mean(rlt))
}

#ab <- c(ab[1:8],ac, ab[9:length(ab)])
save(ab,file=paste0('searchFN_',idx,'.RData',sep=''))