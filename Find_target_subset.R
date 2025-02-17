library(data.table)
library(rrBLUP)

chgo <- as.data.frame(fread('mafgo_gene.raw',header=T))
chgo <- chgo[,7:ncol(chgo)]
colne <- read.csv('go_colnames.csv',header=F)
colnames(chgo) <- colne[,1]
fi <- read.csv('fixedboxDM.csv')
load('cvidDM.RData')
val <- read.csv('pheno31.csv')
Y <- na.omit(val[,6])
snpna <- which(is.na(val[,6]))
if(length(snpna)>0){
  chgo <- chgo[-snpna,]
}

snp_num <- c(44000, 44000, 42000, 32000, 44000, 32000, 40000, 40000, 36000, 44000, 30000, 36000, 40000, 44000, 44000)
idx <- 9
id <- which(fi[,idx]==1)
setwd(paste0('/us/user/yjs/bqf/study31/subset',idx))
datap <- fread(paste0('/us/user/yjs/bqf/study31/output/lmmGRp_',idx,'.assoc.txt',sep=''))
datap <- as.matrix(datap)
snp1 <- datap[which(as.numeric(datap[,12])<0.01),2]
datap <- as.data.frame(datap)
pp <- as.numeric(datap[,12])

write.csv(snp1,file = 'find_alterFN',quote = F,row.names = F)
system('./plink --file /us/user/yjs/bqf/study31/go_gene --extract find_alterFN --make-bed --out tmpFN')
system('./plink --bfile tmpFN --indep-pairwise 50 5 0.9 --out tmpFN')
req1 <- read.table('tmpFN.prune.in',header = F)
snp1 <- unique(req1[,1])
system('rm -f tmpFN*')
rm(req1)

prange <- which(pp >= 0.01 & pp < 0.1)
prange_site <- datap[prange,c(2,12)]
prange_order <- order(prange_site[,2],decreasing=F)
prange_confirm1 <- prange_site[prange_order,1][1:10000]

prange <- which(pp >= 0.1 & pp < 0.2)
prange_site <- datap[prange,c(2,12)]
prange_order <- order(prange_site[,2],decreasing=F)
prange_confirm2 <- prange_site[prange_order,1][1:10000]

prange <- which(pp >= 0.2 & pp < 0.3)
prange_site <- datap[prange,c(2,12)]
prange_order <- order(prange_site[,2],decreasing=F)
prange_confirm3 <- prange_site[prange_order,1][1:10000]

prange <- which(pp >= 0.3 & pp < 0.4)
prange_site <- datap[prange,c(2,12)]
prange_order <- order(prange_site[,2],decreasing=F)
prange_confirm4 <- prange_site[prange_order,1][1:10000]

prange <- which(pp >= 0.4 & pp < 0.5)
prange_site <- datap[prange,c(2,12)]
prange_order <- order(prange_site[,2],decreasing=F)
prange_confirm5 <- prange_site[prange_order,1][1:10000]

prange <- which(pp >= 0.5 & pp < 0.6)
prange_site <- datap[prange,c(2,12)]
prange_order <- order(prange_site[,2],decreasing=F)
prange_confirm6 <- prange_site[prange_order,1][1:10000]

prange <- which(pp >= 0.6 & pp < 0.7)
prange_site <- datap[prange,c(2,12)]
prange_order <- order(prange_site[,2],decreasing=F)
prange_confirm7 <- prange_site[prange_order,1][1:10000]

prange <- which(pp >= 0.7 & pp < 0.8)
prange_site <- datap[prange,c(2,12)]
prange_order <- order(prange_site[,2],decreasing=F)
prange_confirm8 <- prange_site[prange_order,1][1:10000]

prange <- which(pp >= 0.8 & pp < 0.9)
prange_site <- datap[prange,c(2,12)]
prange_order <- order(prange_site[,2],decreasing=F)
prange_confirm9 <- prange_site[prange_order,1][1:10000]

prange <- which(pp >= 0.9)
prange_site <- datap[prange,c(2,12)]
prange_order <- order(prange_site[,2],decreasing=F)
prange_confirm10 <- prange_site[prange_order,1][1:10000]


rrblup2 <- function(snp, Y, id){
  impute1 <- as.matrix(snp[-id,])
  ar_1test <- as.matrix(snp[id,])
  pheno1 <- Y[-id]
  y_anawer <- mixed.solve(pheno1,Z=impute1,K=NULL,SE=FALSE,return.Hinv = FALSE)
  gc()
  rlt <- ar_1test %*% y_anawer$u + c(y_anawer$beta)
  r <- cor(rlt,Y[id])^2
  return(r)
}

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

snp2 <- c()
for(i in 1:10){snp2[[i]] <- get(paste0('prange_confirm',i))}

m = 1000
snp <- c()
n=0
for(i in 1:10){
  for(j in 1:10){
    n=n+1
    snp[[n]] <- snp2[[j]][seq((i-1)*m+1,i*m)]}}


Find_alter <- function(snp3, snp4, idx, num){ 
  P = 1
  m = 10
  n = 0
  find_alter <- c()                                 
  while(P==1){
    alter_site <- c()
    rlt <- c()
    n = n + 1
    for(i in 1:length(snp3)){
      cat(' \n i:\n',i)
      q1 <- snp3[1:10]
      q1 <- unique(unlist(q1)) 
      snp3[[length(snp3)+1]] <- snp3[[1]]
      snp3[[1]] <- NULL
      data2 <- chgo[-id,q1]
      write.csv(q1,file = 'find_alterFN',quote = F,row.names = F)
      system('./plink --file /us/user/yjs/bqf/study31/go_gene --extract find_alterFN --make-bed --out tmpFN')
      system('./plink --bfile tmpFN --indep-pairwise 50 5 0.9 --out tmpFN')
      # Sys.sleep(8)
      req1 <- read.table('tmpFN.prune.in',header = F)
      req1 <- req1[,1]
      req2 <- unique(c(req1,snp4))
      alter_site[[i]] <- req2
      system('rm -f tmpFN*')
      #data2 <- chgo[-id,req2]
      rlt1 <- rrb(chgo[-id,req2], Y[-id], cvid)
      rlt <- c(rlt,rlt1)
    }
    max_rlt <- which(rlt==max(rlt))
    max_rlt <- max_rlt[1]
    find_alter <- c(find_alter,alter_site[[max_rlt]])
    find_alter <- unique(find_alter)
    cat(' len:\n',length(find_alter))
    snp3[[max_rlt]] <- NULL
    snp4 <- alter_site[[max_rlt]]
    if(length(find_alter) >= num){
      cat('\n Stop while:\n',length(find_alter))
      P <- 0
    }
  }
  return(find_alter)
}

site <- Find_alter(snp3 = snp, snp4 = snp1, idx = idx, num = snp_num[idx])
write.csv(site, file = 'find_alterFN',quote = F,row.names = F)
system('./plink --file /us/user/yjs/bqf/study31/go_gene --extract find_alterFN --make-bed --out tmpFN')
system('./plink --bfile tmpFN --indep-pairwise 100 5 0.9 --out tmpFN')
req <- read.table('tmpFN.prune.in',header = F)
system('rm -f tmpFN*')
cat('\n idx:\n', idx)
rlt1 <- rrblup2(chgo[,site], Y, id)
cat('\n no re1:\n',rlt1)
cat('  len no re1:',length(site))
req <- req[,1]
rlt2 <- rrblup2(chgo[,snp1], Y, id)
cat('\n marker1:\n',rlt2)
cat('  len marker1:',length(snp1))
rlt <- rrblup2(chgo[,req], Y, id)
cat('\n rlt:\n',rlt)
cat('   len rlt:',length(req))
save(req, file=paste0('updata01ReqGR_',idx,'.RData'))
