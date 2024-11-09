rm(list=ls())
gc()
######Package loading######
library(igraph)
library(org.Hs.eg.db)
library(clusterProfiler)
library(KEGGgraph)
library(openxlsx)
library(philentropy)
library(GEOquery)
library(readxl)

############GSE30550###############
GSE30550_series_matrix <- read.delim(".../GSE30550_series_matrix.txt", header=FALSE)
GPL9188 <- read_excel(".../GPL9188.xlsx")
anot <- GPL9188
anot <- as.data.frame(anot)
colnames(anot) <- c("ID","ENTREZID")
anot <- as.data.frame(anot)
samples <- GSE30550_series_matrix[1,]
exdata <- GSE30550_series_matrix[-1,]
colnames(exdata) <- exdata[1,]
exdata <- exdata[-1,]
rownames(exdata) <- exdata$ID_REF
exdata <- exdata[,-1]
esetm = exdata[rownames(exdata) %in% anot$ID, ]
esetm$ENTREZID <- as.character(anot$ENTREZID[match(rownames(esetm), anot$ID)])
esetm<-esetm[!is.na(esetm$ENTREZID),]
esetm<-esetm[!duplicated(esetm$ENTREZID),]
rownames(esetm) <- esetm$ENTREZID
esetm <- esetm[,1:268]
samples <- samples[-1]
circsp<-strsplit(as.character(samples),", ")
samples<-do.call(rbind,circsp)
tts <- cbind(samples,colnames(esetm))
colnames(tts) <- c("subject","hour","GSM")

dsubjects <- c("Subject 01","Subject 05","Subject 06","Subject 07","Subject 08","Subject 10","Subject 12","Subject 13","Subject 15")#???? 8 13ȱʧ
nsubjects <- c("Subject 02","Subject 03","Subject 04","Subject 09","Subject 11","Subject 14","Subject 16","Subject 17")#???? 17ȱʧ
subjects <- c(dsubjects,nsubjects)
hours<-c("Baseline","Hour 00","Hour 005","Hour 012",
         "Hour 021","Hour 029","Hour 036","Hour 045",
         "Hour 053","Hour 060","Hour 069","Hour 077",
         "Hour 084","Hour 093","Hour 101","Hour 108")
allgvs <- esetm
allgvs <- apply(allgvs,1,as.numeric)
allgvs <- scale(allgvs)
allgvs <- t(allgvs)
rownames(allgvs) <- rownames(esetm)
colnames(allgvs) <- colnames(esetm)

tts<-rbind(tts,c("Subject 08","Hour 021","GSM758012"))
tts<-rbind(tts,c("Subject 13","Baseline","GSM758088"))
tts<-rbind(tts,c("Subject 13","Hour 036","GSM758092"))
tts<-rbind(tts,c("Subject 17","Hour 036","GSM758155"))


rankscores <- NULL
for(ii in 1:length(subjects)){
  tt1 <- tts[which(tts[,1] == subjects[[ii]]),]
  samsone <- tt1[,3]
  same <- allgvs[,samsone]
  same <- apply(same,2,as.numeric)
  rownames(same) <- rownames(allgvs)
  same1i <- same[,1:2]
  same1irank <- apply(same1i,2,rank)
  samerbase <- same1irank[,1]
  names(samerbase) <- rownames(same)
  expbase <- same[,1]
  
  samern <- same1irank[,2]
  names(samern) <- rownames(same)
  ndiff <- abs(samern-samerbase)
  names(ndiff) <- rownames(same)
  ndiff <- -sort(-ndiff)
  ndiffgenes <- names(ndiff[1:50])
  
  rankscore <- NULL
  for(i in 1:16){
    samei <- as.vector(t(same[,i]))
    expdiff <- abs(samei-expbase)
    names(expdiff) <- rownames(same)
    sameiranks <- rank(samei)
    names(sameiranks) <- rownames(same)
    rankdiff <- abs(sameiranks-samerbase)
    names(rankdiff) <- rownames(same)
    alldiff <- expdiff*rankdiff
    alldiff <- -sort(-alldiff)
    topdiff <- alldiff[1:50]
    ddiffgenes <- names(topdiff)
    if(length(ddiffgenes)==0){
      score <- 0
    }
    else{
      topdiff <- topdiff[ddiffgenes]
      score <- sum(topdiff)/(length(rankdiff)*length(topdiff))
    }
    rankscore <- c(rankscore,score)
  }
  rankscores <- rbind(rankscores,rankscore)
  cat("\n");
  cat(paste("Done  ",ii,sep=""))
}
rownames(rankscores) <- subjects
colnames(rankscores) <- hours



################simulated dataset##############
library(BB)
qs <- c(-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,-0.03,-0.001,0,0.05,0.1)
scores1000 <- NULL
for(j in 1:1000){
  s <- rnorm(18,0,0.05)
  zvs <- NULL
  for(i in 1:length(qs)){
    q=qs[i]
    model1 <- function(z){
      f <- numeric(length(z))
      f[1] <-(8-4*abs(q)*z[2])/(15+15*z[2])-((4+4*abs(q))/15)*z[1]+s[1]
      f[2] <- (4-2*abs(q)*z[1])/(15+15*z[1])-((8+2*abs(q))*z[2])/(15+15*z[2])+s[2]
      f[3] <- (4*abs(q)-10)/15+(5-2*abs(q))/(15+15*z[1])+(5-2*abs(q))/(15+15*z[2])-z[3]+s[3]
      f[4] <- (12-4*abs(q))/15+(2*abs(q)-6)/(15+15*z[1])+(2*abs(q)-6)/(15+15*z[2])-(6/5)*z[4]+s[4]
      f[5] <- (4*abs(q)-14)/15+(7-2*abs(q))/(15+15*z[1])+(7-2*abs(q))/(15+15*z[2])-(7/5)*z[5]+s[5]
      f[6] <- (4*abs(q)-16)/15+(8-2*abs(q))/(15+15*z[1])+(8-2*abs(q))/(15+15*z[2])-(8/5)*z[6]+s[6]
      f[7] <- (18-4*abs(q))/15+(2*abs(q)-9)/(15+15*z[1])+(2*abs(q)-9)/(15+15*z[2])-(9/5)*z[7]+s[7]
      f[8] <- -(2*z[1])/(15+15*z[1])-(2*z[2])/(15+15*z[2])-(2*z[6])/(5+5*z[6])+(2*z[10])/(5+5*z[10])+(3*z[12])/(5+5*z[12])+(z[15])/(5+5*z[15])-(z[16])/(5+5*z[16])-2*z[8]+s[8]
      f[9] <-  -(z[1])/(5+5*z[1])-(z[2])/(5+5*z[2])-(3*z[6])/(5+5*z[6])-(11/5)*z[9]+s[9]
      f[10] <- (3*z[12])/(5+5*z[12])-(12/5)*z[10]+s[10]
      f[11] <- (z[12])/(4+4*z[1])-(13/5)*z[11]+s[11]
      f[12] <- (2*z[15])/(5+5*z[15])-(2*z[16])/(5+5*z[16])-(14/5)*z[12]+s[12]
      f[13] <- -(z[15])/(5+5*z[15])-(19*z[16])/(5+5*z[16])-5*z[13]+s[13]
      f[14] <- -(4*z[10])/(5+5*z[10])-(4*z[12])/(5+5*z[12])-5*z[14]+s[14]
      f[15] <- (z[16])/(10+10*z[16])-(7/2)*z[15]+s[15]
      f[16] <- (z[15])/(10+10*z[15])-(7/2)*z[16]+s[16]
      f[17] <- -(z[15])/(10+10*z[15])+(z[16])/(10+10*z[16])-(19/5)*z[17]+s[17]
      f[18] <- -(z[15])/(10+10*z[15])+(z[16])/(10+10*z[16])-(z[17])/(5+5*z[17])-4*z[18]+s[18]
      f
    }
    startz <- rep(1,18)
    res <- dfsane(startz,model1,control=list(maxit=2500,trace=FALSE))
    zv <- res$par
    zvs <- cbind(zvs,zv)
  }
  zvsf <- zvs[,1:2]
  zvsf <- apply(zvsf,1,mean)
  zvsfrank <- rank(zvsf)
  scores <- NULL
  for(jj in 1:ncol(zvs)){
    zvsranks <- rank(zvs[,jj])
    rd <- abs(zvsranks-zvsfrank)
    ed <- abs(zvs[,jj]-zvsf)
    score <- sum(rd*ed)
    scores <- c(scores,score)
  }
  scores1000 <- rbind(scores1000,scores)
}
colnames(scores1000) <- qs
write.csv(scores1000,"scores1000.csv")