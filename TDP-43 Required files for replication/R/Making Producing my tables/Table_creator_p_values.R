#Create the table standard edition
library(genefilter)
key <- read.csv(file="Bioscreen_key_no_doubles.csv",sep=",",head=T)
colnames(p11) <- c("Container.Name","Lag","GT","Yield")
colnames(p12) <- c("Container.Name","Lag","GT","Yield")
colnames(p13) <- c("Container.Name","Lag","GT","Yield")
colnames(p21) <- c("Container.Name","Lag","GT","Yield")
colnames(p22) <- c("Container.Name","Lag","GT","Yield")
colnames(p23) <- c("Container.Name","Lag","GT","Yield")
colnames(p31) <- c("Container.Name","Lag","GT","Yield")
colnames(p32) <- c("Container.Name","Lag","GT","Yield")
colnames(p33) <- c("Container.Name","Lag","GT","Yield")
p111 <- aggregate(p11[,-c(1)],by=list(p11$Container.Name),mean,na.rm=T)
p121 <- aggregate(p12[,-c(1)],by=list(p12$Container.Name),mean,na.rm=T)
p131 <- aggregate(p13[,-c(1)],by=list(p13$Container.Name),mean,na.rm=T)
p212 <- aggregate(p21[,-c(1)],by=list(p21$Container.Name),mean,na.rm=T)
p222 <- aggregate(p22[,-c(1)],by=list(p22$Container.Name),mean,na.rm=T)
p232 <- aggregate(p23[,-c(1)],by=list(p23$Container.Name),mean,na.rm=T)
p313 <- aggregate(p31[,-c(1)],by=list(p31$Container.Name),mean,na.rm=T)
p323 <- aggregate(p32[,-c(1)],by=list(p32$Container.Name),mean,na.rm=T)
p333 <- aggregate(p33[,-c(1)],by=list(p33$Container.Name),mean,na.rm=T)
p1 <- merge(p111,p121,by.x="Group.1",by.y="Group.1",all.x=T,all.y=T)
p1 <- merge(p1,p131,by.x="Group.1",by.y="Group.1",all.x=T,all.y=T)
p2 <- merge(p212,p222,by.x="Group.1",by.y="Group.1",all.x=T,all.y=T)
p2 <- merge(p2,p232,by.x="Group.1",by.y="Group.1",all.x=T,all.y=T)
p3 <- merge(p313,p323,by.x="Group.1",by.y="Group.1",all.x=T,all.y=T)
p3 <- merge(p3,p333,by.x="Group.1",by.y="Group.1",all.x=T,all.y=T)
a <- nrow(p1)
b <- nrow(p2)
c <- nrow(p3)
a <- a-1
b <- b-1
#Divide each row by the NEGative control formt he same plate for normalization
p <- matrix(c(p1[,c(2)]/p1[a,c(2)],p2[,c(2)]/p2[b,c(2)],p3[,c(2)]/p3[c,c(2)],
              p1[,c(5)]/p1[a,c(5)],p2[,c(5)]/p2[b,c(5)],p3[,c(5)]/p3[c,c(5)],
              p1[,c(8)]/p1[a,c(8)],p2[,c(8)]/p2[b,c(8)],p3[,c(8)]/p3[c,c(8)],
              p1[,c(3)]/p1[a,c(3)],p2[,c(3)]/p2[b,c(3)],p3[,c(3)]/p3[c,c(3)],
              p1[,c(6)]/p1[a,c(6)],p2[,c(6)]/p2[b,c(6)],p3[,c(6)]/p3[c,c(6)],
              p1[,c(9)]/p1[a,c(9)],p2[,c(9)]/p2[b,c(9)],p3[,c(9)]/p3[c,c(9)],
              p1[,c(4)]/p1[a,c(4)],p2[,c(4)]/p2[b,c(4)],p3[,c(4)]/p3[c,c(4)],
              p1[,c(7)]/p1[a,c(7)],p2[,c(7)]/p2[b,c(7)],p3[,c(7)]/p3[c,c(7)],
              p1[,c(10)]/p1[a,c(10)],p2[,c(10)]/p2[b,c(10)],p3[,c(10)]/p3[c,c(10)]),ncol=9)
names <- matrix(c(as.character(p1[,c(1)]),as.character(p2[,c(1)]),as.character(p3[,c(1)])),ncol=1)
p1$Graphs <- apply(p1[,c(3,4,6,7,9,10)],MARGIN = 1,FUN=function(x) length(x[!is.na(x)]) )/2
p2$Graphs <- apply(p2[,c(3,4,6,7,9,10)],MARGIN = 1,FUN=function(x) length(x[!is.na(x)]) )/2
p3$Graphs <- apply(p3[,c(3,4,6,7,9,10)],MARGIN = 1,FUN=function(x) length(x[!is.na(x)]) )/2
v <- matrix(c(p1[,c(11)],p2[,c(11)],p3[,c(11)]),ncol=1)
v <- data.frame(v)
p <- data.frame(p)
names <- data.frame(names)
p$ORF <- names[,c(1)]
v$ORF <- names[,c(1)] 
p <- aggregate(p[,-c(10)],by=list(p$ORF),mean,na.rm=T)
v <- aggregate(v[,-c(2)],by=list(v$ORF),sum,na.rm=T)
p <- merge(p,v,by.x="Group.1",by.y="Group.1",all.x=T,all.y=T)
colnames(p) <- c("ORF","Lag_1","Lag_2","Lag_3","Growth rate set 1","Growth rate set 2","Growth rate set 3","Yield set 1","Yield set 2","Yield set 3","Graphs")

variation <- rowVars(x=matrix(c(p1[a,c(3)],p1[a,c(6)],p1[a,c(9)],p2[b,c(3)],p2[b,c(6)],p2[b,c(9)],p3[c,c(3)],p3[c,c(6)],p3[c,c(9)]),ncol=9),na.rm=T)