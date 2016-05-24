setwd("~/Documents/TDP-43/Data_files/Bioscreen")
library(genefilter)
library(org.Sc.sgd.db)
library(namespace)
library(GO.db)

#1 READ ALL FILES 
key <- read.csv(file=".../Files/Extras/Bioscreen_key_no_doubles_FEN1_NOT_RAD27.csv",head=T,sep=",")
aggregation <- read.csv(file=".../Files/Extras/aggregation score low and high for Erik.csv",sep=",",head=T)
aggregation <- aggregation[,c(1,4,5)]

#Growth data files. 
c1 <- read.csv(file=".../Files/Raffinose/Normalised_Not_Averaged_raffinose.csv",sep=",",head=T)
c2 <- read.csv(file=".../Files/Raffinose 1.2M KCl/Normalised_Not_Averaged_Raffinose_1.2M_KCl.csv",sep=",",head=T)
c3 <- read.csv(file="...Files/Raffinose 1M NaCl/Normalised_Not_Averaged_Raffionse_1M_NaCl_value_removed_p11.csv",sep=",",head=T)
c4 <- read.csv(file=".../Files/Raffinose+Galactose 1.2M KCl/Normalised_Not_Averaged_Raffinose_1.2M_KCl_TDP43.csv",sep=",",head=T)
c5 <- read.csv(file=".../Files/Raffinose+Galactose 1M NaCl/Normalised_Not_Averaged_Raffinose_Value_selected_1M_NaCl_TDP43.csv",sep=",",head=T)
c6 <- read.csv(file=".../Files/Glucose/Normalised_Not_Averaged_Glucose.csv",sep=",",head=T)
c7 <- read.csv(file=".../Files/Glucose 1.2M KCl/Normalised_Not_Averaged_Glucose_1.2M_KCl.csv",sep=",",head=T)
c8 <- read.csv(file=".../Files/Glucose 1M NaCl/Normalised_Not_Averaged_Glucose_1M_NaCl.csv",sep=",",head=T)

#Select Columns 
a <- 9
b <- 10
c <- 11

#Generate one matrix from all growth data with the columns selected
gts <- matrix(c(c1[,c(a)],c1[,c(b)],c1[,c(c)],
                c2[,c(a)],c2[,c(b)],c2[,c(c)],
                c3[,c(a)],c3[,c(b)],c3[,c(c)],
                c4[,c(a)],c4[,c(b)],c4[,c(c)],
                c5[,c(a)],c5[,c(b)],c5[,c(c)],
                c6[,c(a)],c6[,c(b)],c6[,c(c)],
                c7[,c(a)],c7[,c(b)],c7[,c(c)],
                c8[,c(a)],c8[,c(b)],c8[,c(c)]),ncol=24)

gts <- data.frame(gts)

#For each row of triplicates if the number of NAs is smalles than 2 then replace the NA with the mean of the orther 2 values
for(row in 1:dim(gts)[1]) {
  for(column in 1:dim(gts)[2]) {
    if(is.na(gts[row,column])) {
      if(column==1 | column==2 | column==3) {
        t <- gts[row,c(1,2,3)]
        if(apply(t, MARGIN=1, FUN=function(x) length(x[is.na(x)]) )<2) {
          gts[row,column] <- (rowMeans(x=t,na.rm=T))
        }
      }
      if(column==4 | column==5 | column==6) {
        t <- gts[row,c(4,5,6)]
        if(apply(t, MARGIN=1, FUN=function(x) length(x[is.na(x)]) )<2) {
          gts[row,column] <- (rowMeans(x=t,na.rm=T))
        }
      }
      if(column==7 | column==8 | column==9) {
        t <- gts[row,c(7,8,9)]
        if(apply(t, MARGIN=1, FUN=function(x) length(x[is.na(x)]) )<2) {
          gts[row,column] <- (rowMeans(x=t,na.rm=T))
        }
      }
      if(column==10 | column==11 | column==12) {
        t <- gts[row,c(10,11,12)]
        if(apply(t, MARGIN=1, FUN=function(x) length(x[is.na(x)]) )<2) {
          gts[row,column] <- (rowMeans(x=t,na.rm=T))
        }
      }
      if(column==13 | column==14 | column==15) {
        t <- gts[row,c(13,14,15)]
        if(apply(t, MARGIN=1, FUN=function(x) length(x[is.na(x)]) )<2) {
          gts[row,column] <- (rowMeans(x=t,na.rm=T))
        }
      }
      if(column==16 | column==17 | column==18) {
        t <- gts[row,c(16,17,18)]
        if(apply(t, MARGIN=1, FUN=function(x) length(x[is.na(x)]) )<2) {
          gts[row,column] <- (rowMeans(x=t,na.rm=T))
        }
      }
      if(column==19 | column==20 | column==21) {
        t <- gts[row,c(19,20,21)]
        if(apply(t, MARGIN=1, FUN=function(x) length(x[is.na(x)]) )<2) {
          gts[row,column] <- (rowMeans(x=t,na.rm=T))
        }
      }
      if(column==22 | column==23 | column==24) {
        t <- gts[row,c(22,23,24)]
        if(apply(t, MARGIN=1, FUN=function(x) length(x[is.na(x)]) )<2) {
          gts[row,column] <- (rowMeans(x=t,na.rm=T))
        }
      }
    }
  }
}

gts$ORF <- c1[,c(2)]

#Mark each column by its sample and condition

colnames(gts) <- c("Raffinose_1","Raffinose_2","Raffinose_3","Raffinose_1.2M_KCL_1","Raffinose_1.2M_KCL_2","Raffinose_1.2M_KCL_3",
                   "Raffinose_1M_NaCl_1","Raffinose_1M_NaCl_2","Raffinose_1M_NaCl_3","Raffinose_Galactose_1.2M_KCL_1",
                   "Raffinose_Galactose_1.2M_KCL_2","Raffinose_Galactose_1.2M_KCL_3",
                   "Raffinose_Galactose_1M_NaCl_1","Raffinose_Galactose_1M_NaCl_2","Raffinose_Galactose_1M_NaCl_3",
                   "Glucose_1","Glucose_2","Glucose_3","Glucose_1.2M_KCL_1","Glucose_1.2M_KCL_2","Glucose_1.2M_KCL_3",
                   "Glucose_1M_NaCl_1","Glucose_1M_NaCl_2","Glucose_1M_NaCl_3","Identifier")


gts <- gts[,c(25,16,17,18,1,2,3,19,20,21,4,5,6,22,23,24,7,8,9,10,11,12,13,14,15)]

gts <- merge(gts,key,by.x="Identifier",by.y="number")
gts <- gts[,c(1,26,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)]
wt_gts <- gts

#Make a table where each value is normalized to WT

for(column in 3:dim(gts)[2]) {
  for(row in 1:dim(gts)[1]) {
    wt_gts[row,c(column)] <- gts[row,c(column)]/gts[202,c(column)]
  }
}

controls <- as.matrix(c("OLDWT","EMPTY","GFP","NEG","WT"))

control_gts <- c()

#Remove the controls from the tables

for(row in 1:dim(gts)[1]) {
  for(name in 1:dim(controls)[1]) {
    if(is.na(gts[row,c(2)]==TRUE)) {
      next()
    }
    else if(gts[row,c(2)]==controls[name,c(1)]) {
      control_gts <- rbind(control_gts,gts[row,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26)])
      gts <- gts[-c(row), ]
    }
  } 
}

control_wt_gts <- c()

for(row in 1:dim(wt_gts)[1]) {
  for(name in 1:dim(controls)[1]) {
    if(is.na(wt_gts[row,c(2)]==TRUE)) {
      next()
    }
    else if(wt_gts[row,c(2)]==controls[name,c(1)]) {
      control_wt_gts <- rbind(control_wt_gts,wt_gts[row,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26)])
      wt_gts <- wt_gts[-c(row), ]
    }
  } 
}

#merge the growth data table with the aggregation table. 

gts <- merge(gts,aggregation,by.x="gene.",by.y="GENE",all.x=T)

#Give each gene its ORF name.

Names <- c()
for(row in 1:dim(gts)[1]) {
  ID = as.character(gts[row,c(1)])
  if((substr(ID,1,1)=="Y" & ((substr(ID,7,7)=="C") | (substr(ID,7,7)=="W")))==FALSE) {
    Names <- rbind(Names,get(ID,org.Sc.sgdCOMMON2ORF))
  } else {
    Names <- rbind(Names,ID)
  }
}
Names <- data.frame(Names)

gts$ORF <- Names[,c(1)]

#Use the ORF names to get GO-terms for biological process for each mutant.

f <- c()
for(row in 1:dim(gts)[1]) {
  ID <- as.character(gts[row,c(29)])
  x <- get(ID,org.Sc.sgdGO)
  if(is.na(x)==TRUE){
    f <- rbind(f,"NA")
    next()  
  }
  else if(length(x) > 0) {
    # Try the first one
    got <- x[[1]]
    got[[1]]
    got[[1]]
    got[[1]]
    k <-get(got[[1]],GOTERM)
    f <- rbind(f,slot(k,"Term"))
  }
}
f <- data.frame(f)

gts$"Biological Process" <- f[,c(1)]

#Get the molecular function for each candidate

f <- c()
for(row in 1:dim(gts)[1]) {
  ID <- as.character(gts[row,c(29)])
  x <- get(ID,org.Sc.sgdGO)
  if(is.na(x)==TRUE){
    f <- rbind(f,"NA")
    next()  
  }
  else if(length(x) > 0) {
    for(anno in 1:length(x)){
      if(x[[anno]][["Ontology"]] == "MF") {
        got <- x[[anno]]
        got[[1]]
        got[[1]]
        got[[1]]
        k <- get(got[[1]],GOTERM)
        f <- rbind(f,slot(k,"Term"))
        break()
      }
    }
  }
}
f <- data.frame(f)
gts$"Molecular Function" <- f[,c(1)]

#Split the table into KCl and NaCl.

kcl <- gts[,c(1,29,2,27,28,30,31)]
nacl <- gts[,c(1,29,2,27,28,30,31)]

#Insert the average for each triplicate set.

k <- as.data.frame(gts[,c(2)])
kw <- as.data.frame(wt_gts[,c(1)])
colnames(k) <- c("Identifier")
colnames(kw) <- c("Identifier")
k$"Glucose Average" <- rowMeans(x=gts[,c(3,4,5)],na.rm=T)
k$"Raffinose Average" <- rowMeans(x=gts[,c(6,7,8)],na.rm=T)
k$"Glucose Average + 1.2M KCl" <- rowMeans(x=gts[,c(9,10,11)],na.rm=T)
k$"Raffinose Average + 1.2M KCl" <- rowMeans(x=gts[,c(12,13,14)],na.rm=T)
k$"Raffinose Average + 1.2M KCl + TDP43" <- rowMeans(x=gts[,c(21,22,23)],na.rm=T)
kw$"Raffinose wt_Average + 1.2M KCl + TDP43" <- rowMeans(x=wt_gts[,c(21,22,23)],na.rm=T)
k <- merge(k,kw,by.x="Identifier",by.y="Identifier")
kcl <- merge(kcl,k,by.x="Identifier",by.y="Identifier")

n <- as.data.frame(gts[,c(2)])
nw <- as.data.frame(wt_gts[,c(1)])
colnames(n) <- c("Identifier")
colnames(nw) <- c("Identifier")
n$"Glucose Average" <- rowMeans(x=gts[,c(3,4,5)],na.rm=T)
n$"Raffinose Average" <- rowMeans(x=gts[,c(6,7,8)],na.rm=T)
n$"Glucose Average + 1M NaCl" <- rowMeans(x=gts[,c(15,16,17)],na.rm=T)
n$"Raffinose Average + 1M NaCl" <- rowMeans(x=gts[,c(18,19,20)],na.rm=T)
n$"Raffinose Average + 1M NaCl + TDP43" <- rowMeans(x=gts[,c(24,25,26)],na.rm=T)
nw$"Raffinose wt_Average + 1M NaCl + TDP43" <- rowMeans(x=wt_gts[,c(24,25,26)],na.rm=T)
n <- merge(n,nw,by.x="Identifier",by.y="Identifier")
nacl <- merge(nacl,n,by.x="Identifier",by.y="Identifier")

k <- as.data.frame(c1[,c(2)])
n <- as.data.frame(c1[,c(2)])
colnames(n) <- c("Identifier")
colnames(k) <- c("Identifier")

c1 <- c1[,c(a,b,c)]
c2 <- c2[,c(a,b,c)]
c3 <- c3[,c(a,b,c)]
c4 <- c4[,c(a,b,c)]
c5 <- c5[,c(a,b,c)]
c6 <- c6[,c(a,b,c)]
c7 <- c7[,c(a,b,c)]
c8 <- c8[,c(a,b,c)]

#Add the number of values missing 

k$"Glucose Values missing" <- apply(c6, MARGIN=1, FUN=function(x) length(x[is.na(x)]) )
k$"Raffinose Values missing" <- apply(c1, MARGIN=1, FUN=function(x) length(x[is.na(x)]) )
k$"Glucose 1.2M KCl Values missing" <- apply(c7, MARGIN=1, FUN=function(x) length(x[is.na(x)]) )
k$"Raffinose 1.2M KCl Values missing" <- apply(c2, MARGIN=1, FUN=function(x) length(x[is.na(x)]) )
k$"Raffinose 1.2M KCl TDP43  Values missing" <- apply(c4, MARGIN=1, FUN=function(x) length(x[is.na(x)]) )
n$"Glucose Values missing" <- apply(c6, MARGIN=1, FUN=function(x) length(x[is.na(x)]) )
n$"Raffinose Values missing" <- apply(c1, MARGIN=1, FUN=function(x) length(x[is.na(x)]) )
n$"Glucose 1.2M KCl Values missing" <- apply(c8, MARGIN=1, FUN=function(x) length(x[is.na(x)]) )
n$"Raffinose 1.2M KCl Values missing" <- apply(c3, MARGIN=1, FUN=function(x) length(x[is.na(x)]) )
n$"Raffinose 1.2M KCl TDP43 Values missing" <- apply(c5, MARGIN=1, FUN=function(x) length(x[is.na(x)]) )

nacl <- merge(nacl,n,by.x="Identifier",by.y="Identifier")
kcl <- merge(kcl,k,by.x="Identifier",by.y="Identifier")

#Make a rowttest comparing glucose to raffinose conditions mixing NaCl and KCl. and no stress, separated on Glucose and RAffinose conditions

rg <- rowttests(x=as.matrix(gts[,c(3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)]),factor(c(1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0)))

#DO an FDR correction. 

rg[,c(3)] <- p.adjust(rg[,c(3)],method="fdr",n=length(rg[,c(3)]))

#IF a value was missing for a candidate in NaCL or KCl it would have NA, these are replaced by the data for Glucose/Raffinose without stress

for(row in 1:dim(rg)[1]) {
  if(is.na(rg[row,c(3)])==TRUE) {
    t <- rowttests(x=as.matrix(gts[,c(3,4,5,6,7,8)]),factor(c(1,1,1,0,0,0)))
    t[,c(3)] <- p.adjust(t[,c(3)],method="fdr",n=length(t[,c(3)]))
    rg[row,c(3)] <- t[row,c(3)]
  }
}
rg$Identifier <- gts[,c(2)]

#Test each condition
#KCl to No stress
k <- rowttests(x=as.matrix(gts[,c(3,4,5,6,7,8,9,10,11,12,13,14)]),factor(c(1,1,1,1,1,1,0,0,0,0,0,0)))
k[,c(3)] <- p.adjust(k[,c(3)],method="fdr",n=length(k[,c(3)]))
k$Identifier <- gts[,c(2)]

#NaCl to No stress
n <- rowttests(x=as.matrix(gts[,c(3,4,5,6,7,8,15,16,17,18,19,20)]),factor(c(1,1,1,1,1,1,0,0,0,0,0,0)))
n[,c(3)] <- p.adjust(n[,c(3)],method="fdr",n=length(n[,c(3)]))
n$Identifier <- gts[,c(2)]

#KCl to KCl and TDP-43 expression
tk <- rowttests(x=as.matrix(gts[,c(12,13,14,21,22,23)]),factor(c(1,1,1,0,0,0)))
tk[,c(3)] <- p.adjust(tk[,c(3)],method="fdr",n=length(tk[,c(3)]))
tk$Identifier <- gts[,c(2)]

#NaCl to NaCl and TDP-43 expression
tn <- rowttests(x=as.matrix(gts[,c(18,19,20,24,25,26)]),factor(c(1,1,1,0,0,0)))
tn[,c(3)] <- p.adjust(tn[,c(3)],method="fdr",n=length(tn[,c(3)]))
tn$Identifier <- gts[,c(2)]

#KCl to KCl and TDP-43 expression in the WT normalized set
wt_tk <- rowttests(x=as.matrix(wt_gts[,c(12,13,14,21,22,23)]),factor(c(1,1,1,0,0,0)))
wt_tk[,c(3)] <- p.adjust(wt_tk[,c(3)],method="fdr",n=length(wt_tk[,c(3)]))
wt_tk$Identifier <- wt_gts[,c(1)]

#NaCl to NaCl and TDP-43 expression in the WT normalized set
wt_tn <- rowttests(x=as.matrix(wt_gts[,c(18,19,20,24,25,26)]),factor(c(1,1,1,0,0,0)))
wt_tn[,c(3)] <- p.adjust(wt_tn[,c(3)],method="fdr",n=length(wt_tn[,c(3)]))
wt_tn$Identifier <- wt_gts[,c(1)]

#Reduce each statistic to p-value and identifier.
p <- rg[,c(4,3)]
k <- k[,c(4,3)]
n <- n[,c(4,3)]
tk <- tk[,c(4,3)]
tn <- tn[,c(4,3)]
wt_tk <- wt_tk[,c(4,3)]
wt_tn <- wt_tn[,c(4,3)]

#Merge them into one table
p <- merge(p,k,by.x="Identifier",by.y="Identifier")
p <- merge(p,n,by.x="Identifier",by.y="Identifier")
p <- merge(p,tk,by.x="Identifier",by.y="Identifier")
p <- merge(p,tn,by.x="Identifier",by.y="Identifier")
p <- merge(p,wt_tk,by.x="Identifier",by.y="Identifier")
p <- merge(p,wt_tn,by.x="Identifier",by.y="Identifier")
colnames(p) <- c("Identifier","Carbon Source p.value","KCl p.value","NaCl p.value",
                 "Different to expressing TDP43 KCl","Different to expressing TDP43 NaCl",
                 "Better than WT in TDP43 KCl","Better than WT in TDP43 NaCl")

#Split the table to KCl and NaCl
pk <- p[,c(1,2,3,5,7)]
pn <- p[,c(1,2,4,6,8)]

yldkcl <- merge(kcl,pk,by.x="Identifier",by.y="Identifier")
yldnacl <- merge(nacl,pn,by.x="Identifier",by.y="Identifier")