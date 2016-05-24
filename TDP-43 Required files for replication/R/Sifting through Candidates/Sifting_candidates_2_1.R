k <- merge(gtkcl,yldkcl,by.x="Identifier",by.y="Identifier")
n <- merge(gtnacl,yldnacl,by.x="Identifier",by.y="Identifier")

#Remove all values with 0 or 1000 for their aggregation phenotype. in both KCl and NaCl conditions
ckcl <- c()
for(row in 1:dim(k)[1]) {
  if(is.na(k[row,c(5)]) | k[row,c(5)]==0 | k[row,c(5)]==10000) {
    next()
  }
  else {
    ckcl <- rbind(ckcl,k[row, ])
  }
}

cncl <- c()
for(row in 1:dim(n)[1]) {
  if(is.na(n[row,c(5)]) | n[row,c(5)]==0 | n[row,c(5)]==10000) {
    next()
  }
  else {
    cncl <- rbind(cncl,n[row, ])
  }
}

#Move columns in which candidates had no growth at all to a separate table.
k <- c()
not.growing.k <- c()
for(row in 1:dim(ckcl)[1]) {
  if(is.na(ckcl[row,c(40)]) & is.na(ckcl[row,c(41)]) & is.na(ckcl[row,c(42)]) & is.na(ckcl[row,c(43)])) {
    not.growing.k <- rbind(not.growing.k,ckcl[row, ])
  }
  else {
    k <- rbind(k,ckcl[row, ])
  }
}

n <- c()
not.growing.n <- c()
for(row in 1:dim(cncl)[1]) {
  if(is.na(cncl[row,c(40)]) & is.na(cncl[row,c(41)]) & is.na(cncl[row,c(42)]) & is.na(cncl[row,c(43)])) {
    not.growing.n <- rbind(not.growing.n,cncl[row, ])
  }
  else {
    n <- rbind(n,cncl[row, ])
  }
}

#Select candidates on their p-value for growing in TDP-43 and no TDP-43 Wt-normalized and NEG-normalized
ckcl <- c()
for(row in 1:dim(k)[1]) {
  if(!is.na(k[row,c(21)])==TRUE & !is.na(k[row,c(22)])==TRUE & !is.na(k[row,c(42)])==TRUE & !is.na(k[row,c(43)])==TRUE) {
    if(k[row,c(22)]<0.1 & k[row,c(43)]>0.1) {
      ckcl <- rbind(ckcl,k[row, ])
    }
    else if(k[row,c(21)]<0.1 & k[row,c(22)]>0.1 & k[row,c(43)]<0.1) {
      ckcl <- rbind(ckcl,k[row, ])
    }
    else if(k[row,c(21)]<0.1 & k[row,c(22)]>0.1 & k[row,c(42)]<0.1 & k[row,c(43)]>0.1) {
      ckcl <- rbind(ckcl,k[row, ])
    }
    else if(k[row,c(21)]>0.1 & k[row,c(22)]>0.1 & k[row,c(42)]>0.1 & k[row,c(43)]>0.1) {
      ckcl <- rbind(ckcl,k[row, ])
    }
  }
  else if(is.na(k[row,c(42)])==TRUE & is.na(k[row,c(43)])==TRUE & !is.na(k[row,c(21)])==TRUE & !is.na(k[row,c(22)])==TRUE) {
    if(k[row,c(21)]<0.1 & k[row,c(22)]<0.1) {
      ckcl <- rbind(ckcl,k[row, ])
    }
  }
}

cncl <- c()
for(row in 1:dim(n)[1]) {
  if(!is.na(n[row,c(21)])==TRUE & !is.na(n[row,c(22)])==TRUE & !is.na(n[row,c(42)])==TRUE & !is.na(n[row,c(43)])==TRUE) {
    if(n[row,c(22)]<0.1 & n[row,c(43)]>0.1) {
      cncl <- rbind(cncl,n[row, ])
    }
    else if(n[row,c(21)]<0.1 & n[row,c(22)]>0.1 & n[row,c(43)]<0.1) {
      cncl <- rbind(cncl,n[row, ])
    }
    else if(n[row,c(21)]<0.1 & n[row,c(22)]>0.1 & n[row,c(42)]<0.1 & n[row,c(43)]>0.1) {
      cncl <- rbind(cncl,n[row, ])
    }
    else if(n[row,c(21)]>0.1 & n[row,c(22)]>0.1 & n[row,c(42)]>0.1 & n[row,c(43)]>0.1) {
      cncl <- rbind(cncl,n[row, ])
    }
  }
  else if(is.na(n[row,c(42)])==TRUE & is.na(n[row,c(43)])==TRUE & !is.na(n[row,c(21)])==TRUE & !is.na(n[row,c(22)])==TRUE) {
    if(n[row,c(21)]<0.1 & n[row,c(22)]<0.1) {
      cncl <- rbind(cncl,n[row, ])
    }
  }
}

#Select candidates on their p-value comparing Glucose and Raffinose, Glucose conditions can't have a worse growth than raffinose conditions.
k <- c()
for(row in 1:dim(ckcl)[1]) {
  if(!is.na(ckcl[row,c(8)]) & !is.na(ckcl[row,c(9)]) & !is.na(ckcl[row,c(10)]) & !is.na(ckcl[row,c(11)])) {
    if((ckcl[row,c(19)]<0.1 & ckcl[row,c(8)]< ckcl[row,c(9)] & ckcl[row,c(10)]< ckcl[row,c(11)])==FALSE) {
      k <- rbind(k,ckcl[row, ])
    }
  }
}

n <- c()
for(row in 1:dim(cncl)[1]) {
  if(!is.na(cncl[row,c(8)]) & !is.na(cncl[row,c(9)]) & !is.na(cncl[row,c(10)]) & !is.na(cncl[row,c(11)])) {
    if((cncl[row,c(19)]<0.1 & cncl[row,c(8)]< cncl[row,c(9)] & cncl[row,c(10)]< cncl[row,c(11)])==FALSE) {
      n <- rbind(n,cncl[row, ])
    }
  }
}

#Select candidates on their p-value for comparing stress and no stress, stressed values can't be larger than non-stressed.
ckcl <- c()
for(row in 1:dim(k)[1]) {
  if(!is.na(k[row,c(8)]) & !is.na(k[row,c(9)]) & !is.na(k[row,c(10)]) & !is.na(k[row,c(11)])) {
    if((k[row,c(20)]<0.1 & k[row,c(8)]< k[row,c(10)] & k[row,c(9)]< k[row,c(11)])==FALSE) {
      ckcl <- rbind(ckcl,k[row, ])
    }
  }
}

cncl <- c()
for(row in 1:dim(n)[1]) {
  if(!is.na(n[row,c(8)]) & !is.na(n[row,c(9)]) & !is.na(n[row,c(10)]) & !is.na(n[row,c(11)])) {
    if((n[row,c(20)]<0.1 & n[row,c(8)]< n[row,c(10)] & n[row,c(9)]< n[row,c(11)])==FALSE) {
      cncl <- rbind(cncl,n[row, ])
    }
  }
}

#Merge KCl and NaCl data.
ca <- merge(ckcl,cncl,by.x="Identifier",by.y="Identifier")

#Select based on p-values or specific criterias for the mean Generation Time and Mean Yield of a set. 
cand <- c()
for(row in 1:dim(ca)[1]) {
  if(ca[row,13]<1 & ca[row,c(55)]<1 & ca[row,c(34)]>1 & ca[row,c(76)]>1) {
    cand <- rbind(cand,ca[row, ])
  }
  else if(ca[row,13]>1 & ca[row,c(55)]>1 & ca[row,c(34)]<1 & ca[row,c(76)]<1) {
    cand <- rbind(cand,ca[row, ])
  }
  else if(ca[row,21]<0.1 & ca[row,c(63)]<0.1 & ca[row,c(42)]<0.1 & ca[row,c(84)]<0.1) {
    cand <- rbind(cand,ca[row, ])
  }
  else if(ca[row,22]<0.1 & ca[row,c(64)]<0.1 & ca[row,c(43)]>0.1 & ca[row,c(85)]>0.1) {
    cand <- rbind(cand,ca[row, ])
  }
}