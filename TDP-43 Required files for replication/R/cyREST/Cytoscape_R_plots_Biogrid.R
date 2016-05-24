setwd(".../cy-rest-R/")
# Basic setup
library(igraph)
library(RJSONIO)
library(httr)
library(RColorBrewer)
library(GO.db)
library(org.Sc.sgd.db)
source('utility/cytoscape_util.R')

# Subgraphs of original graph object written by K. Ono
sendGraph <- function(g, name, collection) {
  cyjs <- toCytoscape(g)
  urlparam = paste("networks?title=", name, "&collection=", collection, sep="")
  network.url = paste(base.url, urlparam, sep="/")
  res <- POST(url=network.url, body=cyjs, encode="json")
  network.suid = unname(fromJSON(rawToChar(res$content)))
  
  # Layout
  apply.layout.url = paste(base.url, "apply/layouts/force-directed", toString(network.suid), sep="/")
  GET(apply.layout.url)
  
}

original.list <- read.csv(file="...Files/Extras/Bioscreen_key_no_doubles_FEN1_NOT_RAD27.csv",sep=",",head=T)

original.list <- as.character(unique(original.list[,c(2)]))
original.list <- original.list[1:192]


for(i in 1:length(original.list)) {
  x <- mget(as.character(original.list[i]),org.Sc.sgdCOMMON2ORF,ifnotfound=NA)
  if(is.na(x)) {
    next()
  }
  else {
    original.list[i]<- x[[1]]
  }
}

a <- list("BRE5","UBP11","YDJ1","KSP1","RIC1","RAD55","TIF4632","HGH1","RAV1","YPT7","GCN2")


orf <- mget(as.character(a),org.Sc.sgdCOMMON2ORF)
c <- c()
for(i in 1:length(a)) {
  x <- a[[i]]
  for(j in 1:length(orf)) {
    if(x==names(orf[j])) {
      c <- rbind(c,data.frame(GENENAME = x, orf= orf[[j]]))
    }
  }
}

#cyREST Yeast Interactome, from YeastNet.com. YeastNet_v3. 
# Used are Two ready made files:
#   1. YeastNet.v3.txt ( The Yeastnet.com full yeast interaction network with LLS scores. )
#   2. yeast_protein_annotation_data.csv  Mine own prepared file with annotations for all yeast genes included in the
#   Yeast interactome. Can be reproduced via another script (insertnamehere)
# Read network from edge list file
yeast.table <- read.table(file=".../Files/Extras/Your_chosen_Biogrid_network in tab2 format", comment.char = "!",sep = "\t", fill=TRUE, stringsAsFactors = FALSE )
yeast.table <- yeast.table[,c(6,7,19,12,13,14,18,20,21,22,23)]
colnames(yeast.table) <- c("gene1","gene2","score",yeast.table[1,c(4)],yeast.table[1,c(5)],
                           yeast.table[1,c(6)],yeast.table[1,c(7)],yeast.table[1,c(8)],yeast.table[1,c(9)],
                           yeast.table[1,c(10)],yeast.table[1,c(11)])
yeast.table <- yeast.table[2:nrow(yeast.table),]

#Fetch all genes interacting with genes in "c"
yeast.table.2 <- c()
for(i in 1:nrow(c)[1]) {
  incomming <- yeast.table[yeast.table$gene1 %in% c[i,c(2)],]
  outgoing <- yeast.table[yeast.table$gene2 %in% c[i,c(2)],]
  x <- rbind(incomming,outgoing)
  x$orf <- c[i,c(2)]
  yeast.table.2 <- rbind(yeast.table.2,x)
}

#Split the table by interaction with orf in "c" 
orf.searched <- split(yeast.table.2,yeast.table.2$orf)

#remove repeated targets in one orf
yeast.table.3 <- c()
for(i in 1:length(orf.searched)) {
  x <- as.character(unique(orf.searched[[i]]$orf))
  for(j in 1:length(orf.searched)) {
    y <- as.character(unique(orf.searched[[j]]$orf))
    if(x==y) {
      next()
    }
    else {
      a <- orf.searched[[i]][orf.searched[[i]]$gene1 %in% orf.searched[[j]]$gene1,]
      b <- orf.searched[[i]][orf.searched[[i]]$gene1 %in% orf.searched[[j]]$gene2,]
      d <- orf.searched[[i]][orf.searched[[i]]$gene2 %in% orf.searched[[j]]$gene1,]
      e <- orf.searched[[i]][orf.searched[[i]]$gene2 %in% orf.searched[[j]]$gene2,]
      z <- rbind(a,b)
      y <- rbind(d,e)
      z <- rbind(z,y)
      z <- unique(z)
      yeast.table.3 <- rbind(yeast.table.3,z)
    }
  }
}

#Unique interactions with LLS interaction > 2 included in the graphing. 
yeast.table.3 <- yeast.table.3[,c(1:11)]
yeast.table.3 <- unique(yeast.table.3)

g <- graph.data.frame(yeast.table.3, directed=T)

# Read Annotation table and add more information about each node in the table
annotations <- read.csv(file="data/yeast_protein_annotation_data.csv",head=T,sep=",")

annotations <- annotations[,c(2,3,4,5,6,7)]

column.names <- colnames(annotations)
filtered <- annotations[1:(length(column.names))]

V(g)$GeneName <- as.character(filtered$GENENAME[match(V(g)$name, filtered$ORF)])
V(g)$Description <- as.character(filtered$DESCRIPTION[match(V(g)$name, filtered$ORF)])
V(g)$KEGGpath <- as.character(filtered$PATH[match(V(g)$name, filtered$ORF)])
V(g)$Entrez <- as.character(filtered$ENTREZID[match(V(g)$name, filtered$ORF)])
V(g)$Chromosome <- as.character(filtered$CHR[match(V(g)$name, filtered$ORF)])

# First: Find which genes in yeast.table.2 that have 2 or more interactions, but remove my genes of interest. 
genes.1 <- as.list(yeast.table.3$gene1)
genes.2 <- as.list(yeast.table.3$gene2)
genes.all <-  unique(c(genes.1, genes.2))
two.interactions.orf <- as.factor(sapply(genes.all, toString))

z <- c()
for(i in 1:length(two.interactions.orf)) {
  x <- yeast.table.3[yeast.table.3$gene1 %in% two.interactions.orf[i],]
  y <- yeast.table.3[yeast.table.3$gene2 %in% two.interactions.orf[i],]
  x <- rbind(x,y)
  genes.1 <- as.list(x$gene1)
  genes.2 <- as.list(x$gene2)
  genes <-  unique(c(genes.1, genes.2))
  if(length(genes)>2) {
    z <- rbind(z,as.data.frame(two.interactions.orf[i]))
  }
}

u <- z[z$`two.interactions.orf[i]` %in% c$orf,]

# Two. Find these gene interactions in yeast.table.2 and make another interaction table. 

yeast.table.4 <- c()
for(i in 1:length(z$`two.interactions.orf[i]`)) {
  if((z[i,c(1)] %in% u)==T) {
    next()
  }
  else {
    x <- yeast.table.3[yeast.table.3$gene1 %in% z[i,c(1)],] 
    y <- yeast.table.3[yeast.table.3$gene2 %in% z[i,c(1)],] 
    x <- rbind(x,y)
    for(j in 1:dim(x)[1]) {
      if(as.character(z[i,c(1)])==x[j,c(1)]) {
        yeast.table.4 <- rbind(yeast.table.4,x)
      }
      else if(as.character(z[i,c(1)])==x[j,c(2)]) {
        yeast.table.4 <- rbind(yeast.table.4,x)
      }
    }
  }
}

################################################################################################
# Include genes with an LLS >3.5 or 4 etc. 

genes.1 <- as.list(yeast.table.4$gene1)
genes.2 <- as.list(yeast.table.4$gene2)
genes.all <-  unique(c(genes.1,genes.2))
orfs <- as.factor(sapply(genes.all, toString))
orfs <- orfs[orfs %in% genes.list]
reduced.graph <- induced.subgraph(g, levels(factor(orfs)))

# Post it to Cytoscape
cyjs <- toCytoscape(reduced.graph)
network.url = paste(base.url, "networks?title=Interactome1&collection=BIOGRID_Network", sep="/")
res <- POST(url=network.url, body=cyjs, encode="json")
network.suid = unname(fromJSON(rawToChar(res$content)))

# Apply style
apply.style.url = paste(base.url, "apply/styles", style.name , toString(network.suid), sep="/")
GET(apply.style.url)


# Tweak Layout parameters
layout.params = list(
  name="unweighted",
  value=TRUE
)

layout.params.url = paste(base.url, "apply/layouts/kamada-kawai/parameters", sep="/")
PUT(layout.params.url, body=toJSON(list(layout.params)), encode = "json")

# Apply layout
params <- paste(toString(network.suid), "?column=community.greedy", sep="")
apply.layout.url = paste(base.url, "apply/layouts/kamada-kawai", params, sep="/")
GET(apply.layout.url)

# Perform Edge Bundling
apply.bundling.url = paste(base.url, "apply/edgebundling", toString(network.suid), sep="/")
GET(apply.bundling.url)

# Toggle graphics details
lod.url = paste(base.url, "ui/lod", sep="/")
PUT(lod.url)