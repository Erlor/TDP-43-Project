setwd(".../cy-rest-R/")
# Basic setup
library(igraph)
library(RJSONIO)
library(httr)
library(RColorBrewer)
library(GO.db)
library(org.Sc.sgd.db)
source('utility/cytoscape_util.R')

# Subgraphs of original graph object. Written by K. Ono
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

original.list <- read.csv(file=".../Files/Extras/Bioscreen_key_no_doubles_FEN1_NOT_RAD27.csv",sep=",",head=T)

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

a <- list("BRE5","UBP11","YDJ1","KSP1","RIC1","RAD55","TIF4632","HGH1","RAV1","YPT7")


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
yeast.table <- read.table(file=".../Files/ExtrasYeastNet.v3.txt", comment.char = "!",sep = "\t", fill=TRUE, stringsAsFactors = FALSE )

colnames(yeast.table) <- c("gene1","gene2","LLS Interaction Score")

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
yeast.table.3 <- yeast.table.3[,c(1,2,3)]
yeast.table.3 <- unique(yeast.table.3)
yeast.table.3 <- yeast.table.3[yeast.table.3$`LLS Interaction Score`>2,]

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

# Create a new style.
# Name of this new style
style.name = "PathVisualization"

# Delete the existing style for fresh start...
style.url = paste(base.url, "styles", sep="/")
style.delete.url = paste(style.url, style.name, sep="/")
DELETE(url=style.delete.url)

# Define default values
def.node.color <- list(
  visualProperty = "NODE_FILL_COLOR",
  value = "#FF0000"
)

def.node.size <- list(
  visualProperty = "NODE_SIZE",
  value = 12
)

def.node.border.width <- list(
  visualProperty = "NODE_BORDER_WIDTH",
  value = 0.5
)

def.edge.width <- list(
  visualProperty = "EDGE_WIDTH",
  value = 2
)

def.edge.color <- list(
  visualProperty = "EDGE_STROKE_UNSELECTED_PAINT",
  value = "#000000"
)

def.edge.target.arrow.color = list(
  visualProperty="EDGE_TARGET_ARROW_UNSELECTED_PAINT",
  value = "#000000"
)

def.edge.transparency = list(
  visualProperty="EDGE_TRANSPARENCY",
  value = 100
)

def.node.transparency = list(
  visualProperty="NODE_TRANSPARENCY",
  value = 100
)

def.node.label.transparency = list(
  visualProperty="NODE_LABEL_TRANSPARENCY",
  value = 100
)

def.node.label.color = list(
  visualProperty="NODE_LABEL_COLOR",
  value="000000"
)

def.node.labelposition <- list(
  visualProperty = "NODE_LABEL_POSITION",
  value = "S,NW,c,7.00,0.00"  
)

def.edge.target.arrow <- list(
  visualProperty="EDGE_TARGET_ARROW_SHAPE",
  value="ARROW"
)

defaults <- list(def.node.color, def.node.size, def.node.labelposition, 
                 def.edge.color, def.node.border.width, def.edge.target.arrow, 
                 def.edge.width, def.edge.target.arrow.color,
                 def.node.transparency, def.node.label.transparency,
                 def.node.label.color,def.edge.transparency)

# Visual Mappings
mappings = list()

# This mapping object is also reusable!
discrete.mappings = list()

node.color = list(
  mappingType="discrete",
  mappingColumn="path",
  mappingColumnType="String",
  visualProperty="NODE_FILL_COLOR",
  map = discrete.mappings
)

node.label.color = list(
  mappingType="discrete",
  mappingColumn="path",
  mappingColumnType="String",
  visualProperty="NODE_LABEL_COLOR",
  map = discrete.mappings
)

edge.color = list(
  mappingType="discrete",
  mappingColumn="path",
  mappingColumnType="String",
  visualProperty="EDGE_STROKE_UNSELECTED_PAINT",
  map = discrete.mappings
)

edge.target.arrow.color = list(
  mappingType="discrete",
  mappingColumn="path",
  mappingColumnType="String",
  visualProperty="EDGE_TARGET_ARROW_UNSELECTED_PAINT",
  map = discrete.mappings
)

node.label = list(
  mappingType="passthrough",
  mappingColumn="GeneName",
  mappingColumnType="String",
  visualProperty="NODE_LABEL"
)

mappings = list(
  node.color, node.label, node.label.color, 
  edge.color, edge.target.arrow.color)

style <- list(title=style.name, defaults = defaults, mappings = mappings)
style.JSON <- toJSON(style)

POST(url=style.url, body=style.JSON, encode = "json")

# Post it to Cytoscape
cyjs <- toCytoscape(g)
network.url = paste(base.url, "networks?title=Interactome&collection=YeastNet_v3", sep="/")
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

#ONE:
#FIND GENES THAT INTERACT WITH MY GENES AND THEN ADD THEM TO THE NETWORK
genes.1 <- as.list(yeast.table.3$gene1)
genes.2 <- as.list(yeast.table.3$gene2)
genes.all <-  unique(c(genes.1, genes.2))
two.interactions.orf <- as.factor(sapply(genes.all, toString))

two.interactions.orf <- two.interactions.orf[(!two.interactions.orf %in% c$orf)]

interactions.orfs <- c()
for(j in 1:length(two.interactions.orf)) {
  interactions <- c()
  interactor <- as.character(two.interactions.orf[[j]])
  interactions.1 <- yeast.table[yeast.table$gene1 %in% interactor,]
  interactions.2 <- yeast.table[yeast.table$gene2 %in% interactor,]
  interactions <- rbind(interactions.1,interactions.2)
  interacting.with.orfs.1 <- interactions[interactions$gene1 %in% c$orf,]
  interacting.with.orfs.2 <- interactions[interactions$gene2 %in% c$orf,]
  interactions.with.orfs <- rbind(interacting.with.orfs.1,interacting.with.orfs.2)
  #Since some of the candidates have alot of interactions, This step reduces the number of interactions 
  if(nrow(interactions.with.orfs) > 50) {
    interactions.with.orfs <- interactions.with.orfs[interactions.with.orfs$`LLS Interaction Score` >3.5,]
  }
  else {
    interactions.with.orfs <- interactions.with.orfs[interactions.with.orfs$`LLS Interaction Score` >2,]
  }
  if(nrow(interactions.with.orfs) > 2) {
    interactions.orfs <- rbind(interactions.orfs,interactions.with.orfs)
  }
}

genes.1 <- as.list(interactions.orfs$gene1)
genes.2 <- as.list(interactions.orfs$gene2)
genes.all <-  unique(c(genes.1,genes.2))
orfs <- as.factor(sapply(genes.all, toString))
reduced.graph <- induced.subgraph(g, levels(factor(orfs)))

# Post to cytoscape using the same subroutine that created my gene based networks. 
sendGraph(reduced.graph, "Two_interactions_list", "YeastNet_v3")
