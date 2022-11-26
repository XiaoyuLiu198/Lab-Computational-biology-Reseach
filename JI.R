#Jaccard Index
library(igraph)
library(dplyr)
library(tidyverse)
library(topGO)
library(GO.db)
library(biomaRt)
library(enrichDO)
library(enrichKEGG)

filenames <- list.files("/Users/molly1998/Desktop/rshinny/datasets", pattern="*.txt")
files<-paste("/Users/molly1998/Desktop/rshinny/datasets", filenames, sep="/")
get_data<-function(path){
  edgelist<-read.table(path,skip = 1, header = F)
  headers<- read.table(path,header = F, nrow = 1, as.is = T)
  colnames(edgelist)= headers
  return(edgelist)
}
edgelists<-lapply(files, get_data)

name_list=c()
i=1
for (name in filenames){
  name_list[i]=sub(".txt", "",toString(name))
  i=i+1
}
names(edgelists) <- name_list

### all data
set=list()
for (i in c(1:length(name_list))){
  before_unique=unlist(c(edgelists[[i]]$regulatoryGene,edgelists[[i]]$targetGene))
  set[[i]]<-unique(before_unique)
}

### calculate intersection
calculate_intersection<-function(data1,data2){
  intersect(data1,data2)
}

### calculate union
calculate_union<-function(data1,data2){
  union(data1,data2)
}

##calculate index
options(digits=5)
for_matrix=c()
for (i in c(1:length(set))){
  for (j in c(1:length(set))){
    union_num=length(calculate_union(set[[i]],set[[j]]))
    intersection_num=length(calculate_intersection(set[[i]],set[[j]]))
    for_matrix[(i-1)*length(set)+j]=intersection_num/union_num
  }
}

##transform to matrix&draw
matrixs=matrix(for_matrix,nrow = length(set))
rownames(matrixs)<-name_list
colnames(matrixs)<-name_list
heatmap(matrixs)

