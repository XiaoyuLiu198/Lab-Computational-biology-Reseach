#rm(list=ls())
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

name_list=list()
for (name in filenames){
  name_list=append(name_list,sub(".txt", "",toString(name)))
}
names(edgelists) <- name_list

##universe set
all_genes_list=list()
for (name in name_list){
  all_genes_list=append(all_genes_list,unique(unlist(c(edgelists[[name]]$regulatoryGene,edgelists[[name]]$targetGene))))
}
all_genes_list=unique(all_genes_list)

##test gene list
in_name=name_list[1]
edgelist<-edgelists[[in_name[[1]]]]
subset_for_show<-edgelist[order(-edgelist$weight),]
nodelist<-unique(unlist(c(subset_for_show$regulatoryGene,subset_for_show$targetGene)))
graphset<-graph_from_data_frame(d = subset_for_show, vertices = nodelist, directed = TRUE)
method="Hubs & Authorities"
if (method=="Hubs & Authorities"){
  subset_for_tbl<-subset_for_show
  Hub <- hub.score(graphset)$vector
  Authority <- authority.score(graphset)$vector
  li_degree_h=rep(0,length(subset_for_tbl$regulatoryGene))
  li_degree_a=rep(0,length(subset_for_tbl$regulatoryGene))
  for (i in c(1:length(Hub))){
    li_degree_h[i]=Hub[[i]]
    li_degree_a[i]=Authority[[i]]
  }
  subset_for_tbl$centrality_measure<-li_degree_h
  subset_for_tbl$centrality_measure_1<-li_degree_a
  subset_ord<-subset_for_tbl[order(subset_for_tbl$centrality_measure),]
  cut_no=0.3*length(subset_ord$targetGene)
  subset_en<-subset_ord %>% slice_head(n=cut_no)
  gene_list<-unique(unlist(c(subset_en$regulatoryGene,subset_en$targetGene)))
}

##transform
#UNI=keys(org.Hs.eg.db, keytype ='SYMBOL')
db= useMart('ensembl',dataset='hsapiens_gene_ensembl', host="www.ensembl.org")
all_genes_list=as.character(all_genes_list)
go_ids=getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003'), filters='external_gene_name', values=all_genes_list, mart=db)

candidate_list= as.character(gene_list)
# build the gene 2 GO annotation list (needed to create topGO object)
gene_2_GO=unstack(go_ids[,c(1,2)])

# remove any candidate genes without GO annotation
keep = candidate_list %in% go_ids[,2]
keep =which(keep==TRUE)
candidate_list=candidate_list[keep]

# make named factor showing which genes are of interest
geneList=factor(as.integer(all_genes_list %in% candidate_list))
names(geneList)= all_genes_list

GOdata=new('topGOdata', ontology='BP', allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene_2_GO)

classic_fisher_result=runTest(GOdata, algorithm='classic', statistic='fisher')
resultKS = runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim=runTest(GOdata, algorithm = "elim", statistic = "ks")
allRes <- GenTable(GOdata, classicFisher = classic_fisher_result,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "classicFisher", topNodes = 10)


allRes$classicFisher<-as.numeric(allRes$classicFisher)
allRes$classicFisher=-1*log10(allRes$classicFisher)
allRes

library(ensembldb)
library(EnsDb.Hsapiens.v79)
library("AnnotationDbi")
library("org.Hs.eg.db")

cols <- c("SYMBOL","ENTREZID")
Id_universe<-AnnotationDbi::select(org.Hs.eg.db, keys=all_genes_list, columns=cols, keytype="SYMBOL")
Id_candidate<-AnnotationDbi::select(org.Hs.eg.db, keys=candidate_list, columns=cols, keytype="SYMBOL")
Id_candidate<-Id_candidate[order(Id_candidate$ENTREZID,decreasing = TRUE),]

##KEGG
KEGG_res<-enrichKEGG(as.vector(Id_candidate$ENTREZID), organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", as.vector(Id_universe$ENTREZID), minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
KEGG_res<-as.data.frame(KEGG_res)
KEGG_res<-KEGG_res[order(KEGG_res$p.adjust,decreasing = FALSE),]
P1<-ggplot(KEGG_res, aes(x=Description, y=-log10(as.numeric(p.adjust))))+
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("terms") +
  ylab("Ajusted p-value") +
  ggtitle("KEGG Enrichment Analysis") +
  theme_bw(base_size=12) +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=16, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=12, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=12, face="bold", vjust=0.5),
    axis.title=element_text(size=16, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=16),  #Text size
    title=element_text(size=16)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()
P1

##Disease Ontology
DO_res<-enrichDO(as.vector(Id_candidate$ENTREZID), ont = "DO", pvalueCutoff = 0.05, pAdjustMethod = "BH", as.vector(Id_universe$ENTREZID), minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, readable = FALSE)

DO_res<-as.data.frame(DO_res)
DO_res<-DO_res[order(DO_res$p.adjust,decreasing = FALSE),]
P1<-ggplot(DO_res, aes(x=Description, y=-log10(as.numeric(p.adjust))))+
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("terms") +
  ylab("Ajusted p-value") +
  ggtitle("KEGG Enrichment Analysis") +
  theme_bw(base_size=12) +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=16, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=12, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=12, face="bold", vjust=0.5),
    axis.title=element_text(size=16, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=16),  #Text size
    title=element_text(size=16)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()
P1
