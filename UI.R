#rm(list=ls())
library(igraph)
library(dplyr)
library(tidyverse)
library(topGO)
library(GO.db)
library(biomaRt)
library(org.Hs.eg.db)
library(ensembldb)
library(EnsDb.Hsapiens.v79)
library("AnnotationDbi")
library(DOSE)
library(clusterProfiler)
options(scipen = 30)
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
db= useMart('ensembl',dataset='hsapiens_gene_ensembl', host="www.ensembl.org")
all_genes_list=list()
for (name in name_list){
  all_genes_list=append(all_genes_list,unique(unlist(c(edgelists[[name]]$regulatoryGene,edgelists[[name]]$targetGene))))
}
all_genes_list=unique(all_genes_list)
all_genes_list=as.character(all_genes_list)
go_ids=getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003'), filters='external_gene_name', values=all_genes_list, mart=db)

cen_methods<-list("Eigenvector Centrality","Hubs & Authorities","Closeness","Betweenness") #Degree Centrality
enrichment_methods<-list("Disease Ontology","topGO","KEGG") #enrichment analysis methods

##centrality datasets

library(shiny)
#install.packages("DT")
library(shinyWidgets)
library(shinycssloaders)
library(shinythemes)
library(plotly)
library(visNetwork)
library(DT)

UI=shinyUI(
  fluidPage(
    #shinyjs::useShinyjs(),
    id="main",
    title="Network Demo",
    theme = shinytheme("simplex"),
    tags$head(tags$style(
      HTML('
                       body, label, input, button, select { 
                       font-family: "Helvetica","Arial";
                       
                       } 
                       .btn {
                        border-color: grey; 
                        background: grey;
                       }
                       
                      .bttn-bordered.bttn-sm {
                          width: 200px;
                          text-align: left;
                          margin-bottom : 0px;
                          margin-top : 20px;
                       }
                       #refresh{
                       	
                       	margin-top: 40px;
                       }
                       .div{
                       	margin-top : 20px;
                       }
    					table.dataTable th {
    						background-color: #555555 !important;
    						color: white;
    					}
    					#netColumn{
    						border: 1.5px solid #F0F0F0;
    						border-top:0px;
    						margin-left:5px;
    					}
                       '
      )
    )),
    titlePanel("Demo"),
    
    tabsetPanel(
      tabPanel("Network query",  
               icon = icon("object-group"),
               fluidRow(
                 column(2,
                        textInput("gene", label = strong("Enter Gene List"),value="MEF2C,CUX2"),
                        sliderInput(inputId = "threshold", label = strong("Percentage of Nodes"),
                                    min = 1, max = 100, value = 50, step = 1),
                        selectInput(inputId = "cell", label = strong("Cell Type"),choices = c(as.character(name_list))),
                        downloadButton('download_plot', 'Download Network Plot') ),#end column
                 column(9, 
                        visNetworkOutput(outputId = "network_plot",width = "100%", height = "250px")%>% withSpinner(color="black", type=6, size=1)
                          ) 
                        ),#fluidRow
                 #column(1),
               fluidRow(	
                 div(id="table",hr(),
                       column(2,
                              downloadButton('download_tbl', 'Download Table')    
                       ),
                     column(9, dataTableOutput("gene_tbl")	        
                     )
                 )	
               )#fluidRow
      ),#end tab panel
      
      tabPanel("Centrality analysis",  
               icon = icon("object-group"),
               fluidRow(
                 column(2,
                        selectInput(inputId = "cell_enr", label = strong("Cell Type"),choices = c(as.character(name_list))),
                        selectInput(inputId = "central_method", label = strong("Centrality Measure Method"),choices = c(as.character(cen_methods))),
                        selectInput(inputId = "enrichment_method", label = strong("Enrichment Analysis Method"),choices = c(as.character(enrichment_methods))),
                        downloadButton('download_barplot', 'Download Enrichment Plot') ),#end column
                 column(9,
                        textOutput(outputId = "threshold_or_not"),
                        sliderInput(inputId = "enrichment_threshold", label = strong("Percentage of genes in Enrichment"),
                                    min = 1, max = 100, value = 50, step = 1),
                        plotOutput(outputId = "enrich_barplot",width = "100%", height = "250px")
                 ) 
               ),#fluidRow
               #column(1),
               fluidRow(	
                 div(id="table",hr(),
                     column(2,
                            downloadButton('download_enr_tbl', 'Download Enrichment Table')    
                     ),
                     column(9, dataTableOutput("enrich_tbl")	        
                     )
                 )	
               )#fluidRow
      ),#end tab panel
      tabPanel("Logic gate analysis", icon = icon("object-group")
               #includeMarkdown("about.Rmd")
      ),
      tabPanel("About", icon = icon("info")
               #includeMarkdown("about.Rmd")
      ),
      tabPanel("Help", icon = icon("book"))
               #includeMarkdown("docs.Rmd")
      )#end of tabset panel  
    ) )#end of shinyUI	

server = function(input, output) { 	

    #nodelist<-unique(unlist(c(subset_for_show$regulatoryGene,subset_for_show$targetGene)))
    #graphset<-graph_from_data_frame(d = subset_for_show, vertices = nodelist, directed = TRUE)
  data_tbl<-reactive({
    gene_text=as.vector(scan(text = input$gene, sep = ",", what = ""))
    in_name<-input$cell
    threshold<-input$threshold
    edgelist<-edgelists[[in_name[[1]]]]
    edgelist<-edgelist[order(-edgelist$weight),]
    subset<-edgelist %>% filter_all(any_vars(. %in% gene_text))
    thresholdn=floor(length(subset$regulatoryGene)*(threshold/1000))
    subset_for_show<-subset %>% slice_head(n=thresholdn)
    nodeset<-unique(unlist(c(subset_for_show$regulatoryGene,subset_for_show$targetGene)))
    graphset<-graph_from_data_frame(d = subset_for_show, vertices = nodeset, directed = TRUE)
    #if (method=="Degree Centrality"){
    #  subset_for_tbl<-subset_for_show
    #  Degree <- degree(graphset)
    #  li_degree=rep(0,length(subset_for_show$regulatoryGene))
    #  for (i in c(1:length(Degree))){
    #    li_degree[i]=Degree[[i]]
    #  }
    #  subset_for_tbl$centrality_measure_degree<-li_degree
   # }
    if (method=="Eigenvector Centrality"){
      subset_for_tbl<-subset_for_show
      Eig <- evcent(graphset)$vector
      li_degree=rep(0,length(subset_for_show$regulatoryGene))
      for (i in c(1:length(Eig))){
        li_degree[i]=Eig[[i]]
      }
      subset_for_tbl$centrality_measure_eigenvector<-li_degree
    }
    if (method=="Hubs & Authorities"){
      subset_for_tbl<-subset_for_show
      Hub <- hub.score(graphset)$vector
      Authority <- authority.score(graphset)$vector
      li_degree_h=rep(0,length(subset_for_show$regulatoryGene))
      li_degree_a=rep(0,length(subset_for_show$regulatoryGene))
      for (i in c(1:length(Hub))){
        li_degree_h[i]=Hub[[i]]
        li_degree_a[i]=Authority[[i]]
      }
      subset_for_tbl$centrality_measure_hub<-li_degree_h
      subset_for_tbl$centrality_measure_authority<-li_degree_a
    }
    if (method=="Closeness"){
      subset_for_tbl<-subset_for_show
      Closeness <- closeness(graphset)
      li_degree=rep(0,length(subset_for_show$regulatoryGene))
      for (i in c(1:length(Closeness))){
        li_degree[i]=Closeness[[i]]
      }
      subset_for_tbl$centrality_measure_closeness<-li_degree
    }
    if (method=="Betweenness"){
      subset_for_tbl<-subset_for_show
      Betweenness <- betweenness(graphset)
      li_degree=rep(0,length(subset_for_show$regulatoryGene))
      for (i in c(1:length(Betweenness))){
        li_degree[i]=Betweenness[[i]]
      }
      subset_for_tbl$centrality_measure_betweenness<-li_degree
    }
    subset_for_tbl
  })
  
  graph_object<-reactive({
    gene_text=as.vector(scan(text = input$gene, sep = ",", what = ""))
    in_name<-input$cell
    threshold<-input$threshold
    edgelist<-edgelists[[in_name[[1]]]]
    edgelist<-edgelist[order(-edgelist$weight),]
    subset<-edgelist %>% filter_all(any_vars(. %in% gene_text))
    thresholdn=floor(length(subset$regulatoryGene)*(threshold/1000))
    subset_for_show<-subset %>% slice_head(n=thresholdn)
    nodeset<-unique(unlist(c(subset_for_show$regulatoryGene,subset_for_show$targetGene))) 
    graphset<-graph_from_data_frame(d = subset_for_show, vertices = nodeset, directed = TRUE)
    V(graphset)$name<- nodeset
    V(graphset)$color  <- ifelse(V(graphset)$name %in% gene_text, "grey","green")
    g1<-plot(graphset,edge.arrow.size = 0.2,edge.color ="grey",vertex.frame.color="grey",vertex.size=8,vertex.label.cex=0.3,vertex.label.color="black",vertex.label.font=2,vertex.label.family="Times",vertex.color=V(graphset)$color)
    g1
  }) 
  
  output$network_plot <- renderVisNetwork({
    gene_text=as.vector(scan(text = input$gene, sep = ",", what = ""))
    in_name<-input$cell
    threshold<-input$threshold
    edgelist<-edgelists[[in_name[[1]]]]
    edgelist<-edgelist[order(-edgelist$weight),]
    subset<-edgelist %>% filter_all(any_vars(. %in% gene_text))
    thresholdn=floor(length(subset$regulatoryGene)*(threshold/1000))
    subset_for_show<-subset %>% slice_head(n=thresholdn)
    nodeset<-unique(unlist(c(subset_for_show$regulatoryGene,subset_for_show$targetGene))) 
    #graphset<-graph_from_data_frame(d = subset_for_show, vertices = nodeset, directed = TRUE)
    #plot(graphset, edge.arrow.size = 0.2,edge.color ="grey",vertex.size=6,vertex.label.cex=0.3,vertex.label.color="black",vertex.label.font=2,vertex.label.family="Times",vertex.color="grey")
    # visNetwork(nodeset, edges)
    edges<-as.data.frame(cbind(subset_for_show$regulatoryGene,subset_for_show$targetGene))
    nodes<- as.data.frame(cbind(nodeset,nodeset,'circle'))
    colnames(nodes)=c("id","label","shape")
    colnames(edges)=c("from","to")
    vis.links=edges#[,c("V1","V2")]
    nodes$color <- ifelse(nodes$id %in% gene_text,
                          "slategrey","green")
    visnet<-visNetwork(nodes, vis.links) %>% visEdges(arrows = 'to')
    visnet
    })
  
  output$gene_tbl <- renderDataTable({
    gene_text=as.vector(scan(text = input$gene, sep = ",", what = ""))
    in_name<-input$cell
    threshold<-input$threshold
    method<-input$central_method
    edgelist<-edgelists[[in_name[[1]]]]
    edgelist<-edgelist[order(-edgelist$weight),]
    subset<-edgelist %>% filter_all(any_vars(. %in% gene_text))
    thresholdn=floor(length(subset$regulatoryGene)*(threshold/1000))
    subset_for_show<-subset %>% slice_head(n=thresholdn)
    
    df <- datatable(subset_for_show,rownames = FALSE, options = list(pageLength = 15)) %>%
      # First Column Border Left
      formatStyle(c(1),`border-left` = '1px solid black') #%>% 
    # Rest Alternative Bordering
    #formatStyle(alt_vector,`border-right` = '1px solid black')
    df
  })
  
  enr_object<-reactive({
    in_name<-input$cell_enr
    method<-input$central_method
    edgelist<-edgelists[[in_name[[1]]]]
    enrichment_method<-input$enrichment_method
    subset_for_show<-edgelist[order(-edgelist$weight),]
    nodeset<-unique(unlist(c(subset_for_show$regulatoryGene,subset_for_show$targetGene))) 
    graphset<-graph_from_data_frame(d = subset_for_show, vertices = nodeset, directed = TRUE)
   if (method=="Eigenvector Centrality"){
    subset_for_tbl<-subset_for_show
    Eig <- evcent(graphset)$vector
    li_degree=rep(0,length(subset_for_show$regulatoryGene))
    for (i in c(1:length(Eig))){
      li_degree[i]=Eig[[i]]
    }
    subset_for_tbl$centrality_measure<-li_degree
    subset_ord<-subset_for_tbl[order(subset_for_tbl$centrality_measure),]
    cut_no=0.3*length(subset_ord$targetGene)
    subset_en<-subset_ord %>% slice_head(n=cut_no)
    gene_list<-unique(unlist(c(subset_en$regulatoryGene,subset_en$targetGene)))
  }
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
   if (method=="Closeness"){
    subset_for_tbl<-subset_for_show
    Closeness <- closeness(graphset)
    li_degree=rep(0,length(subset_for_tbl$regulatoryGene))
    for (i in c(1:length(Closeness))){
      li_degree[i]=Closeness[[i]]
    }
    subset_for_tbl$centrality_measure<-li_degree
    subset_ord<-subset_for_tbl[order(subset_for_tbl$centrality_measure),]
    cut_no=0.3*length(subset_ord$targetGene)
    subset_en<-subset_ord %>% slice_head(n=cut_no)
    gene_list<-unique(unlist(c(subset_en$regulatoryGene,subset_en$targetGene)))
  }
   if (method=="Betweenness"){
    subset_for_tbl<-subset_for_show
    Betweenness <- betweenness(graphset)
    li_degree=rep(0,length(subset_for_show$regulatoryGene))
    for (i in c(1:length(Betweenness))){
      li_degree[i]=Betweenness[[i]]
    }
    subset_for_tbl$centrality_measure<-li_degree
    subset_ord<-subset_for_tbl[order(subset_for_tbl$centrality_measure),]
    cut_no=0.3*length(subset_ord$targetGene)
    subset_en<-subset_ord %>% slice_head(n=cut_no)
    gene_list<-unique(unlist(c(subset_en$regulatoryGene,subset_en$targetGene)))
   }
  candidate_list= as.character(gene_list)
  if (enrichment_method=="topGO"){
  gene_2_GO=unstack(go_ids[,c(1,2)])
  keep = candidate_list %in% go_ids[,2]
  keep =which(keep==TRUE)
  candidate_list=candidate_list[keep]
  geneList=factor(as.integer(all_genes_list %in% candidate_list))
  names(geneList)= all_genes_list
  GOdata=new('topGOdata', ontology='BP', allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene_2_GO)
  classic_fisher_result=runTest(GOdata, algorithm='classic', statistic='fisher')
  resultKS = runTest(GOdata, algorithm = "classic", statistic = "ks")
  allRes <- GenTable(GOdata, classicFisher = classic_fisher_result,
                     classicKS=resultKS, orderBy = "classicFisher", topNodes = 10)
  }
  if(enrichment_method=="KEGG"){
    threshold_enr<-input$enrichment_threshold
    thresholdn_enr=floor(length(candidate_list)*(threshold_enr/100))
    cols <- c("SYMBOL","ENTREZID")
    Id_universe<-AnnotationDbi::select(org.Hs.eg.db, keys=all_genes_list, columns=cols, keytype="SYMBOL")
    Id_candidate<-AnnotationDbi::select(org.Hs.eg.db, keys=candidate_list, columns=cols, keytype="SYMBOL")
    #Id_candidate<-Id_candidate[order(Id_candidate$ENTREZID,decreasing = TRUE),]
    Id_candidate_after_cut<-Id_candidate %>% slice_head(n=thresholdn_enr)
    allRes<-enrichKEGG(as.vector(Id_candidate_after_cut$ENTREZID), organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, 
                         pAdjustMethod = "BH", as.vector(Id_universe$ENTREZID), minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
  }
  if (enrichment_method=="Disease Ontology"){
    threshold_enr<-input$enrichment_threshold
    thresholdn_enr=floor(length(candidate_list)*(threshold_enr/100))
    cols <- c("SYMBOL","ENTREZID")
    Id_universe<-AnnotationDbi::select(org.Hs.eg.db, keys=all_genes_list, columns=cols, keytype="SYMBOL")
    Id_candidate<-AnnotationDbi::select(org.Hs.eg.db, keys=candidate_list, columns=cols, keytype="SYMBOL")
    #Id_candidate<-Id_candidate[order(Id_candidate$ENTREZID,decreasing = TRUE),]
    Id_candidate_after_cut<-Id_candidate %>% slice_head(n=thresholdn_enr)
    allRes<-enrichDO(as.vector(Id_candidate_after_cut$ENTREZID), ont = "DO", pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                     as.vector(Id_universe$ENTREZID), minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, readable = FALSE)
    
  }
  allRes1<-as.data.frame(allRes)
  allRes1
  })
  
  output$enrich_tbl <- renderDataTable({
    df <- datatable(enr_object(),rownames = FALSE, options = list(pageLength = 15)) %>%
      # First Column Border Left
      formatStyle(c(1),`border-left` = '1px solid black') #%>% 
    # Rest Alternative Bordering
    #formatStyle(alt_vector,`border-right` = '1px solid black')
    df
  })
  
  barplot_object<-reactive({
    allRes<-enr_object()
    enrichment_method<-input$enrichment_method
    if (enrichment_method=="topGO"){
    P1<-ggplot(allRes, aes(x=Term, y=-log10(as.numeric(classicFisher))))+
      stat_summary(geom = "bar", fun = mean, position = "dodge") +
      xlab("terms") +
      ylab("Enrichment Significant") +
      ggtitle("topGO Enrichment Analysis") +
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
    }
    if(enrichment_method=="KEGG" || enrichment_method=="Disease Ontology"){
      allRes<-allRes[order(allRes$p.adjust,decreasing = FALSE),]
      P1<-ggplot(allRes, aes(x=Description, y=-log10(as.numeric(p.adjust))))+
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
    }
    P1
  })
  
  output$enrich_barplot <- renderPlot({
    barplot_object()
  })
  
  output$threshold_or_not=renderText({
    enrichment_method<-input$enrichment_method
    if (enrichment_method=="topGO"){
      print_content<-"No threshold needed.Include all genes in this cell type."
    }
    else{
      print_content<-"Please choose percentage of genes to include in enrichment analysis."
    }
    print_content
    })
  
  output$download_tbl=downloadHandler(
    filename = function() {
      paste("Gene_network_table",input$gene,"csv", sep = ".")
    },
    content = function(file) {
      write.csv(data_tbl(), file, row.names = FALSE, col.names=TRUE, quote=FALSE, sep="\t")
    })
  
  output$download_plot=downloadHandler(
    filename = function() {
      paste("Gene_network_plot",input$gene,"png", sep = ".")
    },
    content = function(file) {
      device <- function(..., width, height) grDevices::png(..., width = width, height = height, res = 300, units = "in")
      ggsave(file, plot = graph_object(), device = device)
    })
  
  output$download_enr_tbl=downloadHandler(
   filename = function() {
      paste("Enrichment_table",input$cell,"csv", sep = ".")
    },
    content = function(file) {
      write.csv(enr_object(), file, row.names = FALSE, col.names=TRUE, quote=FALSE, sep="\t")
    })
  
  output$download_barplot=downloadHandler(
    filename = function() {
      paste("Enrichment_Barplot",input$cell,"png", sep = ".")
    },
    content = function(file) {
      device <- function(..., width, height) grDevices::png(..., width = width, height = height, res = 300, units = "in")
      ggsave(file, plot = barplot_object(), device = device)
    })
  
  
}
shinyApp(ui = UI, server = server)
