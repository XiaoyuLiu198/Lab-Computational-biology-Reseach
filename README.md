# Lab-Computational-biology-Reseach
Build ML pipeline of genomics data and reactive website to display

## UI: user interaction part
### Overview
This analysis platform include three main functions: Visualization, Genetic Anlaysis and Machine learning. 

##### Visualization:
Visualize edgelists with igraph and visNet. You can adjust threshold of edges showing and choose cells.
##### Genetic Anlysis:
1. Calculate centrality degree with chosen centrality measure.
2. Enrichment analysis with chosen enrichment method.(Based on filtered data)
##### Machine Learning:
1. Transform nodes into vectors with node2vec(https://snap.stanford.edu/node2vec/)
2. Classify genes with ML models.
![UI Overview](https://user-images.githubusercontent.com/65391883/119729599-75488680-be3a-11eb-80ac-faffcf5abc37.png)

#### Visualize the network 
In the visualization part, you can amplify the image. Different types of nodes are shown in different colors.
![Screen Shot 2021-05-26 at 3 44 36 PM](https://user-images.githubusercontent.com/65391883/119729966-da9c7780-be3a-11eb-9f17-f4141c2f46cf.png)
#### Genetic Anlaysis
Enrichment analysis based on centrality measures
![Screen Shot 2021-05-26 at 3 45 46 PM](https://user-images.githubusercontent.com/65391883/119730198-11728d80-be3b-11eb-8a7f-3da437aaebc0.png)

### Method: Machine learning methods involved
