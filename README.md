# Lab-Computational-biology-Reseach
Build ML pipeline of genomics data and reactive website to display

## Machine learning pipeline
1. Pre-process dataset by filtering and transposing.
2. Transform nodes into vectors with node2vec(https://snap.stanford.edu/node2vec/)
3. Classify genes with ML models including XGBoost, random forest, linear reression and deep learning models. Tuned hyper parameter and validated classification with cross validation.

## Website overview
This analysis platform include two main functions: Visualization and Genetic Anlaysis.

#### Network visualization based on Igraph and VisNet(R)
Visualize genomics edgelists with igraph and visNet. The framework enables adjusting threshold of edges shown, and choose cells to display.
#### Genetic Anlysis: interactive simulation with R
1. Calculate centrality degree with chosen centrality measure methods.
2. Enrichment analysis with chosen enrichment methods.(Based on filtered data)

## Demo
![UI Overview](https://user-images.githubusercontent.com/65391883/119729599-75488680-be3a-11eb-80ac-faffcf5abc37.png)

#### Visualize the network 
In the visualization part, you can amplify the image. Different types of nodes are shown in different colors to display the relation between genes in network.
![Screen Shot 2021-05-26 at 3 44 36 PM](https://user-images.githubusercontent.com/65391883/119729966-da9c7780-be3a-11eb-9f17-f4141c2f46cf.png)
#### Genetic Anlaysis
Enrichment analysis based on centrality measures.
![Screen Shot 2021-05-26 at 3 45 46 PM](https://user-images.githubusercontent.com/65391883/119730198-11728d80-be3b-11eb-8a7f-3da437aaebc0.png)
