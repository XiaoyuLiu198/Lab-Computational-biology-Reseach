#!/usr/bin/env python
# coding: utf-8

# In[ ]:


class pipeline:
    def __init__(self,file1,file2):
        ##load disease-gene matching file and edgelist
        self.edgelist=pd.read_table(file1, sep=' ')
        self.graph=nx.read_edgelist(file2)
        self.disease_gene=pd.read_csv('all_gene_disease_associations.tsv', sep='\t')
        
    def node2vec(self,training):
        if training==1:
            import networkx as nx
            from node2vec import Node2Vec

            # Precompute probabilities and generate walks - **ON WINDOWS ONLY WORKS WITH workers=1**
            node2vec = Node2Vec(graph, dimensions=64, walk_length=30, num_walks=200, workers=4)  # Use temp_folder for big graphs

            # Embed nodes
            model = node2vec.fit(window=10, min_count=1, batch_words=4)  # Any keywords acceptable by gensim.Word2Vec can be passed, `dimensions` and `workers` are automatically passed (from the Node2Vec constructor)
            # Save embeddings for later use
            model.wv.save_word2vec_format("save_embedding")
            vector_table=pd.read_table("save_embedding",sep=" ",skiprows=1,header=None)
        else:
            vector_table=pd.read_table("save_embedding",sep=" ",skiprows=1,header=None)
            return vector_table
        
    def preprocess(self,training):
        ##collect the gene of AD in complete set
        vector_table=self.node2vec(training)
        AD_list=self.disease_gene[self.disease_gene['diseaseName']=="Alzheimer's Disease"]['geneSymbol'].unique()
        ##match
        self.geneset['genes'] = self.geneset.apply(lambda x: list([x['regulatoryGene'],
                                                x['targetGene']]),axis=1)      

        self.geneset['label']=self.geneset['genes'].map(lambda x: 1 if x[0] in AD_list or x[1] in AD_list else 0)
        ##transform vectors to standard format
        col_list=["geneName"]
        for i in range(1,65):
            col_list.append(str("vec"+str(i)))
        vector_table.columns=col_list
        ##merge edgelist with vectors
        dataset=self.geneset.merge(vector_table,left_on='regulatoryGene',right_on='geneName')
        dataset=self.dataset.merge(vector_table,left_on='targetGene',right_on='geneName')
        dataset=self.dataset.drop(columns=['geneName_y'],axis=1)
        return dataset
    
    def cross_validation(self,training,parameters1,parameters2):
        from sklearn.model_selection import train_test_split
        dataset=self.preprocess(training)
        X_train, X_test, y_train, y_test = train_test_split(trainset.iloc[:,6:134], 
                                                            trainset['label'], test_size=0.33, random_state=42)
        from sklearn.svm import SVC
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.model_selection import GridSearchCV
        rf=RandomForestClassifier()
        rf_grid=GridSearchCV(rf,parameters1,cv=3,n_jobs=3)
        rf_grid.fit(X_train,y_train)
        svm=SVC()
        svm_grid=GridSearchCV(svm,parameters2,cv=3,n_jobs=3)
        svm_grid.fit(X_train,y_train)
        return rf_grid.best_params_,svm_grid.best_params_
    


# In[ ]:


parameters1={'n_estimators':[200,225,250,275],'max_depth':[30,40,50],'max_features':['sqrt']}
parameters2={'kernel':['poly','rbf','sigmoid'],'degree':[2,3,4,5]}
ml=pipeline(ex_genie3.txt,edgelist.txt)
p1,p2=ml.cross_validation(0,parameters1,parameters2)
print(p1,p2)


# In[ ]:



from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_auc_score
def evaluate(model, test_features, test_labels):
    predictions = model.predict(test_features)
    avgPrecision=average_precision_score(test_labels, predictions)
    ROCAUC=roc_auc_score(test_labels, predictions)
    print('Model Performance')
    print('Average Precision: {:0.4f} degrees.'.format(avgPrecision))
    print('ROC AUC: {:0.4f} degrees.'.format(ROCAUC))
    return avgPrecision,ROCAUC

rf_model = RandomForestClassifier(n_estimators = 200,max_features='sqrt',
                                   max_depth=40)
rf_model.fit(X_train,y_train)
rf_accuracy,rf_roc_auc_score = evaluate(rf_model,X_test,y_test)

svm_model = SVC(kernel = 'poly',degree=3)
svm_model.fit(X_train,y_train)
svm_accuracy,svm_roc_auc_score = evaluate(svm_model,X_test,y_test)

