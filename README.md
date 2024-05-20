# IDRBlocks

This repo contains scripts and utilities to predict the partition coefficients of proteins into Mediator condensates as measured in this paper: https://www.sciencedirect.com/science/article/pii/S0092867422015264?via%3Dihub

The raw data of proteins and their partition coefficients are present in the directory ``raw_data``. The directory ``scripts`` contains code to run iprscan bioinformatics analysis, extract sequences of intrinsically disordered regions (IDRs) of proteins, and extract domain-knowledge informed features from these IDR sequences. 

The ``Top_200_analysis`` directory contains analysis of the top-200 proteins with the highest and lowest partition coefficients, just to check whether the way of featurizing sequences contains any information about partition ratios at all, and it seems that there is some information there. 

Therefore, I went ahead and started exploring the full dataset more, the analysis of which is present in ``notebooks``

Finally, the ``.ipynb`` files in the main directory are the Colab notebooks that I have used to train the models and explore how using domain knowledge-based features and embeddings from LLMs can predict partitioning. 


