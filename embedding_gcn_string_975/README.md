# A Biological Knowledge Graph for Representational Learning 
This repository provides the code for our research project "A Biological Knowledge Graph for Representational Learning".

## Data resources
The different dataset and KG used in this project are located in data directory. These files include:

-) The data about pathways from https://reactome.org/download/current/ReactomePathways.txt, relationships between pathways from https://reactome.org/download/current/ReactomePathwaysRelation.txt and pathway-protein relations from https://reactome.org/download/current/NCBI2Reactome.txt on 24 March 2024.

-) The built knowledge graph including pathway-pathway and pathway-protein relationships.

ppi

## Scripts
The code directory contains the following scripts:

-)The script for processing data download from Reactome

-)The script for building KG and save to Neo4j Aura.


## Setup
-)conda create -n kg python=3.10 -y

-)conda activate kg

-)pip install -r requirements.txt


## Get start

python embedding_pathway/gat_embedding.py --in_feats 256 --out_feats 256 --num_layers 2 --num_heads 2 --batch_size 1 --lr 0.0001 --num_epochs 5
python scripts/reactome/kg_reactome.py

python gat/edges_prediction.py --in-feats 128 --out-feats 32 --num-layers 2 --num-heads 2 --lr 0.0001 --epochs 20002


