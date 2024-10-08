# Bioinformatics Scripts Repository

This repository contains scripts for performing various bioinformatics analyses. 

## Current Scripts

### Protein Sequence Conservation Analysis

The program utilizes applications from EMBOSS (European Molecular Biology Open Software Suite), including Plotcon, Patmatmotifs, and Garnier, and the package PhyML. Given a user-defined subset of the taxonomic tree and the name of a protein family, the program will determine the level of conservation in the protein sequences in the dataset using Plotcon. The protein sequences included in the dataset are retrieved from the NCBI Protein database and exclude predicted and hypothetical sequences. Plotcon takes the sequence alignment created using Clustal Omega and plots the similarity found in each position in a given window of the protein sequence. Next, the sequences will be scanned using Patmatmotifs. This is to establish if motifs from the PROSITE database are associated with the dataset. Following this, the secondary structure prediction of the protein sequences will be performed with Garnier. The program applies the original Garnier Osguthorpe Robson algorithm for the predictions. Lastly, a phylogenetic tree will be constructed with PhyML. The phylogenetic tree will be restricted to a maximum of 100 sequences. The program is written for a Python3 interpreter.

## Dependencies

This repository requires the following dependencies:

- **EMBOSS**
- **PHYLIP**
- **Python 3**

Please ensure that all dependencies are installed before running the scripts.
