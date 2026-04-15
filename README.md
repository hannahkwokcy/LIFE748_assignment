# LIFE748_assignment
This repository contains the intermediate and final outputs of the LIFE748 assignment 2 for genome assembly and annotation, machine learning and structural bioinformatics on E. Coli sequencing data. 

# Description
### `part1` Directory
* `annotation`: Contains all the files needed for genome annotation.
  * `bakta` and `prokka` : Contain the genome annotation pipeline outputs for each Flye and SPAdes assembly.
  * `annotations` : Contains copies of the resulting `.faa`, `.fna`, and `.tsv` files from the four annotation outputs for easier access.
* * `assembly`: Contains all the files needed for genome assembly.
  * `flye` and `spades`: Contain the outputs of the respective genome assembly pipelines.
  * `assemblies`: Contains copies of the final FASTA files for easy access.
  * `references`: Contains the E. Coli reference genome required for Quast analysis.
  * `quast`: Contains the results of the Quast analysis.
* part1_graphs.r: Contains the R code used for generating summary graphs.
* part1_code.txt: Contains all the commands and codes used to run the pipeline tools.

### `part2` Directory
* Assessment2_ML_2.qmd: Contains the code used to perform both unsupervised and supervised machine learning.
* significant_genes_with_FC.txt: Contains the extracted genes that have significant fold changes.

### `part3` Directory
* New-SmartTable---2026-03-17T13_46_15--8.txt: Contains the extracted genes paired with their corresponding gene names and GO terms (generated using the EcoCyc E. coli database) for downstream analysis.
* sequence_extraction.py: Script used to extract sequences from all Bakta and Prokka assemblies generated in Part 1.
* Comparative_DEG_Proteins.fasta: Contains the results of the sequence extraction.
* `{gene_name}` : Contain all AlphaFold2 prediction results based on the extracted sequences, along with their relevant PyMOL files. 
* Part3_pymol_code.txt: Contains the specific commands used in PyMOL to generate the desired 3D structures and images.
