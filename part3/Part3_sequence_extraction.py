#Generative AI was not used in the preparation of this script.

import pandas as pd
from Bio import SeqIO
import os

#target gene symbols from EcoCyc
target_degs = ["nikA", "arsC", "stpA", "lrp", "argR"] 

#extract relevant files from part1
pipelines = {
    "Flye_Prokka": {"tsv": "/home/hannah/LIFE748asm/part1/annotation/annotations/prokka_flye.tsv", 
                    "faa": "/home/hannah/LIFE748asm/part1/annotation/annotations/prokka_flye.faa"},
    "SPAdes_Prokka": {"tsv": "/home/hannah/LIFE748asm/part1/annotation/annotations/prokka_spades.tsv", 
                      "faa": "/home/hannah/LIFE748asm/part1/annotation/annotations/prokka_spades.faa"},
    "Flye_Bakta": {"tsv": "/home/hannah/LIFE748asm/part1/annotation/annotations/bakta_flye.tsv", 
                   "faa": "/home/hannah/LIFE748asm/part1/annotation/annotations/bakta_flye.faa"},
    "SPAdes_Bakta": {"tsv": "/home/hannah/LIFE748asm/part1/annotation/annotations/bakta_spades.tsv",
                    "faa": "/home/hannah/LIFE748asm/part1/annotation/annotations/bakta_spades.faa"}
}

#create empty list to store the records
extracted_records = []
#loop through each file to ensure correct extraction for each file format
for name, files in pipelines.items():
    if "Bakta" in name:
        #skip the first 5 lines of bakta to avoid errors
        df = pd.read_csv(files['tsv'], sep='\t', header=5)
    else:
        #standard read-in for prokka files
        df = pd.read_csv(files['tsv'], sep='\t')
    
    #clean the headers for bakta files
    df.columns = [c.lstrip('#').strip().lower().replace(' ', '_') for c in df.columns]
    
    #map locus tag to gene name for target DEGs
    subset = df[df['gene'].astype(str).str.lower().isin([g.lower() for g in target_degs])]
    tag_to_gene = dict(zip(subset['locus_tag'], subset['gene']))
    
    #extract sequences from fasta file
    for record in SeqIO.parse(files['faa'], "fasta"):
        if record.id in tag_to_gene:
            gene_symbol = tag_to_gene[record.id]
            #create unique id using the gene name
            record.id = f"{gene_symbol}"
            record.description = f"Source_Locus:{name}"
            extracted_records.append(record)

#save the results
extracted_records.sort(key=lambda record: record.id)
SeqIO.write(extracted_records, "Comparative_DEG_Proteins.fasta", "fasta")
print(f"Successfully extracted {len(extracted_records)} sequences.")