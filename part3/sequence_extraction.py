import pandas as pd
from Bio import SeqIO
import os

# 1. Your target gene symbols from EcoCyc
target_degs = ["nikA", "arsC", "stpA", "lrp", "argR"] 

# 2. Define pairs of (TSV, FAA) for each pipeline
pipelines = {
    "Flye_Prokka": {"tsv": "/home/hannah/LIFE748asm/part1/annotations/prokka_flye.tsv", 
                    "faa": "/home/hannah/LIFE748asm/part1/annotations/prokka_flye.faa"},
    "SPAdes_Prokka": {"tsv": "/home/hannah/LIFE748asm/part1/annotations/prokka_spades.tsv", 
                      "faa": "/home/hannah/LIFE748asm/part1/annotations/prokka_spades.faa"},
    "Flye_Bakta": {"tsv": "/home/hannah/LIFE748asm/part1/annotations/bakta_flye.tsv", 
                   "faa": "/home/hannah/LIFE748asm/part1/annotations/bakta_flye.faa"},
    "SPAdes_Bakta": {"tsv": "/home/hannah/LIFE748asm/part1/annotations/bakta_spades.tsv",
                    "faa": "/home/hannah/LIFE748asm/part1/annotations/bakta_spades.faa"}
}

extracted_records = []
for name, files in pipelines.items():
    # 1. Load the TSV based on the tool's specific format
    if "Bakta" in name:
        # Skips the first 5 lines of metadata
        df = pd.read_csv(files['tsv'], sep='\t', header=5)
    else:
        # Standard Prokka format (header on line 1)
        df = pd.read_csv(files['tsv'], sep='\t')
    
    # 2. Clean the headers (Crucial: Bakta's header starts with #)
    df.columns = [c.lstrip('#').strip().lower().replace(' ', '_') for c in df.columns]
    
    # 3. Map Locus Tag to Gene Name for your target DEGs
    # Using .astype(str) handles cases where 'gene' might be NaN
    subset = df[df['gene'].astype(str).str.lower().isin([g.lower() for g in target_degs])]
    tag_to_gene = dict(zip(subset['locus_tag'], subset['gene']))
    
    # 4. Extract sequences from the .faa file
    for record in SeqIO.parse(files['faa'], "fasta"):
        if record.id in tag_to_gene:
            gene_symbol = tag_to_gene[record.id]
            # Unique ID for the final FASTA: e.g., nikA_Flye_Bakta
            record.id = f"{gene_symbol}"
            record.description = f"Source_Locus:{name}"
            extracted_records.append(record)

# 5. Save the results
extracted_records.sort(key=lambda record: record.id)
SeqIO.write(extracted_records, "Comparative_DEG_Proteins.fasta", "fasta")
print(f"Successfully extracted {len(extracted_records)} sequences.")