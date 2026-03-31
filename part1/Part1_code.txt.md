Genome assembly

* flye --pacbio-hifi GN6\_hifix30.fastq --out-dir GN6\_flye
* spades.py --isolate --only-assembler   -t 8 -m 12   -s GN6\_hifix30.fastq   -o GN6\_spades
* quast.py assemblies/\*.fasta \\

&#x20; 	-r references/EcoliK12\_GCA\_000005845.2\_ASM584v2\_genomic.fna \\

&#x20; 	-g references/EcoliK12\_GCA000005845.gff \\

&#x20; 	--threads 4 \\

&#x20; 	-o quast/quast\_comparison\_report



Genome annotation

* bakta --db \~/tmp/db-light assemblies/flye\_assembly.fasta --output bakta/flye
* bakta --db \~/tmp/db-light assemblies/spades\_contigs.fasta --output bakta/spades
* prokka --outdir prokka/flye --prefix flye assemblies/flye\_assembly.fasta
* prokka --outdir prokka/spades --prefix spades assemblies/spades\_contigs.fasta



Preparations for data analysis

* All bakta and prokka output files were copied into the "annotations" folder and renamed accordingly

