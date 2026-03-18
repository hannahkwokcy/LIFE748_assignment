library(readr)
library(dplyr)
library(Biostrings)
library(tidyr)
library(ggplot2)

setwd('//wsl.localhost/Ubuntu/home/hannah/LIFE748asm/annotations')

# 1. Identify common sample names
prokka_tsv <- list.files(pattern = "^prokka_.*\\.tsv$")
bakta_tsv  <- list.files(pattern = "^bakta_.*\\.tsv$")

prokka_samples <- sub("prokka_", "", tools::file_path_sans_ext(prokka_tsv))
bakta_samples  <- sub("bakta_", "", tools::file_path_sans_ext(bakta_tsv))

samples <- intersect(prokka_samples, bakta_samples)
all_results <- list()

for (s in samples) {
  message("Processing sample: ", s)
  
  p_path <- paste0("prokka_", s, ".tsv")
  b_path <- paste0("bakta_", s, ".tsv")
  f_path <- paste0("bakta_", s, ".fna") 
  
  # --- PROCESS PROKKA ---
  # Prokka is usually cleaner, but we'll use the same robust method
  p_data <- read.table(p_path, sep = "\t", header = TRUE, comment.char = "#", 
                       fill = TRUE, quote = "", stringsAsFactors = FALSE)
  colnames(p_data) <- tolower(colnames(p_data))
  
  p_n_cds  <- sum(p_data$ftype == "CDS", na.rm = TRUE)
  p_n_trna <- sum(p_data$ftype == "tRNA", na.rm = TRUE)
  p_n_rrna <- sum(p_data$ftype == "rRNA", na.rm = TRUE)
  p_func   <- sum(!is.na(p_data$product) & !grepl("hypothetical", p_data$product, ignore.case = TRUE))
  p_func_pct <- p_func / p_n_cds
  p_hypo = p_n_cds - p_func
  
  # --- PROCESS BAKTA ---
  b_data <- read.delim(b_path, sep = "\t", header = TRUE, comment.char = "#", 
                       quote = "", fill = TRUE, check.names = FALSE)
  colnames(b_data) <- c(
    "sequence_id", "type", "start", "stop",
    "strand", "locus_tag", "gene", "product"
  )
  
  b_n_cds  <- sum(b_data$type == "cds", na.rm = TRUE)
  b_n_trna <- sum(b_data$type == "tRNA", na.rm = TRUE)
  b_n_rrna <- sum(b_data$type == "rRNA", na.rm = TRUE)
  
  b_cds <- b_data[toupper(b_data$type) == "CDS", ]
  
  b_func <- sum(!is.na(b_cds$product) & 
                  !grepl("hypothetical", b_cds$product, ignore.case = TRUE)) 
  b_func_pct <- b_func / b_n_cds
  b_hypo  = b_n_cds - b_func
  # --- PROCESS FASTA ---
  fa <- readDNAStringSet(f_path)
  total_bp <- sum(width(fa))
  gc_fraction <- sum(letterFrequency(fa, c("G", "C"))) / total_bp
  
  # Store results
  all_results[[s]] <- tibble(
    sample        = s,
    genome_size   = total_bp,
    gc_content    = gc_fraction,
    prokka_CDS    = p_n_cds,
    bakta_CDS     = b_n_cds,
    prokka_tRNA   = p_n_trna,
    bakta_tRNA    = b_n_trna,
    prokka_rRNA   = p_n_rrna,
    bakta_rRNA    = b_n_rrna,
    prokka_func   = p_func,   # STORED AS COLUMN
    prokka_hypo   = p_hypo,   # STORED AS COLUMN
    bakta_func    = b_func,   # STORED AS COLUMN
    bakta_hypo    = b_hypo,   # STORED AS COLUMN
    prokka_func_pct = p_func_pct,
    bakta_func_pct  = b_func_pct
  ) 
}
print(all_results)
summary_df <- bind_rows(all_results)
write_csv(summary_df, "annotation_summary.csv")

#preparation for plotting
ann <- read_csv("annotation_summary.csv")

gene_long <- ann |>
  select(sample, prokka_CDS, bakta_CDS, prokka_tRNA, bakta_tRNA) |>
  pivot_longer(-sample, names_to = "feature", values_to = "count")

#functional annotation coverage plot
ggplot(gene_long, aes(x = sample, y = count, fill = feature)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    x = "Genome",
    y = "Feature count",
    fill = "Feature type"
  )

#assembly size vs gene count plot
plot_data <- ann |>
  select(sample, genome_size, prokka_CDS, bakta_CDS) |>
  pivot_longer(
    cols = c(prokka_CDS, bakta_CDS),
    names_to = "tool",
    values_to = "cds"
  )

ggplot(plot_data, aes(x = genome_size / 1e6, y = cds, colour = tool)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Assembly size (Mb)",
    y = "Number of CDS",
    colour = "Annotation tool"
  )

#contiguity vs functional coverage
gene_long <- gene_long |>
  separate(feature, into = c("tool", "gene_type"), sep = "_")

ggplot(gene_long,
       aes(x = sample, y = count, fill = tool)) +
  geom_col(position = "dodge") +
  facet_wrap(~ gene_type, scales = "free_y") +
  labs(
    x = "Sample",
    y = "Gene count",
    fill = "Annotation tool"
  ) +
  theme_minimal()

#functional vs hypothetical counts
# 1. Transform to Long Format
plot_data <- summary_df %>%
  select(sample, p_func, p_hypo, b_func, b_hypo) %>%
  pivot_longer(cols = -sample, 
               names_to = c("tool", "category"), 
               names_sep = "_") %>%
  mutate(category = recode(category, "func" = "Functional", "hypo" = "Hypothetical"),
         tool = str_to_title(tool))

# 2. Create Stacked Bar Plot
ggplot(plot_data, aes(x = tool, y = value, fill = category)) +
  geom_col(position = "stack", color = "white") +
  facet_wrap(~sample) +
  scale_fill_manual(values = c("Functional" = "#2c7fb8", "Hypothetical" = "#a1dab4")) +
  labs(
    title = "Annotation Benchmarking: Functional vs. Hypothetical CDS",
    subtitle = "Comparing predicted gene quality between Prokka and Bakta",
    x = "Annotation Tool",
    y = "Total CDS Count",
    fill = "Annotation Type"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold", size = 11))
