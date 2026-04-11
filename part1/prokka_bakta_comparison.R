library(readr)
library(dplyr)
library(Biostrings)
library(tidyr)
library(ggplot2)

setwd('//wsl.localhost/Ubuntu/home/hannah/LIFE748asm/annotations')

#
prokka_tsv <- list.files(pattern = "^prokka_.*\\.tsv$")
bakta_tsv  <- list.files(pattern = "^bakta_.*\\.tsv$")

prokka_samples <- sub("prokka_", "", tools::file_path_sans_ext(prokka_tsv))
bakta_samples  <- sub("bakta_", "", tools::file_path_sans_ext(bakta_tsv))

samples <- intersect(prokka_samples, bakta_samples)
all_results <- list()
print(all_results)

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
    prokka_func   = p_func,  
    prokka_hypo   = p_hypo,  
    bakta_func    = b_func,  
    bakta_hypo    = b_hypo,  
    prokka_func_pct = p_func_pct,
    bakta_func_pct  = b_func_pct
  ) 
}
print(all_results)
summary_df <- bind_rows(all_results)
write_csv(summary_df, "annotation_summary.csv")

#preparation for plotting
ann <- read_csv("annotation_summary.csv")

plot_func <- summary_df %>%
  select(sample, prokka_func_pct, bakta_func_pct) %>%
  pivot_longer(
    cols = ends_with("_pct"),
    names_to = "tool",
    values_to = "frac_function"
  ) %>%
  mutate(tool = stringr::str_to_title(sub("_func_pct", "", tool)),
         sample = stringr::str_to_title(sample))

p1 <- ggplot(gene_long, aes(x = sample, y = count, fill = tool)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(
                aes(label = count),
                position = position_dodge(width = 0.9),
                vjust = -0.5,    # Pushes the text slightly above the bar**
              size = 3.5,      # Adjust text size as needed**
              fontface = "bold"
                ) +
  
  facet_wrap(~type, scales = "free_y", ncol = 1) + 
  scale_fill_manual(values = c("Prokka" = "#0072B2", "Bakta" = "#D55E00")) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "Assembly Sample",
    y = "Feature Count",
    fill = "Annotation Tool",
    title = "Genomic Feature Comparison: Bakta vs. Prokka"
  )
p1
# -------------------------
# 2. Assembly size vs CDS
# -------------------------
plot_data_cds <- ann %>%
  select(sample, genome_size, prokka_CDS, bakta_CDS) %>%
  pivot_longer(
    cols = c(prokka_CDS, bakta_CDS),
    names_to = "tool",
    values_to = "cds"
  )%>%
  mutate(tool = recode(tool,
                       "prokka_CDS" = "Prokka",
                       "bakta_CDS"  = "Bakta"))

p2 <- ggplot(plot_data_cds, aes(x = genome_size / 1e6, y = cds, colour = tool)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_colour_manual(values = c("Prokka" = "#0072B2", "Bakta" = "#D55E00")) +
  theme_classic() +
  labs(
    title = "Assembly Size vs CDS Count",
    x = "Assembly Size (Mb)",
    y = "Number of CDS",
    colour = "Tool"
  )

# -------------------------
# 3. Functional vs hypothetical
# -------------------------
plot_data_comp <- summary_df %>%
  select(sample, prokka_func, prokka_hypo, bakta_func, bakta_hypo) %>%
  pivot_longer(
    cols = -sample,
    names_to = c("tool", "category"),
    names_sep = "_"
  ) %>%
  mutate(
    category = recode(category, "func" = "Functional", "hypo" = "Hypothetical"),
    tool = stringr::str_to_title(tool),
    sample = stringr::str_to_title(sample)
  )

p3 <- ggplot(plot_data_comp, aes(x = tool, y = value, fill = category)) +
  geom_col(position = "fill", color = "white") +
  facet_wrap(~sample) +
  scale_fill_manual(values = c("Functional" = "#0072B2", "Hypothetical" = "#D55E00")) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  labs(
    title = "Functional vs Hypothetical CDS",
    x = "Tool",
    y = "Proportion",
    fill = "Category"
  )

p1
p2
p3
