suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(tidyverse)
  library(tibble)
  library(ggrepel)
})

option_list <- list(
  make_option("--counts", type = "character"),
  make_option("--metadata", type = "character"),
  make_option("--case", type = "character"),
  make_option("--control", type = "character"),
  make_option("--output", type = "character"),
  make_option("--volcano", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

counts <- read.csv(opt$counts, row.names = 1, check.names = FALSE)
counts_t <- t(counts)
metadata <- read.csv(opt$metadata, stringsAsFactors = FALSE)
metadata <- metadata %>% tibble::column_to_rownames("sample_id")
metadata <- metadata[colnames(counts_t), , drop = FALSE]

if (any(metadata$condition == opt$case) == FALSE || any(metadata$condition == opt$control) == FALSE) {
  stop("Case or control labels not found in metadata")
}

dds <- DESeqDataSetFromMatrix(countData = counts_t,
                              colData = metadata,
                              design = ~ condition)

dds$condition <- relevel(dds$condition, ref = opt$control)

dds <- DESeq(dds)
coef_name <- paste0("condition_", opt$case, "_vs_", opt$control)
if (!(coef_name %in% resultsNames(dds))) {
  message("Coefficient ", coef_name, " not found. Falling back to Wald test with explicit contrast.")
  res <- results(dds, contrast = c("condition", opt$case, opt$control))
} else {
  res <- lfcShrink(dds, coef = coef_name, type = "apeglm")
}
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df %>% arrange(padj)

readr::write_csv(res_df, opt$output)

volcano <- res_df %>%
  mutate(sig = if_else(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 0.5, "Significant", "NS"))

g <- ggplot(volcano, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.8, size = 1.6) +
  scale_color_manual(values = c("Significant" = "#d73027", "NS" = "#bdbdbd")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") +
  labs(title = paste("Pseudobulk DE:", opt$case, "vs", opt$control),
       x = "log2 Fold Change",
       y = "-log10 adj. p-value",
       color = "") +
  theme_minimal()

ggsave(filename = opt$volcano, plot = g, width = 6, height = 5, dpi = 300)
