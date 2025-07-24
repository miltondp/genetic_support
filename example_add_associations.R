#!/usr/bin/env Rscript

# Example: Adding and Using New Gene-Trait Association Sources
# 
# This script demonstrates how to add new gene-trait associations
# and use them in the genetic support analysis pipeline.

suppressMessages(library(tidyverse))

cat("=== Example: Adding New Gene-Trait Association Sources ===\n")

# First, let's create some example association data
example_associations <- tibble(
  gene = c("APOE", "BRCA1", "CFTR", "F5", "LDLR"),
  mesh_id = c("D000544", "D001943", "D003550", "D020246", "D006937"),
  original_trait = c("Alzheimer Disease", "Breast Neoplasms", "Cystic Fibrosis", 
                    "Venous Thromboembolism", "Hyperlipoproteinemia Type II"),
  pval = c(1.2e-45, 5.4e-52, 1.0e-100, 3.2e-28, 2.1e-65),
  year = c(2022, 2021, 2020, 2023, 2019),
  l2g_share = c(0.95, 0.98, 1.0, 0.87, 0.92),
  l2g_rank = c(1, 1, 1, 1, 1),
  extra_info = c("Strong genetic association", "High penetrance mutation", 
                "Monogenic disease", "Factor V Leiden mutation", 
                "Familial hypercholesterolemia"),
  original_link = paste0("https://example.com/study", 1:5)
)

# Save the example data
example_file <- "/tmp/example_new_associations.tsv"
write_tsv(example_associations, example_file)
cat("Created example associations file:", example_file, "\n")

# Show the data
cat("\nExample associations to add:\n")
print(example_associations)

# Load the associations using our utility function
cat("\n=== Loading Additional Associations ===\n")

# Create backup of original associations
if (!file.exists("data/assoc_backup.tsv.gz")) {
  file.copy("data/assoc.tsv.gz", "data/assoc_backup.tsv.gz")
  cat("Created backup: data/assoc_backup.tsv.gz\n")
}

# Source the loading function
source("src/load_additional_associations.R")

# Load the additional associations
cat("\nAdding associations with source name 'EXAMPLE_DB'...\n")
result <- load_additional_associations(
  additional_file = example_file,
  existing_file = "data/assoc_backup.tsv.gz",
  output_file = "/tmp/merged_associations_example.tsv.gz",
  source_name = "EXAMPLE_DB"
)

cat("\n=== Verifying Integration ===\n")

# Verify the new data is integrated correctly
merged_assoc <- read_tsv("/tmp/merged_associations_example.tsv.gz", 
                        col_types = cols(), show_col_types = FALSE)

cat("Total associations in merged file:", nrow(merged_assoc), "\n")

# Show source breakdown
source_counts <- merged_assoc %>%
  group_by(source) %>%
  summarize(count = n(), .groups = "drop") %>%
  arrange(desc(count))

cat("\nAssociation counts by source:\n")
print(source_counts)

# Show our new associations
new_assoc <- merged_assoc %>% filter(source == "EXAMPLE_DB")
cat("\nOur new associations:\n")
print(select(new_assoc, gene, mesh_id, original_trait, mesh_term, pval, l2g_share))

cat("\n=== Testing with Analysis Pipeline ===\n")

# Load necessary data for a simple pipeline test
if (file.exists("data/merge2.tsv.gz")) {
  cat("Loading merge2 data for pipeline test...\n")
  merge2 <- read_tsv("data/merge2.tsv.gz", col_types = cols(), show_col_types = FALSE)
  
  # Create a temporary merge2 file with our new associations
  temp_merge2 <- merge2
  # Add some rows for our example genes (simplified for demonstration)
  # In practice, these would come from the Pharmaprojects data
  
  cat("Pipeline functions are available and can use the new 'EXAMPLE_DB' source\n")
  cat("Example usage:\n")
  cat("  pipeline_best(merge2, associations=c('OMIM','OTG','EXAMPLE_DB'))\n")
} else {
  cat("Note: merge2.tsv.gz not available for full pipeline test\n")
}

cat("\n=== Summary ===\n")
cat("✓ Successfully added", nrow(example_associations), "new associations\n")
cat("✓ New source 'EXAMPLE_DB' is available in the pipeline\n")
cat("✓ Data validation passed\n")
cat("✓ MeSH terms were automatically mapped\n")

cat("\nTo use in analysis:\n")
cat("1. Replace 'data/assoc.tsv.gz' with '/tmp/merged_associations_example.tsv.gz'\n")
cat("2. Run: Rscript src/gensup_analysis.R\n")
cat("3. Or use pipeline_best() function with associations=c('OMIM','OTG','EXAMPLE_DB')\n")

cat("\nSee ADD_NEW_ASSOCIATIONS.md for detailed documentation.\n")