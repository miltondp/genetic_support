#!/usr/bin/env Rscript

#' Load Additional Gene-Trait Associations
#' 
#' This script provides functionality to load additional sources of gene-trait 
#' associations and integrate them with the existing association data.
#'
#' @param additional_file Path to additional associations file (TSV format)
#' @param existing_file Path to existing associations file (default: 'data/assoc.tsv.gz')
#' @param output_file Path for merged output file (default: existing_file)
#' @param source_name Name for the new source (will be added to 'source' column)
#'
#' Expected schema for additional associations file:
#' - gene: Gene symbol (required)
#' - mesh_id: Medical Subject Heading ID (required) 
#' - original_trait: Original trait/phenotype name (required)
#' - pval: P-value (optional, numeric)
#' - year: Year of association (optional, integer)
#' - l2g_share: Locus-to-gene share score (optional, numeric 0-1)
#' - l2g_rank: Locus-to-gene rank (optional, integer)
#' - extra_info: Additional information (optional)
#' - original_link: Link to original study/database (optional)
#' - beta: Effect size beta (optional, numeric)
#' - odds_ratio: Odds ratio (optional, numeric)
#' - pic_qtl_pval: PiCCoLO QTL p-value (optional, numeric)
#' - pic_h4: PiCCoLO H4 score (optional, numeric)
#' - af_gnomad_nfe: Allele frequency in gnomAD NFE (optional, numeric)

suppressMessages(library(tidyverse))

#' Validate additional associations data
#' @param data Tibble with additional associations
#' @param source_name Name of the source for error messages
validate_associations <- function(data, source_name) {
  required_cols <- c("gene", "mesh_id", "original_trait")
  missing_cols <- setdiff(required_cols, colnames(data))
  
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns for source", source_name, ":", 
               paste(missing_cols, collapse = ", ")))
  }
  
  # Check for missing values in required columns
  for (col in required_cols) {
    if (any(is.na(data[[col]]) | data[[col]] == "")) {
      warning(paste("Source", source_name, "has missing values in required column:", col))
    }
  }
  
  # Validate data types for optional numeric columns
  numeric_cols <- c("pval", "l2g_share", "l2g_rank", "year", "beta", "odds_ratio", 
                   "pic_qtl_pval", "pic_h4", "af_gnomad_nfe")
  
  for (col in intersect(numeric_cols, colnames(data))) {
    if (col %in% colnames(data)) {
      # Try to convert to numeric, warn if issues
      suppressWarnings({
        numeric_vals <- as.numeric(data[[col]])
        if (sum(is.na(numeric_vals) & !is.na(data[[col]])) > 0) {
          warning(paste("Some values in column", col, "for source", source_name, 
                       "could not be converted to numeric"))
        }
      })
    }
  }
  
  # Validate l2g_share range if present
  if ("l2g_share" %in% colnames(data)) {
    invalid_shares <- !is.na(data$l2g_share) & (data$l2g_share < 0 | data$l2g_share > 1)
    if (any(invalid_shares)) {
      warning(paste("Source", source_name, "has", sum(invalid_shares), 
                   "l2g_share values outside valid range [0,1]"))
    }
  }
  
  TRUE
}

#' Load and merge additional associations
#' @param additional_file Path to additional associations file
#' @param existing_file Path to existing associations file  
#' @param output_file Path for output file
#' @param source_name Name for the new source
load_additional_associations <- function(additional_file, 
                                       existing_file = "data/assoc.tsv.gz",
                                       output_file = NULL,
                                       source_name = "NEW_SOURCE") {
  
  if (is.null(output_file)) {
    output_file <- existing_file
  }
  
  cat("Loading existing associations from:", existing_file, "\n")
  if (!file.exists(existing_file)) {
    stop(paste("Existing associations file not found:", existing_file))
  }
  
  # Load existing associations
  existing_assoc <- read_tsv(existing_file, col_types = cols(), show_col_types = FALSE)
  cat("Loaded", nrow(existing_assoc), "existing associations from", 
      length(unique(existing_assoc$source)), "sources\n")
  
  # Load additional associations
  cat("Loading additional associations from:", additional_file, "\n")
  if (!file.exists(additional_file)) {
    stop(paste("Additional associations file not found:", additional_file))
  }
  
  additional_assoc <- read_tsv(additional_file, col_types = cols(), show_col_types = FALSE)
  cat("Loaded", nrow(additional_assoc), "additional associations\n")
  
  # Validate additional associations
  validate_associations(additional_assoc, source_name)
  
  # Add source column
  additional_assoc$source <- source_name
  
  # Get column names from existing data to ensure compatibility
  existing_cols <- colnames(existing_assoc)
  
  # Add missing columns to additional data with appropriate defaults
  for (col in existing_cols) {
    if (!col %in% colnames(additional_assoc)) {
      if (col %in% c("year", "l2g_rank")) {
        additional_assoc[[col]] <- as.integer(NA)
      } else if (col %in% c("pval", "l2g_share", "beta", "odds_ratio", 
                            "pic_qtl_pval", "pic_h4", "af_gnomad_nfe")) {
        additional_assoc[[col]] <- as.numeric(NA)
      } else {
        additional_assoc[[col]] <- as.character(NA)
      }
    }
  }
  
  # Ensure additional data has same column order and types
  additional_assoc <- additional_assoc[existing_cols]
  
  # Add arow column (unique identifier)
  max_arow <- max(existing_assoc$arow, na.rm = TRUE)
  additional_assoc$arow <- seq(max_arow + 1, max_arow + nrow(additional_assoc))
  
  # Add mesh_term by joining with mesh names if available
  mesh_names_file <- "data/mesh_best_names.tsv.gz"
  if (file.exists(mesh_names_file)) {
    cat("Adding MeSH term names from:", mesh_names_file, "\n")
    mesh_names <- read_tsv(mesh_names_file, col_types = cols(), show_col_types = FALSE)
    additional_assoc <- additional_assoc %>%
      select(-mesh_term) %>%  # Remove if already present
      left_join(mesh_names %>% select(id, labeltext), by = c("mesh_id" = "id")) %>%
      rename(mesh_term = labeltext)
    
    # Put mesh_term in correct position
    additional_assoc <- additional_assoc[existing_cols]
  }
  
  # Combine datasets
  combined_assoc <- bind_rows(existing_assoc, additional_assoc)
  
  cat("Combined dataset has", nrow(combined_assoc), "associations from", 
      length(unique(combined_assoc$source)), "sources\n")
  
  # Summary by source
  source_summary <- combined_assoc %>%
    group_by(source) %>%
    summarize(count = n(), .groups = "drop") %>%
    arrange(desc(count))
  
  cat("\nAssociations by source:\n")
  print(source_summary)
  
  # Write output
  cat("\nWriting merged associations to:", output_file, "\n")
  write_tsv(combined_assoc, output_file, na = "")
  
  cat("Successfully merged additional associations!\n")
  return(combined_assoc)
}

# Command line interface
if (!interactive() && length(commandArgs(trailingOnly = TRUE)) > 0) {
  library(optparse)
  
  option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL,
                help = "Input file with additional associations (TSV format)", 
                metavar = "FILE"),
    make_option(c("-e", "--existing"), type = "character", 
                default = "data/assoc.tsv.gz",
                help = "Existing associations file [default %default]", 
                metavar = "FILE"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Output file (default: overwrite existing file)", 
                metavar = "FILE"),
    make_option(c("-s", "--source"), type = "character", default = "NEW_SOURCE",
                help = "Name for the new source [default %default]", 
                metavar = "NAME")
  )
  
  parser <- OptionParser(usage = "%prog [options]", option_list = option_list,
                        description = "Load additional gene-trait associations into genetic support pipeline")
  opt <- parse_args(parser)
  
  if (is.null(opt$input)) {
    cat("Error: Input file is required\n")
    print_help(parser)
    quit(status = 1)
  }
  
  # Run the function
  result <- load_additional_associations(
    additional_file = opt$input,
    existing_file = opt$existing, 
    output_file = opt$output,
    source_name = opt$source
  )
}