# Adding New Gene-Trait Association Sources

This guide explains how to add new sources of gene-trait associations to the genetic support pipeline.

## Overview

The genetic support pipeline currently includes associations from multiple sources:
- **OMIM**: Online Mendelian Inheritance in Man
- **OTG**: Open Targets Genetics
- **PICCOLO**: Phenome-wide association studies
- **Genebass**: UK Biobank genetic associations
- **intOGen**: Cancer driver gene associations

You can add additional sources using the `load_additional_associations.R` script.

## Quick Start

1. **Prepare your data** in TSV format with required columns:
   ```
   gene	mesh_id	original_trait	pval	year	l2g_share
   APOE	D000544	Alzheimer Disease	1.2e-45	2022	0.95
   BRCA1	D001943	Breast Neoplasms	5.4e-52	2021	0.98
   ```

2. **Load the new associations**:
   ```bash
   Rscript src/load_additional_associations.R \
     --input your_associations.tsv \
     --source "YOUR_DATABASE_NAME"
   ```

3. **Run the analysis** with your new source:
   ```bash
   Rscript src/gensup_analysis.R
   ```

## Data Format Requirements

### Required Columns
- **gene**: Gene symbol (e.g., "APOE", "BRCA1")
- **mesh_id**: Medical Subject Heading ID (e.g., "D000544" for Alzheimer Disease)
- **original_trait**: Original trait/phenotype name

### Optional Columns
- **pval**: P-value (numeric)
- **year**: Year of association (integer)
- **l2g_share**: Locus-to-gene share score (numeric, 0-1)
- **l2g_rank**: Locus-to-gene rank (integer)
- **extra_info**: Additional information (text)
- **original_link**: Link to original study/database (URL)
- **beta**: Effect size beta (numeric)
- **odds_ratio**: Odds ratio (numeric)
- **pic_qtl_pval**: PiCCoLO QTL p-value (numeric)
- **pic_h4**: PiCCoLO H4 score (numeric)
- **af_gnomad_nfe**: Allele frequency in gnomAD NFE (numeric)

### Notes
- Missing optional columns will be filled with appropriate defaults (NA)
- MeSH terms will be automatically added if the mesh_id exists in the reference database
- Associations will be assigned unique row identifiers automatically

## Command Line Usage

```bash
Rscript src/load_additional_associations.R [options]

Options:
  -i FILE, --input=FILE       Input file with additional associations (TSV format)
  -e FILE, --existing=FILE    Existing associations file [default: data/assoc.tsv.gz]
  -o FILE, --output=FILE      Output file (default: overwrite existing file)
  -s NAME, --source=NAME      Name for the new source [default: NEW_SOURCE]
  -h, --help                  Show help message and exit
```

### Examples

**Basic usage** (overwrites existing associations file):
```bash
Rscript src/load_additional_associations.R \
  -i new_associations.tsv \
  -s "MyDatabase"
```

**Save to new file** (preserves original):
```bash
Rscript src/load_additional_associations.R \
  -i new_associations.tsv \
  -o data/assoc_with_mydb.tsv.gz \
  -s "MyDatabase"
```

**Use different existing file**:
```bash
Rscript src/load_additional_associations.R \
  -i new_associations.tsv \
  -e data/custom_assoc.tsv.gz \
  -s "MyDatabase"
```

## Using New Sources in Analysis

Once you've added new associations, you can include them in the genetic support analysis by specifying the source name in the `associations` parameter of the `pipeline_best()` function.

For example, to analyze only your new source:
```r
# In R
library(tidyverse)
source('src/gensup_analysis.R')
merge2 <- read_tsv('data/merge2.tsv.gz', col_types=cols())

# Run analysis with your new source
result <- pipeline_best(merge2, 
                       phase='combined', 
                       basis='ti', 
                       associations=c('MyDatabase'))
```

Or combine with existing sources:
```r
result <- pipeline_best(merge2, 
                       phase='combined', 
                       basis='ti', 
                       associations=c('OMIM', 'OTG', 'MyDatabase'))
```

## Validation

The script performs several validation checks:
- **Required columns**: Ensures gene, mesh_id, and original_trait are present
- **Data types**: Validates numeric columns have appropriate values
- **Range checks**: Ensures l2g_share values are between 0 and 1
- **Missing values**: Warns about missing values in required columns

## MeSH ID Mapping

The script automatically maps MeSH IDs to standardized terms using `data/mesh_best_names.tsv.gz`. If your MeSH IDs are not found in this reference, the `mesh_term` will be left as NA, but the association will still be included.

For best results, ensure your `mesh_id` values follow the standard MeSH format (e.g., "D000544").

## Troubleshooting

**Error: Missing required columns**
- Ensure your input file has `gene`, `mesh_id`, and `original_trait` columns
- Check that column names are spelled correctly (case-sensitive)

**Warning: Values could not be converted to numeric**
- Check that numeric columns (pval, year, l2g_share, etc.) contain valid numbers
- Missing values should be empty or "NA"

**File not found errors**
- Ensure file paths are correct relative to the repository root
- Use absolute paths if working from a different directory

## Integration with Existing Pipeline

The new associations are seamlessly integrated with the existing pipeline. The main analysis script `src/gensup_analysis.R` will automatically detect and use the new associations when they are present in the data files.

The pipeline preserves all existing functionality while adding support for your new sources.