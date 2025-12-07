# =============================================================================
# Data Conversion Helper Functions
# Converts sylph/kraken results to phyloseq objects
# =============================================================================

#' Read Sylph Profile Data
#' @param file_path Path to sylph_profile.tsv or directory with .sylphmpa files
#' @return List with abundance matrix and taxonomy table
read_sylph_data <- function(file_path) {

  if (dir.exists(file_path)) {
    # Read individual .sylphmpa files
    files <- list.files(file_path, pattern = "\\.sylphmpa$", full.names = TRUE)
    if (length(files) == 0) {
      stop("No .sylphmpa files found in directory")
    }

    # Read all files (skip first line with #SampleID)
    all_data <- lapply(files, function(f) {
      read.table(f, header = TRUE, sep = "\t", comment.char = "#",
                 stringsAsFactors = FALSE, quote = "", skip = 1)
    })

    # Extract sample names from filenames
    sample_names <- gsub("\\.sylphmpa$", "", basename(files))
    sample_names <- gsub("_fastp.*", "", sample_names)  # Clean up

    # Combine data
    result <- combine_sylphmpa_files(all_data, sample_names)

  } else if (file.exists(file_path)) {
    # Read sylph_profile.tsv
    data <- read.table(file_path, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "")
    result <- parse_sylph_profile(data)

  } else {
    stop("File or directory not found: ", file_path)
  }

  return(result)
}

#' Parse sylph profile.tsv format
parse_sylph_profile <- function(data) {

  # Extract unique samples and genomes
  samples <- unique(data$Sample_file)
  samples <- gsub("_fastp.*\\.fastq\\.gz$", "", samples)

  # Create abundance matrix (Taxonomic_abundance)
  abundance_list <- split(data, data$Sample_file)

  # Get all unique genomes
  all_genomes <- unique(data$Genome_file)

  # Create matrix
  otu_mat <- matrix(0, nrow = length(all_genomes), ncol = length(samples))
  rownames(otu_mat) <- all_genomes
  colnames(otu_mat) <- samples

  # Fill matrix
  for (i in seq_along(samples)) {
    sample_data <- abundance_list[[unique(data$Sample_file)[i]]]
    match_idx <- match(sample_data$Genome_file, all_genomes)
    otu_mat[match_idx, i] <- sample_data$Taxonomic_abundance
  }

  # Create taxonomy table (extract from genome names)
  tax_mat <- create_taxonomy_from_genomes(all_genomes)

  return(list(
    abundance = otu_mat,
    taxonomy = tax_mat,
    metadata_extra = data[, c("Genome_file", "Adjusted_ANI", "Eff_cov")]
  ))
}

#' Combine multiple .sylphmpa files
combine_sylphmpa_files <- function(data_list, sample_names) {

  # Extract clade names from first file
  all_clades <- unique(data_list[[1]]$clade_name)

  # Get all unique clades across all samples
  for (i in seq_along(data_list)) {
    all_clades <- union(all_clades, data_list[[i]]$clade_name)
  }

  # Create abundance matrix
  otu_mat <- matrix(0, nrow = length(all_clades), ncol = length(sample_names))
  rownames(otu_mat) <- all_clades
  colnames(otu_mat) <- sample_names

  # Fill matrix with relative abundances
  for (i in seq_along(data_list)) {
    df <- data_list[[i]]
    match_idx <- match(df$clade_name, all_clades)
    otu_mat[match_idx, i] <- df$relative_abundance
  }

  # Create taxonomy table from clade names
  tax_mat <- parse_sylph_taxonomy(all_clades)

  return(list(
    abundance = otu_mat,
    taxonomy = tax_mat
  ))
}

#' Parse sylph MetaPhlAn-style taxonomy
parse_sylph_taxonomy <- function(clade_names) {

  tax_ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax_mat <- matrix(NA, nrow = length(clade_names), ncol = length(tax_ranks))
  rownames(tax_mat) <- clade_names
  colnames(tax_mat) <- tax_ranks

  for (i in seq_along(clade_names)) {
    parts <- strsplit(clade_names[i], "\\|")[[1]]

    for (j in seq_along(parts)) {
      # Extract rank prefix (d__, p__, c__, etc.)
      rank_prefix <- substr(parts[j], 1, 3)
      taxon_name <- substr(parts[j], 4, nchar(parts[j]))

      rank_idx <- switch(rank_prefix,
                        "d__" = 1,
                        "p__" = 2,
                        "c__" = 3,
                        "o__" = 4,
                        "f__" = 5,
                        "g__" = 6,
                        "s__" = 7,
                        NA)

      if (!is.na(rank_idx)) {
        tax_mat[i, rank_idx] <- taxon_name
      }
    }
  }

  # Fill in missing values with higher rank
  for (i in 1:nrow(tax_mat)) {
    for (j in 2:ncol(tax_mat)) {
      if (is.na(tax_mat[i, j])) {
        tax_mat[i, j] <- paste0("Unknown_", tax_mat[i, j-1])
      }
    }
  }

  return(tax_mat)
}

#' Create taxonomy from genome file names (GTDB format)
create_taxonomy_from_genomes <- function(genome_names) {

  tax_ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax_mat <- matrix(NA, nrow = length(genome_names), ncol = length(tax_ranks))
  rownames(tax_mat) <- genome_names
  colnames(tax_mat) <- tax_ranks

  # Extract species names from GCF/GCA identifiers
  for (i in seq_along(genome_names)) {
    # Simple extraction - can be improved with GTDB mapping
    acc <- gsub("_genomic\\.fna\\.gz$", "", genome_names[i])
    acc <- gsub("\\.\\d+", "", acc)  # Remove version numbers

    tax_mat[i, "Species"] <- acc
    tax_mat[i, "Genus"] <- paste0("Genus_", acc)
    tax_mat[i, "Family"] <- "Unknown"
  }

  return(tax_mat)
}

#' Read Kraken2/Bracken Report Data
#' @param dir_path Directory containing kraken2_bracken_species.report files
#' @return List with abundance matrix and taxonomy table
read_kraken_data <- function(dir_path) {

  files <- list.files(dir_path, pattern = "_kraken2_bracken_species\\.report$",
                     full.names = TRUE, recursive = TRUE)

  if (length(files) == 0) {
    stop("No kraken2_bracken_species.report files found")
  }

  # Read all files
  all_data <- lapply(files, read_kraken_report)

  # Extract sample names
  sample_names <- gsub("_kraken2_bracken_species\\.report$", "", basename(files))

  # Combine data
  result <- combine_kraken_reports(all_data, sample_names)

  return(result)
}

#' Read single kraken report
read_kraken_report <- function(file_path) {

  data <- read.table(file_path, header = FALSE, sep = "\t",
                    stringsAsFactors = FALSE, quote = "")

  colnames(data) <- c("percentage", "reads_clade", "reads_taxon",
                     "rank", "taxid", "name")

  # Clean up name (remove leading spaces)
  data$name <- trimws(data$name)

  # Extract taxonomy level from rank code
  data$rank_name <- sapply(data$rank, function(r) {
    switch(r,
           "D" = "Domain",
           "P" = "Phylum",
           "C" = "Class",
           "O" = "Order",
           "F" = "Family",
           "G" = "Genus",
           "S" = "Species",
           "U" = "Unclassified",
           "R" = "Root",
           "Other")
  })

  return(data)
}

#' Combine multiple kraken reports
combine_kraken_reports <- function(data_list, sample_names) {

  # Focus on species-level (rank = "S")
  species_data <- lapply(data_list, function(df) {
    df[df$rank == "S", ]
  })

  # Get all unique species
  all_species <- unique(unlist(lapply(species_data, function(df) df$name)))

  # Create abundance matrix
  otu_mat <- matrix(0, nrow = length(all_species), ncol = length(sample_names))
  rownames(otu_mat) <- all_species
  colnames(otu_mat) <- sample_names

  # Fill with percentages (relative abundance)
  for (i in seq_along(species_data)) {
    df <- species_data[[i]]
    match_idx <- match(df$name, all_species)
    otu_mat[match_idx, i] <- df$percentage
  }

  # Build taxonomy table from full kraken reports
  tax_mat <- build_kraken_taxonomy(data_list[[1]], all_species)

  return(list(
    abundance = otu_mat,
    taxonomy = tax_mat
  ))
}

#' Build taxonomy table from kraken report
build_kraken_taxonomy <- function(kraken_df, species_names) {

  tax_ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax_mat <- matrix(NA, nrow = length(species_names), ncol = length(tax_ranks))
  rownames(tax_mat) <- species_names
  colnames(tax_mat) <- tax_ranks

  # Simple extraction - would need full lineage tracking for complete accuracy
  for (i in seq_along(species_names)) {
    sp_name <- species_names[i]

    # Extract genus from species name
    genus <- strsplit(sp_name, " ")[[1]][1]

    tax_mat[i, "Species"] <- sp_name
    tax_mat[i, "Genus"] <- genus
    # Other ranks would need lineage tracking
  }

  return(tax_mat)
}

#' Create phyloseq object from parsed data
#' @param data_list List with abundance and taxonomy matrices
#' @param metadata Sample metadata data frame
#' @param group_column Column name for grouping variable
#' @return phyloseq object
create_phyloseq_from_data <- function(data_list, metadata, group_column) {

  # Create OTU table
  otu_mat <- otu_table(data_list$abundance, taxa_are_rows = TRUE)

  # Create taxonomy table
  tax_mat <- tax_table(as.matrix(data_list$taxonomy))

  # Handle metadata rownames
  # If metadata doesn't have rownames set, try to find sample ID column
  if (is.null(rownames(metadata)) || all(rownames(metadata) == as.character(1:nrow(metadata)))) {
    # Try common sample ID column names
    id_cols <- c("SampleID", "Sample_ID", "sample_id", "accession_used_in_analysis",
                 "sample", "Sample", "ID")
    found_col <- NULL
    for (col in id_cols) {
      if (col %in% colnames(metadata)) {
        found_col <- col
        break
      }
    }

    if (!is.null(found_col)) {
      rownames(metadata) <- metadata[[found_col]]
    } else {
      stop("Could not identify sample ID column in metadata. Please ensure metadata has rownames or a column named 'SampleID'")
    }
  }

  # Match metadata to samples in abundance data
  common_samples <- intersect(colnames(data_list$abundance), rownames(metadata))

  if (length(common_samples) == 0) {
    stop("No matching samples between abundance data and metadata")
  }

  # Subset and reorder
  otu_mat <- otu_mat[, common_samples]
  metadata <- metadata[common_samples, , drop = FALSE]

  # Ensure group column exists and is factor
  if (!group_column %in% colnames(metadata)) {
    stop(paste("Group column", group_column, "not found in metadata"))
  }

  metadata[[group_column]] <- as.factor(metadata[[group_column]])

  # Create sample data
  sample_data_obj <- sample_data(metadata)

  # Create phyloseq object
  ps <- phyloseq(otu_mat, tax_mat, sample_data_obj)

  # Convert to relative abundance (if not already)
  if (max(otu_table(ps)) > 100) {
    ps <- transform_sample_counts(ps, function(x) x / sum(x) * 100)
  }

  return(ps)
}
