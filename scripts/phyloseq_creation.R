
library("phyloseq"); packageVersion("phyloseq")
library(dplyr)
#  ‘1.46.0’ 



phyloseq_creation <- function(OTU.tsv, TAX.csv , SAMPLE.csv, sample_column_index) {
  # Read the tables
  OTU <- read.csv(OTU.tsv, sep = '\t', header = TRUE, row.names = 1)
  TAX <- read.csv(TAX.csv, header = TRUE, row.names = 1)
  SAMPLE <- read.csv(SAMPLE.csv, header = TRUE, row.names = sample_column_index,
                     stringsAsFactors = FALSE)
  # Create the phyloseq object
  OTU <- otu_table(as.matrix(OTU), taxa_are_rows = TRUE)
  TAX <- tax_table(as.matrix(TAX))
  SAMPLE <- sample_data(SAMPLE)
  
  # Return the phyloseq object
   phyloseq(OTU, TAX, SAMPLE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Three file paths and sampleID column in metadata table required: 
        OTU taxonomic table, OTU sample table, metadata table, sampleID column in metadata and output_file_path.
        Check metadata format and give me sample ID column.
        \n Output file name is phyloseq_object ")
}

phyloseq_object <- phyloseq_creation(args[1], args[2], args[3], as.integer(args[4]))
# ps_test <- phyloseq_creation("/data2/leehg/OMD3/Taxonomy_interaction/output.tsv" , 
#                              "/data2/leehg/OMD3/Taxonomy_interaction/BrackentaxID_GTDBtaxa_completed.csv",
#                              "/data2/leehg/OMD3/Taxonomy_interaction/PRJEB10878_metadata.csv",3 )


output_file_path <- file.path(args[5], "phyloseq_object.RDS")
saveRDS(phyloseq_object, file = output_file_path)