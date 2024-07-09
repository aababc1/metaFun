library(shiny)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(InteractiveComplexHeatmap)
library(ragg)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
    cat("No arguments provided. Exiting script.\n -i : Input is skani full matrix, -m metadata CSV file, -out output heatmap file, -meta_cols : metadata column numbers like 2,3,4\n")
  quit(save = "no", status = 1)
}
arg_list <- list()
for (i in seq(1, length(args), 2)) {
  arg_list[[args[i]]] <- args[i + 1]
}

input_csv <- arg_list[["-i"]]
metadata_file <- arg_list[["-m"]]
output_file <- arg_list[["-out"]]
meta_col_selected <- arg_list[["-mc"]]
meta_col_selected <- as.integer(meta_col_selected)
meta_cols <- NULL
output_html <- arg_list[["-out_html"]]

if ("-meta_cols" %in% names(arg_list)) {
    meta_cols <- unlist(strsplit(arg_list[["-meta_cols"]], ","))
}


skani_raw <- read.csv(input_csv, header=FALSE,sep='\t',check.names = FALSE)
rownames(skani_raw) <- skani_raw[,1]
skani_raw <- skani_raw[,-1]

metadata <- read.csv(metadata_file, header = TRUE,check.names = FALSE)
selected_colname = colnames(metadata)[meta_col_selected]

if (!is.null(meta_cols)) {
    meta_cols <- as.numeric(meta_cols)
    metadata <- metadata[, meta_cols, drop = FALSE]
}

rownames(metadata) <- metadata[,1]
#skani_heatmap_data <- skani_raw[,2:(ncol(skani_raw))]

if(!all(rownames(skani_raw) %in% rownames(metadata))) {
    stop("Row names of input matrix and metadata do not match.")
}

selected_columns <- c("Completeness", "Contamination", "classification", colnames(metadata[selected_colname])) 
meta_filtered_forvis <- metadata[, selected_columns]
meta_filtered_forvis <- meta_filtered_forvis %>%
  mutate(!!selected_colname := replace(.[[selected_colname]], .[[selected_colname]] == "", "NA"))

classification_split <- strsplit(as.character(meta_filtered_forvis$classification), ";")
classification_df <- do.call(rbind, classification_split)
colnames(classification_df) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
classification_df <- as.data.frame(classification_df)
meta_filtered_forvis <- cbind(meta_filtered_forvis[, -which(colnames(meta_filtered_forvis) == "classification")], classification_df)

meta_filtered_forvis_ordered <- meta_filtered_forvis[match(rownames(skani_raw), rownames(meta_filtered_forvis)),]

top_annotation <- HeatmapAnnotation(df = meta_filtered_forvis_ordered,
                                    annotation_name_gp = gpar(fontsize = 7))

breaks <- c(0, min(skani_raw[skani_raw > 0]), max(skani_raw))
colors <- c("lightgrey", "white", "red")

if (length(unique(as.vector(skani_raw))) == 1) {
  my_color_function <- "red"
} else {
  my_color_function <- colorRamp2(breaks,colors)
}


ht <- Heatmap(
    as.matrix(skani_raw),
    name = "ANI values inferred using skani",
    #col = colorRampPalette(c("white", "red"))(255),
    col = my_color_function,
    top_annotation = top_annotation,
    show_row_names = FALSE,
    show_column_names = FALSE#,
    #column_split  = meta_filtered_forvis_ordered[,selected_colname]
)
output_file <- "heatmap_skani.pdf"

pdf(file = output_file, width = 15, height = 10)
draw(ht)
dev.off()

htShiny(ht, output_ui_float = TRUE, save = output_html)

