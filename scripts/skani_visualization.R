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

set_font <- function() {
  par(family = "Courier")
}
set_font()

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

print(head(rownames(skani_raw)))
if (anyDuplicated(rownames(skani_raw))) {
    stop("Duplicate row names found in skani_fullmatrix.")
}

metadata <- read.csv(metadata_file, header = TRUE,check.names = FALSE)
#metadata <- read.csv(metadata_file, header = TRUE, check.names = FALSE)
head(metadata)

print(head(rownames(metadata)))
if (anyDuplicated(rownames(metadata))) {
    duplicates <- rownames(metadata)[duplicated(rownames(metadata))]
    stop(paste("Duplicate row names found in the metadata file:", paste(duplicates, collapse = ", ")))
}

if (anyDuplicated(metadata[,1])) {
    stop("Duplicate row names found in the metadata file.")
}


selected_colname = colnames(metadata)[meta_col_selected]

if (!is.null(meta_cols)) {
    meta_cols <- as.numeric(meta_cols)
    metadata <- metadata[, meta_cols, drop = FALSE]
}

rownames(metadata) <- metadata[,1]
print(head(rownames(metadata)))


#skani_heatmap_data <- skani_raw[,2:(ncol(skani_raw))]

print(head(metadata))
print(dim(metadata))

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

meta_filtered_forvis <- meta_filtered_forvis[, !colnames(meta_filtered_forvis) %in% c("classification","Domain", "Phylum", "Class", "Order", "Family")]

#meta_filtered_forvis <- cbind(meta_filtered_forvis[, -which(colnames(meta_filtered_forvis) %in% c("Domain", "Phylum", "Class", "Order", "Family"))], classification_df)
#meta_filtered_forvis <- cbind(meta_filtered_forvis[, -which(colnames(meta_filtered_forvis) == "classification")], classification_df)

#selected_columns <- c("Completeness", "Contamination", selected_colname, "Genus", "Species")
#meta_filtered_forvis <- meta_filtered_forvis[, selected_columns]


meta_filtered_forvis_ordered <- meta_filtered_forvis[match(rownames(skani_raw), rownames(meta_filtered_forvis)),]

top_annotation <- HeatmapAnnotation(df = meta_filtered_forvis_ordered,
 annotation_name_gp = gpar(fontsize = 11, fontfamily = "Courier",  fontface="bold"),

  annotation_legend_param = list(
    title_gp = gpar(fontsize = 12, fontfamily = "Courier", fontface = "bold"),
    labels_gp = gpar(fontsize = 11, fontfamily = "Courier")
  )
)

# breaks <- c(0, min(skani_raw[skani_raw > 0]), max(skani_raw))
# colors <- c("lightgrey", "white", "red")

# if (length(unique(as.vector(skani_raw))) == 1) {
#   my_color_function <- "red"
# } else {
#   my_color_function <- colorRamp2(breaks,colors)
# }

 breaks <- c(80, 80.000001, 100)
  colors <- c("lightgrey", "white", "red")
  my_color_function <- colorRamp2(breaks, colors)


ht <- Heatmap(
    as.matrix(skani_raw),
    name = "ANI values inferred using skani",
    #col = colorRampPalette(c("white", "red"))(255),
    col = my_color_function,
    top_annotation = top_annotation,
    show_row_names = FALSE,
    show_column_names = FALSE,#,
    #column_split  = meta_filtered_forvis_ordered[,selected_colname]
    heatmap_legend_param = list(
      at = c(80, 90, 100),
      labels = c("80%", "90%", "100%"),            
      title_gp = gpar(fontsize = 12, fontfamily = "Courier", fontface = "bold"),
      labels_gp = gpar(fontsize = 11, fontfamily = "Courier")
    ),
    row_names_gp = gpar(fontsize = 11, fontfamily = "Courier"),
    column_names_gp = gpar(fontsize = 11, fontfamily = "Courier")
)
output_file <- "heatmap_skani.pdf"


pdf_width <- (2.5 * ncol(skani_raw) + 90 + 110 + 100) / 25.4
pdf_height <- (2.5 * nrow(skani_raw) + 45 + 40 + 40) / 25.4

if (pdf_width >= 20.9) {
  pdf_width = 20.9
}
if (pdf_height >= 15) {
  pdf_height = 15.65
}





pdf(file = output_file, width =pdf_width , height = pdf_height)
draw(ht)
dev.off()

htShiny(ht, output_ui_float = TRUE, save = output_html)

