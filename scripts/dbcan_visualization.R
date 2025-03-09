#  parse command line arguments
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(ragg)
library(shiny)
library(InteractiveComplexHeatmap)

args <- commandArgs(trailingOnly = TRUE)

set_font <- function() {
  par(family = "Courier")
}
set_font()


if (length(args) == 0) {
    cat("No arguments provided. Exiting script.\n -i : Input is dbcan result file. -m metadata CSV file, -out output heatmap file, -meta_cols : metadata column numbers like 2,3,4\n")
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
# use this . get col number
#  group_col_name <- colnames(meco_ps$sample_table)[as.integer(metadata_col_index)]

# prohibit too long legend text. 
wrap_after_two_semicolons <- function(s) {
  parts <- unlist(strsplit(s, split = ";", fixed = TRUE))
  new_text <- ""
  for (i in seq_along(parts)) {
    new_text <- paste(new_text, parts[i], sep = "")
    if (i < length(parts) && (i %% 2 == 0)) {
      new_text <- paste(new_text, ";", "\n", sep = "")
    } else if (i < length(parts)) {
      new_text <- paste(new_text, ";", sep = "")
    }
  }
  new_text
}

clean_column_names <- function(df) {
  names(df) <- gsub(" ", "_", names(df))
  return(df)
}

dbcan_raw <- read.csv(input_csv, header = TRUE,check.names=FALSE)
dbcan_raw<- clean_column_names(dbcan_raw)

#input should be csv. 
metadata <- read.csv(metadata_file, header = TRUE, check.names = FALSE)
selected_colname = colnames(metadata)[meta_col_selected]

data <- read.csv(input_csv, header = TRUE, check.names = FALSE)
if (!is.null(meta_cols)) {
    meta_cols <- as.numeric(meta_cols)
    metadata <- metadata[, meta_cols, drop = FALSE]
}

#####################################
############# set this ##############
# set metadata rownames(the output should be grepped by genome_analysis_accession)
# genome_analysis_accession
# manually create based on scritps and get them . 
############# set this ##############
#####################################
rownames(metadata) <- metadata[,1] # should revise automated metadata  
# free below line 
# rownames(metadata) <- metadata[,colnaems(metadata[,])]

sort_labels <- function(df, label_col) {
  df %>%
    mutate(
      Prefix = str_extract(!!sym(label_col), "^[A-Z]+"),
      Number = as.numeric(str_extract(!!sym(label_col), "[0-9]+$")),
      Prefix_Order = match(Prefix, c("AA", "CBM", "CE", "GH", "GT", "PL")) # Set custom order of prefixes
    ) %>%
    arrange(Prefix_Order, Number) %>%
    select(-c(Prefix, Number, Prefix_Order))
}

dbcan_raw <- sort_labels(dbcan_raw,"HMMER")

head(dbcan_raw)
dbcan_raw <- dbcan_raw %>% 
  mutate(Prefix = str_extract(!!sym("HMMER"), "^[A-Z]+"))
dbcan_heatmap_data<- dbcan_raw[,2:(ncol(dbcan_raw)-1)]

selected_columns <- c("Completeness", "Contamination", "classification", colnames(metadata[selected_colname])) #meta1380 %>% select(ncol(.)))
meta_filtered_forvis <- metadata[, selected_columns]
meta_filtered_forvis <- meta_filtered_forvis %>%
  mutate(!!selected_colname := replace(.[[selected_colname]], .[[selected_colname]] == "", "NA"))

classification_split <- strsplit(as.character(meta_filtered_forvis$classification), ";")
classification_df <- do.call(rbind, classification_split)
colnames(classification_df) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
classification_df <- as.data.frame(classification_df)
meta_filtered_forvis <- cbind(meta_filtered_forvis[, -which(colnames(meta_filtered_forvis) == "classification")], classification_df)

meta_filtered_forvis_ordered <- meta_filtered_forvis[match(colnames(dbcan_heatmap_data), rownames(meta_filtered_forvis)),]

# top_annotation <- HeatmapAnnotation(df = meta_filtered_forvis_ordered,
#                                     annotation_name_gp = gpar(fontsize = 7))

basic_metadata <- meta_filtered_forvis_ordered[, !colnames(meta_filtered_forvis_ordered) %in% c(selected_colname,"Domain", "Phylum", "Class", "Order", "Family")]
user_metadata <- meta_filtered_forvis_ordered[, selected_colname, drop = FALSE]

top_annotation_boxplot <- HeatmapAnnotation(Average_value = anno_boxplot(as.matrix(dbcan_heatmap_data), axis = TRUE, height = unit(1.0, "cm"),
                               axis_param = list(labels_rot = 0)),
  annotation_name_side = "left", annotation_name_rot = 0,
  annotation_name_gp = gpar(fontsize = 11, fontfamily = "Courier", fontface="bold")
)

top_annotation_basic <- HeatmapAnnotation(
  df = basic_metadata,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 11, fontfamily = "Courier", fontface="bold"),
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 12, fontfamily = "Courier", fontface = "bold"),
    labels_gp = gpar(fontsize = 11, fontfamily = "Courier")
  )
)

top_annotation_user <- HeatmapAnnotation(
  df = user_metadata,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 11, fontfamily = "Courier", col = "red", fontface="bold"),
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 12, fontfamily = "Courier", fontface = "bold"),
    labels_gp = gpar(fontsize = 11, fontfamily = "Courier")
  )
)

integrated_data <- cbind(dbcan_heatmap_data, dbcan_raw[, c("Prefix", "HMMER")]) %>%

  mutate(rowsum = rowSums(select(., 1:(ncol(dbcan_heatmap_data)))))

integrated_data_sorted <- integrated_data %>%
  group_by(Prefix) %>%
  arrange(Prefix, desc(rowsum)) %>%
  mutate(row_index = row_number()) %>%
  ungroup() %>%
  arrange(Prefix, row_index)

row_order <- order(integrated_data_sorted$Prefix, integrated_data_sorted$row_index)

#row_order <- integrated_data_sorted$row_index
split_data <- split(as.data.frame(t(integrated_data_sorted[, 1:ncol(dbcan_heatmap_data)])), meta_filtered_forvis_ordered[, selected_colname])
split_matrix <- do.call(cbind, lapply(split_data, function(x) colMeans(x, na.rm = TRUE)))

distribution_colors <- scales::hue_pal()(ncol(split_matrix))
distribution_legend <- Legend(
    labels = colnames(split_matrix),
    legend_gp = gpar(fill = distribution_colors),
    title = "CAZyme_averages",
    title_gp = gpar(fontsize = 12, fontfamily = "Courier", fontface="bold"),
    labels_gp = gpar(fontsize = 11, fontfamily = "Courier")    
)

left_annotation <- rowAnnotation(
  CAZyme_averages = anno_points(
    split_matrix,
    pch = 1:ncol(split_matrix),
    size = unit(1, "mm"),
    title_gp = gpar(fontsize = 12, fontfamily = "Courier", fontface = "bold"),
    gp = gpar(col = scales::hue_pal()(ncol(split_matrix))),
    width = unit(3, "cm")
  ),
  annotation_name_gp = gpar(fontsize = 12, fontfamily = "Courier", fontface = "bold")
)

caz_info <- integrated_data_sorted[, c("Prefix", "HMMER")]
names(caz_info) <- c("CAZymes family", "CAZymes subfamily")

row_annotation <- HeatmapAnnotation(df = caz_info, which='row',
                                                                    show_legend = c("CAZymes family" = TRUE, "CAZymes subfamily" = FALSE),
  annotation_name_gp = gpar(fontsize = 11, fontfamily = "Courier", fontface = "bold"),
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 12, fontfamily = "Courier", fontface = "bold"),
    labels_gp = gpar(fontsize = 11, fontfamily = "Courier")
  )
)
message("which")
head(row_annotation)

length(unique(as.vector(dbcan_heatmap_data)))

breaks <- c(0, min(dbcan_heatmap_data[dbcan_heatmap_data>0]), max(dbcan_heatmap_data))
colors <- c("lightgrey", "white", "#22d62b")

if (length(unique(as.vector(dbcan_heatmap_data))) == 1) {
  my_color_function <- "#22d62b"
} else {
  my_color_function <- colorRamp2(breaks,colors)
}

ht <- Heatmap(
    as.matrix(integrated_data_sorted[, 1:ncol(dbcan_heatmap_data)]),
    name = "Carbohydrate-active enzymes",
    col = my_color_function,
    top_annotation = c(top_annotation_boxplot, top_annotation_basic, top_annotation_user),
    right_annotation = row_annotation,
    left_annotation = left_annotation,
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_column_dend = FALSE,
    show_row_dend = FALSE,
    cluster_rows = FALSE,
    column_split  = meta_filtered_forvis_ordered[,selected_colname],
    heatmap_legend_param = list(
        legend_direction = "horizontal",
        title_gp = gpar(fontsize = 12, fontfamily = "Courier", fontface="bold"),
        labels_gp = gpar(fontsize = 11, fontfamily = "Courier")
    ),
    row_names_gp = gpar(fontsize = 11, fontfamily = "Courier"),
    column_names_gp = gpar(fontsize = 11, fontfamily = "Courier"),
    column_title_gp = gpar(fontsize = 12, fontfamily = "Courier", fontface = "bold", col = "blue")
)

pdf_width <- (2.5*ncol(dbcan_heatmap_data) + 90 + 110 + 40) / 25.4
pdf_height <- (2.5*nrow(dbcan_heatmap_data) + 45 + 40 + 40) / 25.4
if (pdf_width >= 20.9) {
  pdf_width = 20.9
}


if (pdf_height >= 15) {
  pdf_height = 15.65
}

output_file <- "heatmap_dbCAN.pdf"
ht <- draw(ht,
    padding = unit(c(2, 2, 2, 2), "cm"),
    annotation_legend_list = list(distribution_legend),
    annotation_legend_side = "left"
)
pdf(file = output_file, width = pdf_width, height = pdf_height)
draw(ht, padding = unit(c(2, 2, 2, 2), "cm"), annotation_legend_list = list(distribution_legend), annotation_legend_side = "left")
#pdf(file = output_file)

dev.off()

htShiny(ht, output_ui_float = TRUE, save = output_html)
