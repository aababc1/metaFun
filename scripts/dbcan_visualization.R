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


# arg_info <- data[, (ncol(data)-1):ncol(data)]
# headmap_data <- data[,(1:ncol(data)-2)]
# top_annotation <- HeatmapAnnotation(df = metadata)
# row_annotation <- rowAnnotation(
#   CARD_Short_Name = anno_text(resistance_mechanism$CARD_Short_Name),
#   Resistance_Mechanism = anno_text(resistance_mechanism$Resistance_Mechanism)
# )




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

top_annotation <- HeatmapAnnotation(df = meta_filtered_forvis_ordered,
                                    annotation_name_gp = gpar(fontsize = 7))

# row_annotation <- rowAnnotation(
#     df = arg_info,
#     annotation_name_gp = gpar(fontsize = 8)
# )

#arg_info <- card_raw[, (ncol(card_raw) - 2):ncol(card_raw)]

caz_info <- dbcan_raw[,c(1,ncol(dbcan_raw))]
head(caz_info)
caz_info <- caz_info[, c("Prefix","HMMER")]
names(caz_info) <- c("CAZymes family","CAZymes subfamily")


row_annotation <- HeatmapAnnotation(df = caz_info, which='row',
                                                                    show_legend = c("CAZymes family" = TRUE, "CAZymes subfamily" = FALSE),
                                                                    annotation_name_gp = gpar(fontsize = 7))
message("which")
head(row_annotation)
#annotation_width = unit(4, "cm")

length(unique(as.vector(dbcan_heatmap_data)))

breaks <- c(0, min(dbcan_heatmap_data[dbcan_heatmap_data>0]), max(dbcan_heatmap_data))
colors <- c("lightgrey", "white", "#22d62b")

if (length(unique(as.vector(dbcan_heatmap_data))) == 1) {
  my_color_function <- "#22d62b"
} else {
  my_color_function <- colorRamp2(breaks,colors)
}


ht <- Heatmap(
    as.matrix(dbcan_heatmap_data),
    name = "Carbohydrate-active enzymes",
    col = my_color_function,
    top_annotation = top_annotation,
    right_annotation = row_annotation,
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_column_dend = FALSE,
    show_row_dend = FALSE,
    cluster_rows = FALSE,
    column_split  = meta_filtered_forvis_ordered[,selected_colname]
)

output_file <- "heatmap_dbCAN.pdf"

pdf(file = output_file, width = 15, height = 10)
#pdf(file = output_file)
draw(ht)
dev.off()






htShiny(ht, output_ui_float = TRUE, save = output_html)
