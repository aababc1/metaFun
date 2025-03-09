library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(circlize)
library(InteractiveComplexHeatmap)
library(shiny)
library(DT)
library(htmlwidgets)
#  parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

set_font <- function() {
  par(family = "Courier")
}

set_font()

if (length(args) == 0) {
# Handle the case where no arguments are passed
#    Set default values or exit with a message
    cat("No arguments provided. Exiting script.\n -i : KO.csv matrix, -m metadata.csv -out_html shinyout -out_table module_table.tsv,  \n -mc : module completeness - optional  \n -meta_cols : column number like 2,3,4   ")
    quit(save = "no", status = 1)
    }

arg_list <- list()
for (i in seq(1, length(args), 2)) {
  arg_list[[args[i]]] <- args[i + 1]
}

# Assign arguments to variables
input_kegg_matrix <- arg_list[["-i"]]
# metadata <- arg_list[["-m"]]
output_html <- arg_list[["-out_html"]]
output_table <- arg_list[["-out_table"]]
metadata_file <- arg_list[["-m"]]

meta_col_selected <- arg_list[["-n"]]
meta_col_selected <- as.integer(meta_col_selected)
meta_cols <- NULL


#metadata_provided <- "-m" %in% names(arg_list)
#metadata <- NULL
#meta_cols <- NULL
# if ("-meta_cols" %in% names(arg_list)) {
#     meta_cols <- unlist(strsplit(arg_list[["-meta_cols"]], ","))
# }
metadata <- read.csv(metadata_file, header = TRUE, check.names = FALSE)
# default KEGG module threshold is 50%. 
module_completeness_threshold <- 0.5

selected_colname = colnames(metadata)[meta_col_selected]

if ("-mc" %in% names(arg_list)) {
  mc_value <- as.numeric(arg_list[["-mc"]])
  if (!is.na(mc_value) && mc_value >= 0 && mc_value <= 1) {
    module_completeness_threshold <- mc_value
  } else {
    warning("Invalid module completeness threshold provided. Using default value of 0.5 that is 50%.")
  }
}

library(ggkegg)
library(ggfx)
library(tidygraph)
library(dplyr)
library(tidyr)
library(parallel)

# Read the input files

KOtable <- read.csv(input_kegg_matrix, header = TRUE,check.names =FALSE)
KOtable <- KOtable[,-1]


ko_long_format <-  KOtable %>%
  pivot_longer(cols = everything(), names_to = "Genome", values_to = "KOs")
unique(ko_long_format$Genome)


mapper <- data.table::fread("/opt/database/kofam/mapper.tsv", header=TRUE)
allmods <-  gsub("md:", "", mapper$V1) %>%  unique()


# modules downloaded in 2023 Dec 20 
mf <-  list.files("/opt/database/KEGG_modules")
mf <- mf[startsWith(mf, "M")]
annos <- list()

library(parallel)
cl <-  makeCluster(32)
clusterEvalQ(cl, {
  library(ggkegg)
  library(ggfx)
  library(tidygraph)
  library(dplyr)
  library(tidyr)
})
calculate_completeness <-  function(genome) {
  kos <- KOtable[[genome]][!is.na(KOtable[[genome]])]
  mcs <- sapply(mf, function(mod) {
    mc <- module_completeness(module(mod, directory="/opt/database/KEGG_modules/"), 
                              query = kos)
    mean(mc$complete)
  })
  names(mcs) <- mf
  return(mcs)
}
clusterExport(cl, c("KOtable", "mf", "calculate_completeness")) 
results <- parLapply(cl, colnames(KOtable), calculate_completeness)
names(results) <- colnames(KOtable)
stopCluster(cl)

hdf <- data.frame(results,check.names=FALSE)

# for (genome in colnames(KOtable)) {
#   # get KOs in one genome 
#   kos <- KOtable[[genome]]
#   kos <- kos[!is.na(kos)]

#   mcs <-  NULL
#   for (mod in mf) {
#     mc <- module_completeness(module(mod, directory="/scratch/tools/microbiome_analysis/database/KEGG_modules"),
#                               query = kos)
#     mcs <-c(mcs, mc$complete |> mean())
#   }
#   annos[[as.character(genome)]] <- mcs 
# }
#hdf <- data.frame(annos, check.names=FALSE)

row.names(hdf) <- mf
hdf[is.na(hdf)] <- 0
threshold <- 0.5
hdf <- hdf[apply(hdf, 1, function(x) any(x >= module_completeness_threshold)), ]

# module description and CLASS description third columns are most bigger concept 
#moddesc <- data.table::fread("https://rest.kegg.jp/list/module", header=FALSE)
#data.table::fwrite(module_desc,"module_desc.txt",col.names=TRUE)
moddesc<- data.table::fread("/opt/database/kofam/module_desc.txt", header=FALSE)

mod_class <- read.csv("/opt/database/KEGG_modules/modes_class_description.tsv",
         sep='\t',header=FALSE)
#mod_class_name <- read.csv("/opt/database/KEGG_modules/modes_class_Name_description.tsv",
#         sep='\t',header=FALSE)


merged_modedesc <- merge(moddesc, mod_class, by="V1")

library(ragg)
library(ComplexHeatmap)
library(circlize)
library(InteractiveComplexHeatmap)
rownames(hdf)

#generate ordered module table following CLASS
merged_modedesc <- merged_modedesc %>%
  arrange(V2.y, V3, V2.x)
module_order <- merged_modedesc$V1[merged_modedesc$V1 %in% rownames(hdf)]
hdf <- hdf[module_order, ]
annotation_reordered_modedesc <- merged_modedesc %>%
  filter(V1 %in% rownames(hdf)) %>%
  select(V2.x, V3,V2.y)


library(tibble)
## new add, order 
hdf <- hdf %>%
  rownames_to_column("Module") %>%
  left_join(merged_modedesc, by = c("Module" = "V1")) %>%
  mutate(rowsum = rowSums(select(., where(is.numeric)))) %>%
  group_by(V2.y) %>%
  arrange(V2.y, desc(rowsum)) %>%
  ungroup() %>%
  column_to_rownames("Module") %>%
  select(-rowsum, -V2.x, -V3, -V2.y)

# Ensure the order is consistent across all annotations and heatmaps
# This was changed.
# annotation_reordered_modedesc <- annotation_reordered_modedesc %>%
#   arrange(V2.y, desc(rowsum))


annotation_reordered_modedesc <- merged_modedesc %>%
  filter(V1 %in% rownames(hdf)) %>%
  mutate(rowsum = rowSums(select(hdf, where(is.numeric)))) %>%
  arrange(V2.y, desc(rowsum)) %>%
  select(V2.x, V3, V2.y)

names(annotation_reordered_modedesc) <- c("Module_Name", "Pathway_Category", "Pathway_Class")

ha_module <- rowAnnotation(df = annotation_reordered_modedesc[, c("Pathway_Category", "Pathway_Class")],
                           gp = gpar(fontsize = 11, fontfamily = "Courier"),  
                           annotation_name_gp = gpar(fontsize = 12, fontfamily = "Courier", fontface = "bold"),
                           annotation_legend_param = list(
                             title_gp = gpar(fontsize = 12, fontfamily = "Courier", fontface = "bold"),
                             labels_gp = gpar(fontsize = 11, fontfamily = "Courier")
                           )
)

## new add, order 
## previous order
# names(annotation_reordered_modedesc) <- c("Module_Name", "Pathway_Category", "Pathway_Class")

# #    title_gp = gpar(fontsize = 12, fontfamily = "Courier",fontface="bold"),
#     labels_gp = gpar(fontsize = 11, fontfamily = "Courier")    

# ha_module <- rowAnnotation(df = annotation_reordered_modedesc[, c("Pathway_Category", "Pathway_Class")],
# gp = gpar(fontsize = 11, fontfamily = "Courier"),  
# annotation_name_gp = gpar(fontsize = 12, fontfamily = "Courier", fontface="bold"),
# annotation_legend_param = list(
#     title_gp = gpar(fontsize = 12, fontfamily = "Courier", fontface = "bold"),
#     labels_gp = gpar(fontsize = 11, fontfamily = "Courier")
#   )
# )
## previous order

#title_gp = gpar(fontsize = 12, fontfamily = "Courier",fontface="bold"))


# metadata 
# if (!is.null(metadata)) {
#     meta_selected <- metadata %>%  
#       select(country,continent, disease_group, age_group)
#     rownames(meta_selected) <- metadata[,"Name"]

#     meta_selected <- meta_selected[rownames(meta_selected)%in% colnames(hdf),]
#     meta_selected_ordered <- meta_selected[match(colnames(hdf), rownames(meta_selected)), ]

#     ha_meta <-  HeatmapAnnotation(df = meta_selected_ordered) 
# } else { 
#     ha_meta <- NULL 
# } 
#metadata parsing
rownames(metadata) <- metadata[,1]
selected_columns <- c("Completeness", "Contamination", "classification", colnames(metadata[selected_colname]))
meta_filtered_forvis <- metadata[, selected_columns]
head(meta_filtered_forvis)
meta_filtered_forvis <- meta_filtered_forvis %>%
  mutate(!!selected_colname := replace(.[[selected_colname]], .[[selected_colname]] == "", "NA"))

classification_split <- strsplit(as.character(meta_filtered_forvis$classification), ";")
classification_df <- do.call(rbind, classification_split)
colnames(classification_df) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
classification_df <- as.data.frame(classification_df)
meta_filtered_forvis <- cbind(meta_filtered_forvis[, -which(colnames(meta_filtered_forvis) == "classification")], classification_df)
meta_filtered_forvis_ordered <- meta_filtered_forvis[match(colnames(hdf), rownames(meta_filtered_forvis)),]


#split selected one and basic metadata
basic_metadata <- meta_filtered_forvis_ordered[, !colnames(meta_filtered_forvis_ordered) %in% c(selected_colname,"Domain", "Phylum", "Class", "Order", "Family")]
user_metadata <- meta_filtered_forvis_ordered[, selected_colname, drop = FALSE]

top_annotation_boxplot <- HeatmapAnnotation(Average_value = anno_boxplot(as.matrix(hdf), axis = TRUE,height = unit(1.0, "cm"),
                               axis_param = list(labels_rot = 0)),
  annotation_name_side = "left", annotation_name_rot = 0,
  annotation_name_gp = gpar(fontsize = 11, fontfamily = "Courier",fontface="bold")#,
  #height = unit(1.0, "cm")

)

top_annotation_basic <- HeatmapAnnotation(
  df = basic_metadata,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 11, fontfamily = "Courier",fontface="bold"),
  annotation_legend_param = list(
  title_gp = gpar(fontsize = 12, fontfamily = "Courier", fontface = "bold"),
  labels_gp = gpar(fontsize = 11, fontfamily = "Courier")
  ))

#  annotation_legend_param = list(
#    title_gp = gpar(fontsize = 12, fontfamily = "Courier", fontface = "bold"),
#    abels_gp = gpar(fontsize = 11, fontfamily = "Courier")

top_annotation_user <- HeatmapAnnotation(
  df = user_metadata,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 11, fontfamily = "Courier", col = "red", fontface="bold"),
    annotation_legend_param = list(
    title_gp = gpar(fontsize = 12, fontfamily = "Courier", fontface = "bold"),
    labels_gp = gpar(fontsize = 11, fontfamily = "Courier")
  )
)


#    annotation_legend_param = list(
#    title_gp = gpar(fontsize = 12, fontfamily = "Courier", fontface = "bold"),
 #   abels_gp = gpar(fontsize = 11, fontfamily = "Courier")
    #,
  #height = unit(1.0, "cm")

#top_annotation <- HeatmapAnnotation(Average_value = anno_boxplot(as.matrix(hdf), axis = TRUE,height = unit(1.2, "cm")),
#                          df = meta_filtered_forvis_ordered,
#                          annotation_name_side = "left",
#                          annotation_name_gp = gpar(fontsize =10,fontfamily = "Courier"))
#box_width = 0.6,

#left annotation #left annotation #left annotation 
#left annotation 
#left annotation 

split_data <- split(as.data.frame(t(hdf)), meta_filtered_forvis_ordered[, selected_colname])
split_matrix <- do.call(cbind, lapply(split_data, function(x) colMeans(x, na.rm = TRUE)))

print(paste("Dimensions of hdf:", paste(dim(hdf), collapse = " x ")))
print(paste("Dimensions of split_matrix:", paste(dim(split_matrix), collapse = " x ")))

distribution_colors <- scales::hue_pal()(ncol(split_matrix))
distribution_legend <- Legend(
    labels = colnames(split_matrix),
    legend_gp = gpar(fill = distribution_colors),
    title = "Module_averages",
    title_gp = gpar(fontsize = 12, fontfamily = "Courier",fontface="bold"),
    labels_gp = gpar(fontsize = 11, fontfamily = "Courier")    
)

print(paste("Number of groups after split:", length(split_data)))
print(paste("Dimensions of split_matrix:", paste(dim(split_matrix), collapse = " x ")))
print(paste("Dimensions of hdf:", paste(dim(hdf), collapse = " x ")))
#head(hdf)
#print(split_matrix)


if (nrow(split_matrix) != nrow(hdf)) {
  stop("Mismatch in dimensions between split_matrix and hdf.")
}

# left annotation , average module value in each  group of genomes 
left_annotation <- rowAnnotation(
  Module_averages = anno_points(split_matrix, 
                             pch = 1:ncol(split_matrix),
                             size = unit(1, "mm"),
                             title_gp = gpar(fontsize = 12, fontfamily = "Courier",fontface="bold"),

                             gp = gpar(col = scales::hue_pal()(ncol(split_matrix))),
                             width = unit(3, "cm")
  ),  annotation_name_gp = gpar(fontsize = 12, fontfamily = "Courier", fontface="bold")
)

#left annotation #left annotation #left annotation 

col_fun <- colorRamp2(c(0, module_completeness_threshold -0.0000001 ,module_completeness_threshold, 1), c("black","black", "white", scales::muted("red")))
head(hdf)

# cell_width <- unit(2.5, "mm")
# cell_height <- unit(2.5, "mm")
# left_annotation_width <- unit(5, "cm")  
#pdf_width <- (heatmap_width + left_annotation_width + unit(210, "mm") + margin * 2) / 25.4
#left 5cm + 4cm  , right 11cm , top 4.5cm  , margin 20 mm 
margin <- 20 
# margin *2 + bottom top , marig*2 + left right
pdf_width <- (2.5*ncol(hdf) +90 + 110 +  40) / 25.4
pdf_height <- (2.5*nrow(hdf) + 45 +40 + 40 ) / 25.4

# set maximum length
if (pdf_width >= 20.9) {
  pdf_width = 20.9
}
if (pdf_height >= 15) {
  pdf_height = 15.65
}

#width = unit(4, "cm")

ht1 <- Heatmap(as.matrix(hdf), show_column_names = FALSE,
               col = col_fun,
               right_annotation = ha_module,
               left_annotation=left_annotation,
#               top_annotation = top_annotation,
               top_annotation = c(top_annotation_boxplot,top_annotation_basic, top_annotation_user),
               show_column_dend = FALSE,
               cluster_rows = FALSE,
               heatmap_legend_param = list(
                legend_direction = "horizontal",
                title_gp = gpar(fontsize = 12, fontfamily = "Courier",fontface="bold"),
                labels_gp = gpar(fontsize = 11, fontfamily = "Courier")
                ), 
                name = "Module completeness",
                border = TRUE,
                column_split  = meta_filtered_forvis_ordered[,selected_colname],
                row_names_gp = gpar(fontsize = 11, fontfamily = "Courier"),
                column_names_gp = gpar(fontsize = 11, fontfamily = "Courier"),
               column_title_gp = gpar(fontsize = 12, fontfamily = "Courier", fontface = "bold", col = "blue"),  #set font of row title (metadata )

                #height = unit(2.5, "mm"),
                #width = unit(2.5, "mm")#,
                 )
#ht1 <- draw(ht1)
#packLegend(distribution_legend, just = c("left", "top"))
#ht1 <- ht1 + draw(distribution_legend, just = c("left", "top"))
ht1 <- draw(ht1,
padding = unit(c(2, 2, 2, 2), "cm"),
annotation_legend_list = list(distribution_legend),annotation_legend_side = "left")

#draw(lgd, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))

output_file <- "heatmap_KEGG.pdf"
pdf(file = output_file, width = pdf_width, height = pdf_height)
#pdf(file = output_file, width = 20, height = 15)
draw(ht1)
dev.off()

# Save the interactive plot
hover_info <- data.frame(
  Module = rownames(hdf),
  Name = annotation_reordered_modedesc$Module_Name,
  Category = annotation_reordered_modedesc$Pathway_Category,
  Class = annotation_reordered_modedesc$Pathway_Class
)

######################################################

######################################################
# htShiny기존 코드 . 

htShiny(ht1, output_ui_float = TRUE, save = output_html)

# Save the module completeness table
write.csv(hdf, file = output_table)
