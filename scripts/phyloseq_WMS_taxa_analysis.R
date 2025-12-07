library(microeco)
library(phyloseq)
library(file2meco)
library(ggplot2)
library(aplot)
library(magrittr)
library(gridExtra)
library(patchwork)
library(Cairo)
library(paletteer)
library(cowplot)
library(grid)
library(igraph)


run_analysis <- function(phyloseq_path, output_dir, metadata_col_index) {
  # load phyloseq object  and convert format 
  message("Loading the phyloseq object...")
  ps_obj <- readRDS(phyloseq_path)

  meco_ps <- phyloseq2meco(ps_obj)
  meco_ps$tidy_dataset()
  group_col_name <- colnames(meco_ps$sample_table)[as.integer(metadata_col_index)]
  print(group_col_name)
  # Alpha beta diversity calculation and visualize it.
  theme_set(theme_bw())
  meco_ps$cal_abund()
  meco_ps$cal_alphadiv(measures =c("Shannon","Simpson") ,PD = FALSE)
  # alpha_diversity.csv is saved into directory # make dirpath as argument 
  meco_ps$save_alphadiv(dirpath = output_dir)
  # calculate beta div and save : bray and jaccard 
  meco_ps$cal_betadiv(unifrac = FALSE)
  meco_ps$save_betadiv(dirpath = output_dir)
  alpha_meco <- trans_alpha$new(dataset = meco_ps, 
                                group = group_col_name)
                                #colnames(meco_ps$sample_table)[17])  
  alpha_meco$cal_diff(method = "KW")

  
  write.csv(alpha_meco$res_diff, file.path(output_dir, "alpha_group_diff_stat.csv"),
            row.names = FALSE)
  message("alpha_group_diff_stat.csv is saved.")
  #This action is for downstream visualization , prohibiting errors
  alpha_meco$res_diff %<>% base::subset(Significance != "ns")

  #Beta diversity calculation ,distance, visualization 
  abdiv_fig1 <- alpha_meco$plot_alpha(color_values = paletteer_d("ggsci::default_aaas"), 
                                      point_size = 2.0, point_alpha = 0.7,
                                      add_sig_text_size = 3)
  abdiv_fig1 <-  abdiv_fig1+  theme(axis.text.x = element_text(size = 6),
                      axis.text.y = element_text(size = 6),
                      axis.title.y = element_text(size = 7)) +
                scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  message("Alpha and beta diversity figure is created.")
  #Set group as argument 
  beta_meco <- trans_beta$new(dataset = meco_ps, 
                              group = group_col_name,
                              measure="bray")
  # group distance calculation and generate  barplot per group 
  beta_meco$cal_group_distance(within_group = TRUE)
  beta_meco$res_group_distance
  write.csv(beta_meco$res_group_distance,
            file.path(output_dir, "beta_group_distance_raw.csv"))
  beta_meco$cal_group_distance_diff(method = "wilcox")
  beta_meco$res_group_distance_diff
  write.csv(beta_meco$res_group_distance_diff,
            file.path(output_dir, "beta_group_distance_diff_stat.csv"))
  #Intra group, bray curtis distance visualization 
  bcd_plot <- beta_meco$plot_group_distance(boxplot_add = "mean",
                               color_values = paletteer_d("ggsci::default_aaas"),
                               add_sig_text_size = 3) + 
  theme(plot.margin = margin(0.5, 0.3, 0.3, 0.3, "cm")) 

  bcd_plot <-  bcd_plot+ theme(axis.text.x = element_text(size = 6),
                              axis.text.y = element_text(size = 6),
                              axis.title.y =  element_text(size = 7)) +
                                scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  message("beta diversity figure is created.")
  #permanova calculation
  beta_meco$cal_manova(manova_all = TRUE,p_adjust_method = "fdr")
  write.csv(beta_meco$res_manova,
            file.path(output_dir, "beta_group_perMANOVA_stat.csv"))
  R2_value <- beta_meco$res_manova$R2[1]
  Pseudo_F <- beta_meco$res_manova$F[1]
  PrF_value <- beta_meco$res_manova$'Pr(>F)'[1]
  R2_text <- paste("R² =", format(R2_value, digits = 4))
  Pseudo_F_text <- paste("Pseudo-F statistic =", format(Pseudo_F,digits=4))
  PrF_text <- paste("Pr(>F) =", format(PrF_value, digits = 4))
  beta_meco$cal_ordination(ordination ="PCoA")
  beta_p1 <- beta_meco$plot_ordination(plot_color = group_col_name, 
                                      plot_shape = group_col_name,
                                      point_size = 1.0,point_alpha = 0.7,
                                      plot_type = c("point", "ellipse"),
                                      color_values = paletteer_d("ggsci::default_aaas"))
  beta_p1 <- beta_p1 + 
    annotate("text", x = Inf, y = Inf, label = R2_text, hjust = 1.1, vjust = 5.5, size = 1.7) +
    annotate("text", x = Inf, y = Inf, label = Pseudo_F_text, hjust = 1.1, vjust = 3.5, size = 1.7) +
    annotate("text", x = Inf, y = Inf, label = PrF_text, hjust = 1.1, vjust = 1.5, size = 1.7)+
   theme(legend.position = "none", 
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y =  element_text(size = 7),
        axis.title.x =  element_text(size = 7))         
  beta_p1
  # generate x axis and y axis plots and comine PCoA plot.
  pcoa_axis <-  beta_meco$res_ordination$scores
  pcoa_te1 <- trans_env$new(dataset = meco_ps, add_data = pcoa_axis[, 1:2])
  pcoa_te1$cal_diff(group =  group_col_name,
                    method = "KW")

  # groups order in p2 is same with p1; use legend.position = "none" to remove redundant legend
  beta_p2 <- pcoa_te1$plot_diff(measure = "PCo1", add_sig = T,
                              point_size = 2,point_alpha = 0.7,
                              color_values = paletteer_d("ggsci::default_aaas")) + 
  theme_bw() + coord_flip() + theme(legend.position = "none", 
                                    axis.title.x = element_blank(), 
                                    axis.title.y = element_blank(),                                    
                                    axis.text.x  = element_blank(), 
                                    axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank())

  beta_p3 <- pcoa_te1$plot_diff(measure = "PCo2", add_sig = T,
                              color_values = paletteer_d("ggsci::default_aaas"),
                              point_size = 1.0,point_alpha = 0.7) +
  theme_bw() + theme(legend.position = "none", 
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     axis.text.x = element_blank(), 
                     axis.text.y  = element_blank(), 
                     axis.ticks.x = element_blank())

beta_fig <- beta_p1 %>% insert_top(beta_p2, height = 0.2) %>% 
  insert_right(beta_p3, width = 0.2) 


legend <- cowplot::get_legend(beta_p1 + theme(legend.position = "bottom",
                                              strip.text = element_text(size = 7),
                                              legend.title = element_text(size = 7),
                                              legend.text = element_text(size = 6),
                                              legend.key.size = unit(0.3, "cm")))

  # height of the upper figure and width of the right-hand figure are both 0.2-fold of the main figure
  # this is the result fig 
  # convert aplot class (beta div figure) to combine 
  beta_fig_grob <- grid::grid.draw(beta_fig)
  beta_fig_grob <- grid::grid.grab()

layout <- rbind(c(1, 2),
                c(3),
                c(3),
                c(4))

# Adjust the top margin of the second plot
abdiv_fig1 <- abdiv_fig1 + theme(plot.margin = margin(t = 1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))
bcd_plot <- bcd_plot + theme(plot.margin = margin(t = 1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))

combined_diversity <- gridExtra::grid.arrange(abdiv_fig1,bcd_plot, beta_fig_grob,
                                              legend,
                                              layout_matrix = layout,
                                              heights = unit(c(1.2, 1.2, 2, 0.1), "null"),
                                              widths = unit(c(1, 1),"null"))

                                                 
  CairoPDF(file = file.path(output_dir, "Taxonomy_alpha_beta_diversity.pdf"), 
         width = 100 / 25.4, 
         height = 130  / 25.4
)
grid.draw(combined_diversity)
dev.off()
message("Alpha and beta diversity figure was saved.")     

tryCatch({
  # differential analysis 
  
  meco_ps$filter_taxa(rel_abund = 0.0001, freq=3) 
  lefse_group <- trans_diff$new(meco_ps, method = "lefse", 
                      group = group_col_name, taxa_level = "Species",
                      alpha = 0.01,p_adjust_method = "none")

  glefse_1 <- lefse_group$plot_diff_bar(use_number = 1:20,
  color_values = paletteer_d("ggsci::default_aaas"))

  glefse_1 <- glefse_1 + theme(legend.position = "none", 
        axis.text.y = element_text(size = 7),
        axis.text.x =element_text(size = 6) ,
        axis.title.x = element_text(size = 7)) 

  glefse_2 <- lefse_group$plot_diff_abund(use_number = 1:20,
    color_values = paletteer_d("ggsci::default_aaas"), 
              add_sig = TRUE,select_taxa = lefse_group$plot_diff_bar_taxa)
  
glefse_2 <- glefse_2 +  theme(axis.text.y = element_blank(), 
                            axis.ticks.y = element_blank(),
                            axis.text.x =element_text(size = 6) ,
                            axis.title.x = element_text(size = 7),
                            legend.key.size = unit(0.3, "cm"),
                            legend.title = element_text(size = 7),
                            legend.text = element_text(size = 6.5),
                            strip.text = element_text(size = 6)) + 
  geom_vline(xintercept = -Inf, color = "gray", size = 0.5) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

  glefse_1 + glefse_2

  lefse_plot <- gridExtra::grid.arrange(glefse_1,glefse_2, ncol = 2 , nrow = 1 , widths = c(2,1.3)) 
  CairoPDF(file = file.path(output_dir, "Taxonomy_lefse_result.pdf"), 
         width = 20/2.54, 
         height = (0.77*length(lefse_group$plot_diff_bar_taxa)+1)/2.54
)
grid.draw(lefse_plot)
dev.off()
message("Diefferential taxa analysis by lefse : figure is saved.")    
}, error = function(e) {
  message("Differential analysis failed with error: ", e$message)

})


tryCatch({

  #DA heatmap , order is same with DA plots. 
  DA_species_name <- gsub("^s__", "", lefse_group$plot_diff_bar_taxa)
  tt <- trans_abund$new(meco_ps,taxrank = "Species" ,input_taxaname =DA_species_name )
  glefse_3 <- tt$plot_heatmap(facet = group_col_name, xtext_keep = FALSE, 
                withmargin = FALSE,plot_colorscale = "log10",
                min_abundance = 0.0001,,plot_text_size = 6,
                  ytext_size = 6,strip_text = 7) + 
  ggtitle("DA species in relative abundance(%)") +
    theme(title = element_text(size=7.5),
          legend.title = element_text(size=7),
          legend.text  = element_text(size=7),
          legend.key.height = unit(0.2, 'inches'))
  
  CairoPDF(file = file.path(output_dir, "Taxonomy_lefse_heatmap.pdf"),
           width = 20/2.54, 
           height = (0.35*length(lefse_group$plot_diff_bar_taxa)+2)/2.54
  )
  grid.draw(glefse_3)
  dev.off()

  message("Diefferential taxa heatmap by lefse : figure was saved.")    
  # numerical metadata correlation 
  message("calcluating numerical cateogries correlation.")    

  te1 <- trans_env$new(meco_ps, env_cols =c(1:ncol(meco_ps$sample_table)))
    numerical_metadata_stats <- te1$cal_autocor(group_col_name,
    color_values = paletteer_d("ggsci::default_aaas")) +
  theme(strip.text = element_text(size = 7),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        text = element_text(size = 6))

CairoPDF(file = file.path(output_dir, "Taxonomy_metadata_autocorrelation.pdf"),
         height = (sum(sapply(te1$data_env, is.numeric))+1)*4/2.54, 
         width = (sum(sapply(te1$data_env, is.numeric))+1)*4/2.54
)
grid.draw(numerical_metadata_stats)
dev.off()
message("Calcluating numerical cateogries correlation info was saved.")    
}, error = function(e) {
  message("Differential analysis heatmap generation error: ", e$message)
})

tryCatch({
  message("Aligning lefse plot and relative abundance plot .")    
  species_names <- lefse_group$plot_diff_bar_taxa
  full_taxa <- lefse_group$res_diff$Taxa
  matched_full_taxa <- vector("character", length(species_names))
  for (i in seq_along(species_names)) {
  species_pattern <- paste0(".*", species_names[i], "$")  # Create the pattern
  matched_index <- grep(species_pattern, full_taxa)      # Find the index of the match
  matched_full_taxa[i] <- ifelse(length(matched_index) > 0, full_taxa[matched_index], NA)  # Store the matched full taxonomy
  }
  te1$cal_cor(by_group = group_col_name, use_data = "other", 
            p_adjust_method = "fdr", other_taxa = matched_full_taxa)
  
  message("Correlation betwwen Numerical metadata variable and differentially abundant taxa.")    
  message("Saved into !!! Taxonomy_species_meta_correlation.pdf !!!")
  species_meta_correlation <-te1$plot_cor(xtext_size = 6) + 
  theme(plot.margin = margin(b = 0.1, l = 1.5, unit = "cm"),
        axis.text.y = element_text(size=6),
        strip.text = element_text(size=7),
        legend.title = element_text(size=7),
        legend.text = element_text(size=6))

  #ggsave(file.path(output_dir, "species_meta_correlation.pdf"), 
  #        species_meta_correlation, width = 5.5, height= 6)


CairoPDF(file = file.path(output_dir, "Taxonomy_DA_meta_correlation.pdf"),
         height = (length(unique(te1$res_cor$Taxa))*0.35+3)/2.54, 
         width = (length(unique(te1$res_cor$Env))*
                    length(unique(te1$res_cor$Type))*0.35 +
                    (length(unique(te1$res_cor$Type))-1)*0.35 +
                    10)/2.54
)
grid.draw(species_meta_correlation)
dev.off()
}, error = function(e) {
  message("lefse_pCoA combined plot generation error: ", e$message)

})

tryCatch({

t_net1 <- trans_network$new(dataset = meco_ps, cor_method = "sparcc", 
                        use_sparcc_method = "SpiecEasi",  taxa_level = "Species",
                        filter_thres = 0.003)

t_net1$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb")
t_net1$res_network
t_net1$cal_module(method = "cluster_fast_greedy")
#
t_net1$save_network(filepath = file.path(output_dir,"network.gexf"))
net_fig <- t_net1$plot_network(method = "ggraph", node_color = "Phylum")
    ggsave(file.path(output_dir, "Taxonomy_species_correlation_network.pdf"), 
          net_fig, width = 5.5, height= 6)
  message("Saved into !!! Taxonomy_species_meta_correlation.pdf !!!")
  message("All finished")
}, error = function(e) {
  message("Network generation error: ", e$message)

})

}


args <- commandArgs(trailingOnly = TRUE)
phyloseq_path <- args[1]
output_dir <- args[2]
metadata_col_index <- args[3]  # This is the new argument for the metadata column index

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

run_analysis(phyloseq_path, output_dir, metadata_col_index)
