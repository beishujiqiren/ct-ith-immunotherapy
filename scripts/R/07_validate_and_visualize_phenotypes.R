

library(vegan)   
library(Rtsne)   
library(pheatmap)

evaluate_separation <- function(data_matrix, labels, set_name) {
  
  
  centroids <- aggregate(data_matrix, by=list(Group=labels), FUN=mean)
  rownames(centroids) <- centroids$Group
  centroids <- centroids[,-1] 
  
  
  dist_val <- dist(rbind(centroids["High",], centroids["Low",]), method = "euclidean")
  
  
  
  
  perm_res <- adonis2(data_matrix ~ labels, method = "euclidean", permutations = 999)
  r2 <- perm_res$R2[1]       
  p_val <- perm_res$`Pr(>F)`[1] 
  
 
  return(list(dist=dist_val, r2=r2, p=p_val))
}


stats_train <- evaluate_separation(X_train_final, train_final_labels, "Training Set")
stats_test <- evaluate_separation(X_test_final, test_pred_labels, "Validation Set")


library(vegan)


dist_train <- dist(X_train_final, method = "euclidean")


mod_dispersion <- betadisper(dist_train, factor(train_final_labels))



res_dispersion <- permutest(mod_dispersion, pairwise = TRUE, permutations = 999)

print(res_dispersion)



plot(mod_dispersion, main = "Multivariate Dispersion Analysis")






cat("\n>>> Generating PCA Plots...\n")


plot_pca_custom <- function(data_matrix, labels, dataset_name) {
  
  pca_res <- prcomp(data_matrix, scale. = TRUE)
  
  
  pca_df <- data.frame(
    PC1 = pca_res$x[,1],
    PC2 = pca_res$x[,2],
    Group = factor(labels, levels = c("High", "Low"))
  )
  
  
  var_explained <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)
  
  
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
    geom_point(size = 2.5, alpha = 0.7) +
    
    stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95) + 
    scale_color_manual(values = c("High" = "#DC143C", "Low" = "#2E8B57")) +
    scale_fill_manual(values = c("High" = "#DC143C", "Low" = "#2E8B57")) +
    labs(
      title = paste0(dataset_name, " PCA Distribution"),
      x = paste0("PC1 (", var_explained[1], "%)"),
      y = paste0("PC2 (", var_explained[2], "%)")
    ) +
    theme_minimal() +
    theme(legend.position = "top")
  return(p)
}


p_pca_train <- plot_pca_custom(X_train_final, train_final_labels, "Training Set")




p_pca_test <- plot_pca_custom(X_test_final, test_pred_labels, "Validation Set")
library(ggplot2)


library(ggplot2)
library(gridExtra)
pca_combined <- grid.arrange(p_pca_train, p_pca_test, ncol = 2)





library(ggplot2)
library(gridExtra)
library(vegan) 

cat("\n>>> Generating PCA Plots with Statistical Annotation...\n")


plot_pca_custom <- function(data_matrix, labels, dataset_name) {
  
  
  
  perm_res <- adonis2(data_matrix ~ labels, method = "euclidean", permutations = 999)
  p_val <- perm_res$`Pr(>F)`[1]
  r2_val <- perm_res$R2[1]
  
  
  p_text <- ifelse(p_val < 0.001, "P < 0.001", paste0("P = ", sprintf("%.3f", p_val)))
  r2_text <- paste0("R2 = ", sprintf("%.3f", r2_val))
  stats_subtitle <- paste0( p_text, " | ", r2_text)
  
  
  
  pca_res <- prcomp(data_matrix, scale. = TRUE)
  
  
  pca_df <- data.frame(
    PC1 = pca_res$x[,1],
    PC2 = pca_res$x[,2],
    Group = factor(labels, levels = c("High", "Low"))
  )
  
  
  var_explained <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)
  
  
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, fill = Group, shape = Group)) +
    
    geom_point(size = 3, alpha = 0.7) + 
    
    
    stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, show.legend = FALSE) + 
    
    
    scale_color_manual(values = c("High" = "#56B4E9", "Low" = "#E64B35")) +
    scale_fill_manual(values = c("High" = "#56B4E9", "Low" = "#E64B35")) +
    
    scale_shape_manual(values = c("High" = 17, "Low" = 16)) + 
    
    labs(
      title = paste0(dataset_name, " PCA Distribution"),
      subtitle = stats_subtitle,  
      x = paste0("PC1 (", var_explained[1], "%)"),
      y = paste0("PC2 (", var_explained[2], "%)")
    ) +
    
    
    theme_bw() + 
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 12, color = "black", hjust = 0.5, face = "italic"), 
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      axis.title = element_text(size = 12),
      panel.grid.minor = element_blank() 
    )
  
  guides(color = guide_legend(override.aes = list(shape = c(17, 16)))) 
  return(p)
}




p_pca_train <- plot_pca_custom(X_train_final, train_final_labels, "Training cohort")



p_pca_test <- plot_pca_custom(X_test_final, test_pred_labels, "Validation cohort")


pca_combined <- grid.arrange(p_pca_train, p_pca_test, ncol = 2)











library(pheatmap)

view_heatmap_sorted <- function(data_matrix, labels, cr_status, dataset_name) {
  
  
  
  
  
  sort_order <- order(factor(labels, levels = c("High", "Low")))
  
  
  data_sorted <- data_matrix[sort_order, ]
  labels_sorted <- labels[sort_order]
  cr_sorted <- cr_status[sort_order]
  
  
  annotation_col <- data.frame(
    Group = factor(labels_sorted, levels = c("High", "Low")),
    CR_Status = factor(ifelse(cr_sorted == 1, "Yes", "No"), levels = c("Yes", "No"))
  )
  
  
  current_names <- rownames(data_sorted)
  rownames(annotation_col) <- current_names
  
  
  ann_colors <- list(
    Group = c("High" = "#DC143C", "Low" = "#2E8B57"),
    CR_Status = c("Yes" = "#377EB8", "No" = "#D95F02") 
  )
  
  
  
  gap_pos <- sum(labels_sorted == "High")
  
  
  pheatmap(t(data_sorted), 
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           scale = "row",        
           show_colnames = FALSE, 
           
           
           cluster_cols = FALSE,  
           cluster_rows = TRUE,   
           gaps_col = gap_pos,    
           
           
           fontsize_row = 8,
           treeheight_row = 40,
           main = paste0(dataset_name, " Feature pheatmap"),
           filename = NA
  )
}



if(is.null(rownames(X_train_final))) rownames(X_train_final) <- paste0("Train_", 1:nrow(X_train_final))
view_heatmap_sorted(X_train_final, train_final_labels, data_train$CR, "Training Set: ")



if(is.null(rownames(X_test_final))) rownames(X_test_final) <- paste0("Test_", 1:nrow(X_test_final))
view_heatmap_sorted(X_test_final, test_pred_labels, data_test$CR, "Validation Set: ")



library(pheatmap)



library(pheatmap)


view_heatmap_preview <- function(data_matrix, labels, cr_status, dataset_name, 
                                 tree_width = 50) { 
  
  
  sort_order <- order(factor(labels, levels = c("High", "Low")))
  data_sorted <- data_matrix[sort_order, ]
  labels_sorted <- labels[sort_order]
  cr_sorted <- cr_status[sort_order]
  
  
  annotation_col <- data.frame(
    Group = factor(labels_sorted, levels = c("High", "Low")),
    CR_Status = factor(ifelse(cr_sorted == 1, "Yes", "No"), levels = c("Yes", "No"))
  )
  rownames(annotation_col) <- rownames(data_sorted)
  
  
  ann_colors <- list(
    Group = c("High" = "#DC143C", "Low" = "#2E8B57"),
    CR_Status = c("Yes" = "#377EB8", "No" = "#D95F02") 
  )
  
  
  gap_pos <- sum(labels_sorted == "High")
  
  
  pheatmap(t(data_sorted), 
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           scale = "row",        
           show_colnames = FALSE, 
           cluster_cols = FALSE, 
           cluster_rows = TRUE,   
           gaps_col = gap_pos,
           
           
           fontsize_row = 8,
           treeheight_row = tree_width, 
           main = paste0(dataset_name, " Preview"),
           filename = NA  
  )
}




view_heatmap_preview(X_train_final, train_final_labels, data_train$CR, "Training Set", tree_width = 50)

view_heatmap_preview(X_test_final, test_pred_labels, data_test$CR, "Validation Set", tree_width = 50)

save_final_heatmap <- function(data_matrix, labels, cr_status, dataset_name, save_path, tree_width = 50) {
  
  sort_order <- order(factor(labels, levels = c("High", "Low")))
  data_sorted <- data_matrix[sort_order, ]
  labels_sorted <- labels[sort_order]
  cr_sorted <- cr_status[sort_order]
  
  annotation_col <- data.frame(
    Group = factor(labels_sorted, levels = c("High", "Low")),
    CR_Status = factor(ifelse(cr_sorted == 1, "Yes", "No"), levels = c("Yes", "No"))
  )
  rownames(annotation_col) <- rownames(data_sorted)
  
  ann_colors <- list(
    Group = c("High" = "#DC143C", "Low" = "#2E8B57"),
    CR_Status = c("Yes" = "#377EB8", "No" = "#D95F02") 
  )
  gap_pos <- sum(labels_sorted == "High")
  
  pheatmap(t(data_sorted), 
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           scale = "row",        
           show_colnames = FALSE, 
           cluster_cols = FALSE, 
           cluster_rows = TRUE,   
           gaps_col = gap_pos,
           fontsize_row = 8,
           treeheight_row = tree_width, 
           main = paste0(dataset_name, " Feature Heatmap"),
           
           
           filename = save_path, 
           width = 9,    
           height = 6
  )
}


heatmap_dir <- file.path(output_dir, "Visualization")
if(!dir.exists(heatmap_dir)) dir.create(heatmap_dir)


save_final_heatmap(X_train_final, train_final_labels, data_train$CR, "Training Set", 
                   save_path = file.path(heatmap_dir, "Heatmap_Training_Final.pdf"), 
                   tree_width = 50)


save_final_heatmap(X_test_final, test_pred_labels, data_test$CR, "Validation Set", 
                   save_path = file.path(heatmap_dir, "Heatmap_Validation_Final.pdf"), 
                   tree_width = 50)




