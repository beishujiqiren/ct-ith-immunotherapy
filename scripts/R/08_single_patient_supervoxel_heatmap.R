
project_root <- normalizePath(Sys.getenv("ITH_PROJECT_ROOT", unset = getwd()),
                              winslash = "/", mustWork = TRUE)
results_dir <- file.path(project_root, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)



library(oro.nifti)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggnewscale)
library(scales) 




base_dir <- project_root

global_csv_path <- file.path(base_dir, "results/supervoxel_features/RF_AllPatients_merged.csv") 


target_patient_id <-  "REPLACE_WITH_DEIDENTIFIED_ID"


sv_mask_path  <- file.path(base_dir, "results/supervoxels", paste0(target_patient_id, ".SVmask.nii.gz"))
ct_image_path <- file.path(base_dir, "results/preprocessed_ct", paste0(target_patient_id, "_preprocessed.nii.gz"))
sv_csv_path   <- file.path(base_dir, "results/supervoxel_features", paste0(target_patient_id, "_RF_AllSV.csv"))




df_global <- read.csv(global_csv_path)





plot_sv_heatmap_fixed <- function(ct_path, mask_path, local_csv_path, feature_name, 
                                  global_limits, slice_idx = NULL) {
  
  
  ct <- readNIfTI(ct_path, reorient = FALSE)
  mask <- readNIfTI(mask_path, reorient = FALSE)
  df_sv_local <- read.csv(local_csv_path)
  
  
  
  
  df_sv_local$SV_ID_Match <- as.integer( df_sv_local$SV_Label)
  
  
  
  clean_col_name <- gsub("iTED_", "", feature_name)
  
  
  if (!clean_col_name %in% colnames(df_sv_local)) {
    
    clean_col_name_short <- gsub("original_", "", clean_col_name)
    if (clean_col_name_short %in% colnames(df_sv_local)) {
      clean_col_name <- clean_col_name_short
    } else {
      warning(paste("⚠️", feature_name))
      return(NULL)
    }
  }
  
  
  sv_data <- df_sv_local %>% 
    select(SV_ID_Match, all_of(clean_col_name)) %>% 
    rename(SV_ID = SV_ID_Match, Value = all_of(clean_col_name))
  
  
  if(is.null(slice_idx)) {
    slice_areas <- apply(mask, 3, function(x) sum(x > 0))
    slice_idx <- which.max(slice_areas)
  }
  
  
  ct_slice <- ct[,,slice_idx]
  mask_slice <- mask[,,slice_idx]
  
  ct_df <- reshape2::melt(ct_slice) %>% rename(x = Var1, y = Var2, Intensity = value)
  
  mask_df <- reshape2::melt(mask_slice) %>% 
    rename(x = Var1, y = Var2, SV_ID = value) %>% 
    filter(SV_ID > 0) %>% 
    left_join(sv_data, by = "SV_ID") 
  
  
  p <- ggplot() +
    
    geom_raster(data = ct_df, aes(x = x, y = y, fill = Intensity)) +
    scale_fill_gradient(low = "black", high = "white", guide = "none") +
    
    
    ggnewscale::new_scale_fill() +
    geom_raster(data = mask_df, aes(x = x, y = y, fill = Value), alpha = 0.85) +
    
    
    scale_fill_viridis_c(
      option = "turbo",        
      name = "Value",
      limits = global_limits,  
      oob = scales::squish,    
      na.value = "transparent"
    ) +
    
    coord_fixed() +
    theme_void() +
    labs(title = clean_col_name) +
    theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
          legend.position = "right",
          legend.key.height = unit(1, "cm"))
  
  return(p)
}



plot_list <- list()

for(feat in final_features) {
  
  
  global_col <- gsub("iTED_", "", feat)
  if (!global_col %in% colnames(df_global)) {
    global_col_short <- gsub("original_", "", global_col)
    if(global_col_short %in% colnames(df_global)) {
      global_col <- global_col_short
    } else { next }
  }
  
  
  global_vals <- df_global[[global_col]]
  
  limit_low <- 0
  
  
  q985 <- quantile(global_vals, probs = 0.99, na.rm = TRUE) 
  q995 <- quantile(global_vals, probs = 0.995, na.rm = TRUE) 
  val_max <- max(global_vals, na.rm = TRUE)                  
  
  
  
  
  tail_ratio <- val_max / (q985 + 1e-6) 
  
  if (tail_ratio > 1.5) {
    
    limit_high <- q985
    strategy_name <- "Aggressive Clip (98.5%)"
  } else {
    
    
    limit_high <- q995 
    strategy_name <- "Soft Clip (99.5%)"
  }
  
  
  if (limit_high <= limit_low) { limit_high <- val_max }
  
  final_limits <- c(limit_low, limit_high)
  

  
  
  p <- plot_sv_heatmap_fixed(ct_image_path, sv_mask_path, sv_csv_path, feat, 
                             global_limits = final_limits)
  
  if(!is.null(p)) {
    plot_list[[feat]] <- p
  }
}


final_plot <- wrap_plots(plot_list, ncol = 3) + 
  plot_annotation(
    title = paste("High ITH patient"),
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

print(final_plot)

ggsave(
  filename = "heatmap of single patient High.pdf",  
  plot = final_plot,
  device = "pdf",  
  width = 15,      
  height = 12,     
  dpi = 300        
)

