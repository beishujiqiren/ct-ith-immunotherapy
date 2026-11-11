
project_root <- normalizePath(Sys.getenv("ITH_PROJECT_ROOT", unset = getwd()),
                              winslash = "/", mustWork = TRUE)
results_dir <- file.path(project_root, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


library(tidyverse)
library(caret)
library(ConsensusClusterPlus)
library(pheatmap)
library(ggpubr)
library(openxlsx)


work_dir <- project_root
output_dir <- file.path(results_dir, "habitat_analysis")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)



path_cohort1 <- file.path(work_dir, "data", "clinical", "ITH_cohort1.csv")
path_cohort2 <- file.path(work_dir, "data", "clinical", "ITH_cohort2.csv")
path_cohort3 <- file.path(work_dir, "data", "clinical", "ITH_cohort3.csv")

df_cohort1 <- read.csv(path_cohort1, fileEncoding = "GB18030")
df_cohort2 <- read.csv(path_cohort2, fileEncoding = "GB18030")
df_cohort3 <- read.csv(path_cohort3, fileEncoding = "GB18030")

for (i in 3:length(colnames(df_cohort1))) {
  df_cohort1[[i]]<-as.numeric(df_cohort1[[i]])
  df_cohort2[[i]]<-as.numeric(df_cohort2[[i]])
  df_cohort3[[i]]<-as.numeric(df_cohort3[[i]])
}


data_train <- df_cohort1
data_test <- rbind(df_cohort2,df_cohort3)



data_train$CR <- as.numeric(data_train$CR)
data_test$CR <- as.numeric(data_test$CR)


feature_cols <- grep("^iTED_", names(data_train), value = TRUE)


preproc <- preProcess(data_train[, feature_cols], method = c("center", "scale"))


X_train_full <- predict(preproc, data_train[, feature_cols])
X_test_full <- predict(preproc, data_test[, feature_cols])





mads <- apply(X_train_full, 2, mad)

cutoff_mad <- quantile(mads, 0.5) 
feats_step1 <- names(mads)[mads > cutoff_mad]


corr_mat <- cor(X_train_full[, feats_step1])
high_corr <- findCorrelation(corr_mat, cutoff = 0.8)
if(length(high_corr) > 0) {
  feats_step2 <- feats_step1[-high_corr]
} else {
  feats_step2 <- feats_step1
}



p_values <- apply(X_train_full[, feats_step2], 2, function(x) {
  
  fit <- glm(data_train$CR ~ x, family = binomial())
  summary(fit)$coefficients[2, 4] 
})


final_features <- names(p_values)[p_values < 0.2]
print(final_features)




X_train_final <- as.matrix(X_train_full[, final_features])
X_test_final <- as.matrix(X_test_full[, final_features])




cc_dir <- file.path(output_dir, "ConsensusCluster")
if(!dir.exists(cc_dir)) dir.create(cc_dir)






set.seed(1234)
results_cc <- ConsensusClusterPlus(
  d = t(X_train_final), 
  maxK = 4,
  reps = 1000,
  pItem = 0.8,
  pFeature = 1,
  title = cc_dir,
  clusterAlg = "pam",
  distance = "euclidean",
  seed = 1234,
  plot = "pdf"
)


k_best <- 2
train_cluster_labels <- results_cc[[k_best]]$consensusClass


cluster_means <- data.frame(X_train_final) %>%
  mutate(Cluster = train_cluster_labels) %>%
  group_by(Cluster) %>%
  summarise(across(everything(), mean))


cluster_means$Overall_ITH <- rowMeans(cluster_means[, -1])
print(cluster_means[, c("Cluster", "Overall_ITH")])


high_cluster_id <- cluster_means$Cluster[which.max(cluster_means$Overall_ITH)]
low_cluster_id <- cluster_means$Cluster[which.min(cluster_means$Overall_ITH)]

get_ith_label <- function(ids) {
  ifelse(ids == high_cluster_id, "High","Low")
}

train_final_labels <- get_ith_label(train_cluster_labels)



centroid_high <- colMeans(X_train_final[train_cluster_labels == high_cluster_id, ])
centroid_low <- colMeans(X_train_final[train_cluster_labels == low_cluster_id, ])



calc_dist <- function(x, center) { sqrt(sum((x - center)^2)) }

test_pred_labels <- apply(X_test_final, 1, function(row) {
  d_high <- calc_dist(row, centroid_high)
  d_low <- calc_dist(row, centroid_low)
  if (d_high < d_low) return("High") else return("Low")
})



plot_data <- rbind(
  data.frame(Set = "Training Set", ITH = train_final_labels, CR = data_train$CR),
  data.frame(Set = "Validation Set", ITH = test_pred_labels, CR = data_test$CR)
)
plot_data$ITH <- factor(plot_data$ITH, levels = c("High", "Low"))


stats_summary <- plot_data %>%
  group_by(Set, ITH) %>%
  summarise(
    Count = n(),
    CR_Rate = mean(CR),
    CR_Num = sum(CR),
    .groups = 'drop'
  )

print(as.data.frame(stats_summary))


calc_p <- function(d) {
  t <- table(d$ITH, d$CR)
  if(min(t)<5) fisher.test(t)$p.value else chisq.test(t)$p.value
}

p_train <- calc_p(subset(plot_data, Set == "Training Set"))
p_test <- calc_p(subset(plot_data, Set == "Validation Set"))






library(ggstatsplot) 
library(patchwork)   
library(scales)      


plot_summary_ci <- plot_data %>%
  group_by(Set, ITH) %>%
  summarise(
    n = n(),
    cr_count = sum(CR),
    .groups = 'drop'
  ) %>%
  rowwise() %>%
  mutate(
    
    res = list(prop.test(cr_count, n, conf.level = 0.95)),
    Rate = res$estimate,
    Lower = res$conf.int[1],
    Upper = res$conf.int[2]
  ) %>%
  ungroup()


plot_data_factor <- plot_data %>%
  mutate(CR_Status = factor(ifelse(CR == 1, "CR", "Non-CR"), levels = c("Non-CR", "CR")))





p_main <- ggplot(plot_summary_ci, aes(x = ITH, y = Rate, fill = ITH)) +
  
  geom_bar(stat = "identity", width = 0.6, color = "black", size = 0.1) +
  
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.15, size = 0.7, color = "black") +
  
  facet_wrap(~Set, scales = "free_x") +
  
  scale_fill_manual(values = c("High" = "#A5BDB2", "Low" = "#F3BE92")) +
  
  scale_y_continuous(labels = scales::percent_format(), 
                     expand = c(0, 0),limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
  labs(
    y = "Complete Response Rate (%)",
    x = "ITH Group",
    fill = "ITH Level"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 16),
    legend.position = "none",
    axis.line = element_line(size = 0.8),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )




















ann_text <- data.frame(
  Set = c("Training Set", "Validation Set"),
  label = c(paste0("P ", ifelse(p_train < 0.001, "< 0.001", paste0("= ", round(p_train, 3)))),
            paste0("P ", ifelse(p_test < 0.001, "< 0.001", paste0("= ", round(p_test, 3)))))
)

p_final <- p_main + 
  geom_text(data = ann_text, aes(x = 1.5, y = 1.05, label = label), 
            inherit.aes = FALSE, size = 5, fontface = "italic")




print(p_final)

ggsave(file.path(output_dir, "Figure_ITH_CR_Rate.pdf"), p_final, width = 7, height = 6)
library(ggpubr)


get_sig_label <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}


bracket_data <- data.frame(
  Set = c("Training Set", "Validation Set"),
  xmin = c("High"),
  xmax = c("Low"),
  y.position = c(1.02, 1.02), 
  label = c(get_sig_label(p_train), get_sig_label(p_test))
)


p_refined <- p_main + 
  
  geom_bracket(
    data = bracket_data,
    aes(xmin = xmin, xmax = xmax, label = label, y.position = y.position),
    inherit.aes = FALSE,
    label.size = 5,           
    tip.length = 0.02,        
    vjust = 0.5               
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(), 
    limits = c(0, 1.15),      
    breaks = seq(0, 1, 0.2)
  ) +
  
  geom_text(
    data = ann_text, 
    aes(x = 1.5, y = 1.12, label = label), 
    inherit.aes = FALSE, 
    size = 4, 
    fontface = "italic"
  )


print(p_refined)

ggsave(file.path(output_dir, "CR rate compare.pdf"), p_refined, width = 7, height = 6)

