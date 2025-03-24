# DATA SET 1: Import and preprocess the lab-based data set
lab_data_original <- read.csv("./Dissertation/OG_PPQST.csv", header = TRUE)
rownames(lab_data_original) <- lab_data_original$PatPIN

# Compute summary statistics for initial inspection
lab_summary <- summary(lab_data_original)

# Select only numerical columns for analysis
lab_data <- lab_data_original[, sapply(lab_data_original, is.numeric)]

# Remove duplicate rows
lab_data <- lab_data[!duplicated(lab_data), ]

# Count incomplete cases of sensory profiles
incomplete_cases_lab <- sum(!complete.cases(lab_data))

# Replace missing values in each numeric column with the column mean
lab_data_filled <- as.data.frame(lapply(lab_data, function(column) {
  if (is.numeric(column)) {
    column[is.na(column)] <- mean(column, na.rm = TRUE)
  }
  return(column)
}))

# Normalise the dataset
lab_scaled <- as.data.frame(scale(lab_data_filled))

# Determine optimal clustering solution using different distance and linkage methods

# Define distance matrices: Euclidean and Manhattan
distance_methods <- list(
  euclidean = dist(lab_scaled, method = "euclidean"),
  manhattan = dist(lab_scaled, method = "manhattan")
)

# Specify linkage methods of interest
linkage_methods <- c("complete", "ward.D2", "average")

# Perform hierarchical clustering for each distance-linkage combination (6 total)
hclust_results_lab <- list()
for (dist_name in names(distance_methods)) {
  dist_matrix <- distance_methods[[dist_name]]
  hclust_results_lab[[dist_name]] <- lapply(linkage_methods, function(linkage) {
    hclust(dist_matrix, method = linkage)
  })
  names(hclust_results_lab[[dist_name]]) <- linkage_methods
}

# Plot dendrograms for visual evaluation of clustering solutions
library(ggplot2)
library(factoextra)
library(gridExtra)
library(stringr)

plot_list <- list()
for (dist_name in names(hclust_results_lab)) {
  for (linkage in linkage_methods) {
    p <- fviz_dend(hclust_results_lab[[dist_name]][[linkage]],
                   cex = 0.3,
                   show_labels = FALSE,
                   ggtheme = theme_classic()) +
      labs(title = paste(linkage, "linkage\n", dist_name, "distance"),
           x = "Patient Sensory Profiles", y = "Height")
    plot_list <- append(plot_list, list(p))
  }
}
grid.arrange(grobs = plot_list, nrow = 2, ncol = 3)

# Set the number of clusters based on observation
k <- 3

# Replot dendrograms with cluster assignments indicated by colors
plot_list <- list()
for (dist_name in names(hclust_results_lab)) {
  for (linkage in linkage_methods) {
    p <- fviz_dend(hclust_results_lab[[dist_name]][[linkage]],
                   k = k, 
                   k_colors = c("#FC4E07", "#E7B800", "#00AFBB"),
                   show_labels = FALSE,
                   cex = 0.5,
                   ggtheme = theme_classic()) +
      labs(title = paste(str_to_title(dist_name), "-", str_to_title(linkage), "Method"),
           x = "Patient Sensory Profiles", y = "Height") +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 24), 
        axis.title.x = element_text(face = "bold", margin = margin(t = 1), size = 17), 
        axis.title.y = element_text(face = "bold", size = 17, margin = margin(r = 10)),
        axis.text.y = element_text(face = "bold", size = 15),
        plot.margin = margin(t = 10, r = 10, b = 30, l = 10)
      )
    
    plot_list <- append(plot_list, list(p))
  }
}
grid.arrange(grobs = plot_list, nrow = 2, ncol = 3)

# Based on dendrogram inspection, remove the 'average' linkage method from further analysis
linkage_methods <- c("complete", "ward.D2")

# Update clustering results using only the selected linkage methods
for (dist_name in names(distance_methods)) {
  dist_matrix <- distance_methods[[dist_name]]
  hclust_results_lab[[dist_name]] <- lapply(linkage_methods, function(linkage) {
    hclust(dist_matrix, method = linkage)
  })
  names(hclust_results_lab[[dist_name]]) <- linkage_methods
}

# Generate new dendrograms using the refined clustering methods
plot_list_NEW <- list()
for (dist_name in names(hclust_results_lab)) {
  for (linkage in linkage_methods) {
    p <- fviz_dend(hclust_results_lab[[dist_name]][[linkage]],
                   k = k, 
                   k_colors = c("#FC4E07", "#E7B800", "#00AFBB"),
                   show_labels = FALSE,
                   cex = 0.5,
                   ggtheme = theme_classic()) +
      labs(title = paste(str_to_title(dist_name), "-", str_to_title(linkage), "Method"), 
           x = "Patient Sensory Profiles", 
           y = "Height") +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 30), 
        axis.title.x = element_text(face = "bold", margin = margin(t = 10), size = 24), 
        axis.title.y = element_text(face = "bold", size = 24, margin = margin(r = 10)),
        axis.text.y = element_text(face = "bold", size = 20),
        axis.text.x = element_text(face = "bold", size = 20),
        plot.margin = margin(t = 10, r = 10, b = 30, l = 10)
      )
    plot_list_NEW <- append(plot_list_NEW, list(p))
  }
}
grid.arrange(grobs = plot_list_NEW, nrow = 2, ncol = 2)

# Segment clusters for each distance-linkage combination with k=3 clusters
clusters_lab <- data.frame(Patient_ID = lab_data_original[, sapply(lab_data_original, is.character)])
for (dist_name in names(hclust_results_lab)) {
  for (linkage in linkage_methods) {
    hc_model <- hclust_results_lab[[dist_name]][[linkage]]
    hcutree_results <- cutree(hc_model, k = 3)
    col_name <- paste(dist_name, linkage, sep = " ")
    clusters_lab[[col_name]] <- hcutree_results
  }
}

# Evaluate consistency in cluster allocation using Fleiss' Kappa
library("irr")
lab_irr_complete <- kappam.fleiss(clusters_lab[, c(2, 4)], detail = TRUE)
lab_irr_Ward <- kappam.fleiss(clusters_lab[, c(3, 5)], detail = TRUE)

# Adjust cluster labels to improve comparability across methods
clusters_lab$`manhattan complete` <- ifelse(
  clusters_lab$`manhattan complete` == 2, 3, 
  ifelse(clusters_lab$`manhattan complete` == 3, 2, clusters_lab$`manhattan complete`))

clusters_lab$`manhattan ward.D2` <- ifelse(clusters_lab$`manhattan ward.D2` == 1, 3, 
                                           ifelse(clusters_lab$`manhattan ward.D2` == 2, 1, 
                                                  ifelse(clusters_lab$`manhattan ward.D2` == 3, 2, clusters_lab$`manhattan ward.D2`)))

# Recompute Fleiss' Kappa to assess improved consistency
lab_irr_complete <- kappam.fleiss(clusters_lab[, c(2, 4)], detail = TRUE)
lab_irr_Ward <- kappam.fleiss(clusters_lab[, c(3, 5)], detail = TRUE)

# Using Ward's method, employ silhouette scores to determine which distance metric better complements it
library(cluster)
lab_ward_euclidean <- silhouette(clusters_lab$`euclidean ward.D2`, dist = dist(lab_scaled, method = "euclidean"))
lab_ward_euclidean <- mean(lab_ward_euclidean[, 3])
lab_ward_manhattan <- silhouette(clusters_lab$`manhattan ward.D2`, dist = dist(lab_scaled, method = "manhattan"))
lab_ward_manhattan <- mean(lab_ward_manhattan[, 3])

# Plot the dendrogram of the optimal solution (Manhattan distance with Ward's linkage)
dend_lab <- as.dendrogram(hclust_results_lab[["manhattan"]][["ward.D2"]])
fviz_dend(dend_lab,
          k = k, 
          k_colors = c("#FC4E07", "#E7B800", "#00AFBB"), 
          show_labels = FALSE,
          cex = 0.7,  
          ggtheme = theme_classic()) +
  labs(title = "Lab-Data Set Clustering Dendrogram\nUsing the Manhattan-Ward Approach", 
       x = "Patient Sensory Profiles", 
       y = "Height") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 30, margin = margin(b = 5)), 
        axis.title.x = element_text(face = "bold", margin = margin(t = 1), size = 24), 
        axis.title.y = element_text(face = "bold", size = 24, margin = margin(r = 10)),
        axis.text.y = element_text(face = "bold", size = 20),
        axis.text.x = element_text(face = "bold", size = 20),
        plot.margin = margin(t = 10, r = 10, b = 30, l = 10)) +
  geom_hline(yintercept = 130, linetype = "dashed", color = "#333", size = 1)

# Calculate silhouette scores for cluster numbers k = 2 to 10 to find the optimal k
silhouette_summary <- list()
sil_avgs <- sapply(2:10, function(k) {
  clusters <- cutree(hclust_results_lab[["manhattan"]][["ward.D2"]], k = k)
  sil <- silhouette(clusters, distance_methods[[dist_name]])
  mean(sil[, 3])
})
best_k <- which.max(sil_avgs) + 1  
sil_per_k <- data.frame(
  k = c(1, 2:10),
  silhouette_width = c(0, sil_avgs)
)
key <- paste(dist_name, linkage, sep = "_")
silhouette_summary[[key]] <- list(
  best_k = best_k,
  mean_sil_width = mean(sil_avgs),
  sil_per_k = sil_per_k
)

# Create a flextable for the silhouette scores
library(flextable)
sil_table <- flextable(sil_per_k) %>%
  set_header_labels(k = "Number of Clusters (k)", 
                    silhouette_width = "Mean Silhouette Coefficient") %>%
  bold(part = "header") %>%
  bold(i = 2, j = 1) %>%
  bold(i = 3, j = 1) %>%
  bold(i = 2, j = 2) %>%
  bold(i = 3, j = 2) %>%
  autofit()
sil_table <- add_header_lines(sil_table, values = "Lab-Data Set")
sil_table

# Generate a line plot highlighting the optimal number of clusters
ggplot(sil_per_k, aes(x = k, y = silhouette_width)) +
  geom_line(color = "#00AFBB") +
  geom_point(color = "#00AFBB") +
  geom_vline(xintercept = best_k, linetype = "dashed", color = "black", size = 0.6) + 
  ggtitle("Silhouette Analysis of the Lab-Data Set") +
  ylim(0, 0.5) +
  xlab("Number of Clusters (k)") +
  ylab("Mean Silhouette Coefficient") +
  scale_x_continuous(breaks = seq(1, 10, 1)) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 24, margin = margin(b = 10)), 
        axis.title.x = element_text(face = "bold", margin = margin(t = 10), size = 20), 
        axis.title.y = element_text(face = "bold", size = 20, margin = margin(r = 10)),
        axis.text.y = element_text(face = "bold", size = 18),
        axis.text.x = element_text(face = "bold", size = 18),
        plot.margin = margin(t = 10, r = 10, b = 30, l = 10))

# Apply the optimal clustering labels to the original filled data
lab_data_filled$Clusters <- clusters_lab$`manhattan ward.D2`

# Segment data by cluster for further analysis of cluster characteristics
subset_lab_1 <- lab_data_filled[lab_data_filled$Clusters == 1, ]
subset_lab_2 <- lab_data_filled[lab_data_filled$Clusters == 2, ]
subset_lab_3 <- lab_data_filled[lab_data_filled$Clusters == 3, ]

# Prepare lab-based data for plotting parameter means per cluster
lab_pp <- data.frame(
  lab_data_filled[1:17],
  Clusters = lab_data_filled$Clusters
)
lab_qst <- data.frame(
  lab_data_filled[18:ncol(lab_data_filled)]
)
# Reorder and rename columns for clarity in downstream analysis
lab_pp <- lab_pp[, c(1, 2, 12, 8, 9, 11, 14, 15, 16, 17, 5, 4, 3, 10, 13, 7, 6, 18)]
colnames(lab_pp) <- c("Spontaneous burning sensation", 
                      "Spontaneous tingling sensation", 
                      "Spontaneous itching",  
                      "Spontaneous numbness", 
                      "Spontaneous pain in numb areas", 
                      "Spontaneous cold sensation", 
                      "Spontaneous squeezing sensation", 
                      "Deep pressure sensation", 
                      "Swelling sensations", 
                      "Sensation of tense muscles", 
                      "Sudden sharp pain with no known cause", 
                      "Sudden pain caused by moving, or certain positions", 
                      "Pain when brushed against lightly",  
                      "Pain by slight pressure", 
                      "Pain caused by a pointed object touching", 
                      "Pain by something cold", 
                      "Pain by something warm", 
                      "Cluster")

# Reshape data for plotting using melt
library(reshape2)
melted_lab_qst <- melt(lab_qst, id.vars = "Clusters", variable.name = "Parameter", value.name = "Outcome")
melted_lab_pp <- melt(lab_pp, id.vars = "Cluster", variable.name = "Parameter", value.name = "Outcome")

# Plot lab-QST results with lines marking the confidence interval at ±1.96
cluster_colors <- c("1" = "#FC4E07", "3" = "#E7B800", "2" = "#00AFBB")
ggplot(melted_lab_qst, aes(x = Parameter, y = Outcome, group = Clusters, color = as.factor(Clusters))) +
  geom_hline(yintercept = 1.96, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = -1.96, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") + 
  stat_summary(fun.data = mean_se, geom = "line", size = 0.5, colour = "black") +  
  stat_summary(fun.data = mean_se, geom = "point", size = 6) + 
  scale_color_manual(values = cluster_colors) +  
  labs(x = "QST Parameter", y = "Mean Z-score", title = "Lab-QST Outcomes") +
  scale_y_continuous(breaks = seq(-5, 7, 1)) +
  theme_classic() +  
  theme(legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 20),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 30), 
        axis.title.x = element_text(face = "bold", size = 24, margin = margin(t = 15)),
        axis.text.x = element_text(face = "bold", size = 20, angle = 90, hjust = 1, vjust = 0.5),
        axis.title.y = element_text(face = "bold", size = 24, margin = margin(r = 15)),
        axis.text.y = element_text(face = "bold", size = 20),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

# Plot corresponding PainPREDICT outcomes
ggplot(melted_lab_pp, aes(x = Parameter, y = Outcome, group = Cluster, color = as.factor(Cluster))) +
  stat_summary(fun.data = mean_se, geom = "line", size = 0.5, colour = "black") +  
  stat_summary(fun.data = mean_se, geom = "point", size = 6) + 
  scale_color_manual(values = cluster_colors) +  
  labs(x = "PainPREDICT Parameter", y = "Mean NRS-11 Score", title = "PainPREDICT Outcomes") +
  theme_classic() +  
  theme(legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 20),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 30), 
        axis.title.x = element_text(face = "bold", size = 24, margin = margin(t = 15)),
        axis.text.x = element_text(face = "bold", size = 20, angle = 90, hjust = 1, vjust = 0.5),
        axis.title.y = element_text(face = "bold", size = 24, margin = margin(r = 15)),
        axis.text.y = element_text(face = "bold", size = 20),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

# DATA SET 2: Import and preprocess bedside data
bedside_data_original <- read.csv("./Dissertation/OG_PPBedside.csv", header = TRUE)
rownames(bedside_data_original) <- bedside_data_original$PatPIN

# Generate summary statistics for initial review
bedside_summary <- summary(bedside_data_original)

# Select numerical columns only for analysis
bedside_data <- bedside_data_original[, sapply(bedside_data_original, is.numeric)]

# Remove duplicate rows
bedside_data <- bedside_data[!duplicated(bedside_data), ]

# Count incomplete sensory profiles 
incomplete_cases_bed <- sum(!complete.cases(bedside_data))

# Replace missing and outlier values with column means
bedside_filled <- as.data.frame(lapply(bedside_data, function(column) {
  if (is.numeric(column)) {
    column[is.na(column)] <- mean(column, na.rm = TRUE)
  }
  return(column)
}))

# Normalise the dataset
bedside_scaled <- as.data.frame(scale(bedside_filled))

# Derive clustering solution using Manhattan distance and Ward linkage
dist_bedside <- dist(bedside_scaled, method = "manhattan")
hc_bedside <- hclust(dist_bedside, method = "ward.D2")
bedside_clusters <- cutree(hc_bedside, k = 3)
bedside_filled$Clusters <- bedside_clusters

# Combine clustering labels from bedside and lab datasets for agreement analysis
ratings <- cbind(bedside_filled$Clusters, lab_data_filled$Clusters)

# Compute Fleiss' Kappa to measure inter-method consistency
library("irr")
lab_bedside_irr <- kappam.fleiss(ratings, detail = TRUE)

# Current agreement is low (0.09); relabeling of clusters evaluated to improve consistency

# Define candidate relabeling schemes
relabel_schemes <- list(
  "1->2, 2->3, 3->1" = c(2, 3, 1),
  "1->2, 2->1, 3->3" = c(2, 1, 3),
  "1->3, 2->1, 3->2" = c(3, 1, 2),
  "1->3, 2->2, 3->1" = c(3, 2, 1),
  "1->1, 2->3, 3->2" = c(1, 3, 2)
)

# Function to relabel clusters based on a given scheme
relabel_clusters <- function(data, scheme) {
  data$Clusters <- ifelse(data$Clusters == 1, scheme[1],
                          ifelse(data$Clusters == 2, scheme[2],
                                 ifelse(data$Clusters == 3, scheme[3], data$Clusters)))
  return(data)
}

kappa_results <- list()

# Testing each relabeling scheme and record Fleiss' Kappa
for (name in names(relabel_schemes)) {
  scheme <- relabel_schemes[[name]]
  
  relabeled_data <- relabel_clusters(bedside_filled, scheme)
  
  print(paste("Relabeling Scheme:", name))
  print(table(relabeled_data$Clusters))
  
  ratings <- cbind(relabeled_data$Clusters, lab_data_filled$Clusters)
  
  kappa_value <- kappam.fleiss(ratings, detail = TRUE)$value
  
  kappa_results[[name]] <- kappa_value
}

# After testing, retain the original labels as they show the best consistency with lab-QST data
dist_name <- "manhattan"
linakge <- "ward.D2"
dist_bedside <- dist(bedside_scaled, method = dist_name)
hc_bedside <- hclust(dist_bedside, method = linakge)
bedside_clusters <- cutree(hc_bedside, k = 3)
bedside_filled$Clusters <- bedside_clusters
bedside_data_original$Clusters <- bedside_clusters

# Silhouette analysis on the bedside-data set
silhouette_summary_bed <- list()
sil_avgs <- sapply(2:10, function(k) {
  clusters <- cutree(hc_bedside, k = k)
  sil <- silhouette(clusters, dist_bedside)
  mean(sil[, 3])
})
best_k <- which.max(sil_avgs) + 1 

sil_per_k <- data.frame(
  k = c(1, 2:10),
  silhouette_width = c(0, sil_avgs)
)
key <- paste(dist_name, linakge, sep = "_")
silhouette_summary_bed[[key]] <- list(
  best_k = best_k,
  mean_sil_width = mean(sil_avgs),
  sil_per_k = sil_per_k
)

# Create a flextable for silhouette scores
library(flextable)
sil_table <- flextable(sil_per_k) %>%
  set_header_labels(k = "Number of Clusters (k)",
                    silhouette_width = "Mean Silhouette Coefficient") %>%
  bold(part = "header") %>%
  bold(i = 2, j = 1) %>%
  bold(i = 2, j = 2) %>%
  autofit()
sil_table <- add_header_lines(sil_table, values = "Bedside-Data Set")
sil_table

# Plot a line graph of silhouette coefficients to visually identify optimal k
ggplot(sil_per_k, aes(x = k, y = silhouette_width)) +
  geom_line(color = "#00AFBB") +
  geom_point(color = "#00AFBB") +
  geom_vline(xintercept = best_k, linetype = "dashed", color = "black", size = 0.6) +
  ggtitle("Silhouette Analysis of the Bedside-Data Set") +
  ylim(0, 0.5) +
  xlab("Number of Clusters (k)") +
  ylab("Mean Silhouette Coefficient") +
  scale_x_continuous(breaks = seq(1, 10, 1)) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 24, margin = margin(b = 10)),
        axis.title.x = element_text(face = "bold", margin = margin(t = 10), size = 20),
        axis.title.y = element_text(face = "bold", size = 20, margin = margin(r = 10)),
        axis.text.y = element_text(face = "bold", size = 18),
        axis.text.x = element_text(face = "bold", size = 18),
        plot.margin = margin(t = 10, r = 10, b = 30, l = 10))

# Segment bedside data by clusters for subsequent analyses
subset_bed_1 <- bedside_filled[bedside_filled$Clusters == 1, ]
subset_bed_2 <- bedside_filled[bedside_filled$Clusters == 2, ]
subset_bed_3 <- bedside_filled[bedside_filled$Clusters == 3, ]

# Plot the dendrogram for the bedside dataset
dend_bedside <- as.dendrogram(hc_bedside)
dend_bedside <- rev(dend_bedside)  # Reverse for consistency with lab-based dendrogram
fviz_dend(dend_bedside,
          k = k,
          k_colors = c("1" = "#FC4E07", "2" = "#00AFBB", "3" = "#E7B800"),
          show_labels = FALSE,
          cex = 0.7,
          ggtheme = theme_classic()
) +
  labs(title = "Bedside-Data Set Clustering Dendrogram\nUsing the Manhattan-Ward Approach",
       x = "Patient Sensory Profiles",
       y = "Height") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 30, margin = margin(b = 5)),
        axis.title.x = element_text(face = "bold", margin = margin(t = 1), size = 24),
        axis.title.y = element_text(face = "bold", size = 24, margin = margin(r = 10)),
        axis.text.y = element_text(face = "bold", size = 20),
        axis.text.x = element_text(face = "bold", size = 20),
        plot.margin = margin(t = 10, r = 10, b = 30, l = 10)) +
  geom_hline(yintercept = 260, linetype = "dashed", color = "#333", size = 1)

# Prepare bedside data for plotting outcomes (bedside-QST and PainPREDICT)
bedside_pp <- data.frame(
  bedside_filled[1:17],
  Clusters = bedside_filled$Clusters
)
bedside_qst <- data.frame(
  bedside_filled[18:ncol(bedside_filled)]
)

# Reorder and rename columns for downstream processing
bedside_pp <- bedside_pp[, c(1,2,12,8,9,11,14,15,16,17,5,4,3,10,13,7,6,18)]
colnames(bedside_pp) <- c("Spontaneous burning sensation",
                          "Spontaneous tingling sensation",
                          "Spontaneous itching",
                          "Spontaneous numbness",
                          "Spontaneous pain in numb areas",
                          "Spontaneous cold sensation",
                          "Spontaneous squeezing sensation",
                          "Deep pressure sensation",
                          "Swelling sensations",
                          "Sensation of tense muscles",
                          "Sudden sharp pain with no known cause",
                          "Sudden pain caused by moving, or certain positions",
                          "Pain when brushed against lightly",
                          "Pain by slight pressure",
                          "Pain caused by a pointed object touching",
                          "Pain by something cold",
                          "Pain by something warm",
                          "Cluster")

library(reshape2)
melted_bedside_qst <- melt(bedside_qst, id.vars = "Clusters", variable.name = "Parameter", value.name = "Outcome")
melted_bedside_pp <- melt(bedside_pp, id.vars = "Cluster", variable.name = "Parameter", value.name = "Outcome")

# Plot bedside PainPREDICT outcomes
ggplot(melted_bedside_pp, aes(x = Parameter, y = Outcome, group = Cluster, color = as.factor(Cluster))) +
  stat_summary(fun.data = mean_se, geom = "line", size = 0.5, colour = "black") +
  stat_summary(fun.data = mean_se, geom = "point", size = 6) +
  scale_color_manual(values = cluster_colors) +
  labs(x = "PainPREDICT Parameter", y = "Mean NRS-11 Score", title = "PainPREDICT Outcomes") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 20),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 30),
        axis.title.x = element_text(face = "bold", size = 24, margin = margin(t = 15)),
        axis.text.x = element_text(face = "bold", size = 20, angle = 90, hjust = 1, vjust = 0.5),
        axis.title.y = element_text(face = "bold", size = 24, margin = margin(r = 15)),
        axis.text.y = element_text(face = "bold", size = 20),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

# Create descriptive tables for interval-scales bedside-QST parameters
library(flextable)

# Define the desired order for variables
original_order <- c(
  "Metal22_cold", "Metal22_felt", "Metal8_cold", "Metal8_felt",
  "Metal37_warm", "Metal37_felt", "Metal45_warm", "Metal45_felt",
  "Metal22_pain", "Metal8_pain", "Metal37_painintensity", "Metal45_painintensity",
  "Qtip_intensity", "Neurotip_intensity", "Neurotip_interval",
  "Qswipe_intensity", "Pressure_intensity", "Pressure_threshold", "Vibration", "Clusters"
)

original_order <- original_order[original_order %in% names(bedside_data_original)]

# Identify detection/presentation parameters (ending in "_felt")
felt_parameters <- grep("_felt$", original_order, value = TRUE)

# Identify numeric parameters excluding 'felt' parameters
numeric_parameters <- setdiff(original_order[sapply(bedside_data_original[ , original_order], is.numeric)],
                              felt_parameters)

clusters <- unique(bedside_data_original$Clusters)
results_list <- data.frame()

# Compute descriptive statistics for each parameter by cluster
for (param in original_order) {
  for (cluster in clusters) {
    cluster_data <- bedside_data_original[bedside_data_original$Clusters == cluster, param, drop = FALSE]
    if (param %in% felt_parameters) {
      results_list <- rbind(results_list,
                            data.frame(Parameter = param,
                                       Cluster = cluster,
                                       Min = "-",
                                       Max = "-",
                                       `Mean ± SD` = "-"))
    } else if (param %in% numeric_parameters) {
      min_val <- min(cluster_data[[param]], na.rm = TRUE)
      max_val <- max(cluster_data[[param]], na.rm = TRUE)
      mean_val <- round(mean(cluster_data[[param]], na.rm = TRUE), 2)
      sd_val <- round(sd(cluster_data[[param]], na.rm = TRUE), 2)
      mean_sd <- paste0(mean_val, " ± ", sd_val)
      results_list <- rbind(results_list,
                            data.frame(Parameter = param,
                                       Cluster = cluster,
                                       Min = min_val,
                                       Max = max_val,
                                       `Mean ± SD` = mean_sd))
    }
  }
}

final_table <- as.data.frame(results_list)
final_table$Parameter <- factor(final_table$Parameter, levels = original_order)
final_table <- final_table[order(final_table$Parameter), ]

library(magrittr)

ft <- flextable(final_table) %>%
  merge_v(j = "Parameter") %>%
  theme_vanilla() %>%
  autofit() %>%
  width(width = 1) %>%
  add_header_lines(values = "Table 1: Descriptive Analysis of Interval-Scaled Bedside-QST Parameters")
save_as_docx(ft, path = "table.docx", pr_section = sect_properties)

ft

# Calculate detection rates for specific parameters
detection_rates <- list()
calc_rate <- function(df, column, value) {
  total <- sum(!is.na(df[[column]]))
  count <- sum(df[[column]] == value, na.rm = TRUE)
  if (total > 0) {
    return(round((count / total) * 100, 2))
  } else {
    return("-")
  }
}

for (cluster in clusters) {
  cluster_data <- subset(bedside_data_original, Clusters == cluster)
  detection_rates <- rbind(detection_rates, data.frame(Parameter = "Metal22_cold", Cluster = cluster, Result = calc_rate(cluster_data, "Metal22_felt", 2)))
  detection_rates <- rbind(detection_rates, data.frame(Parameter = "Metal22_felt", Cluster = cluster, Result = calc_rate(cluster_data, "Metal22_felt", 3)))
  detection_rates <- rbind(detection_rates, data.frame(Parameter = "Metal8_cold", Cluster = cluster, Result = calc_rate(cluster_data, "Metal8_felt", 2)))
  detection_rates <- rbind(detection_rates, data.frame(Parameter = "Metal8_felt", Cluster = cluster, Result = calc_rate(cluster_data, "Metal8_felt", 3)))
  detection_rates <- rbind(detection_rates, data.frame(Parameter = "Metal37_warm", Cluster = cluster, Result = calc_rate(cluster_data, "Metal37_felt", 3)))
  detection_rates <- rbind(detection_rates, data.frame(Parameter = "Metal37_felt", Cluster = cluster, Result = calc_rate(cluster_data, "Metal37_felt", 2)))
  detection_rates <- rbind(detection_rates, data.frame(Parameter = "Metal45_warm", Cluster = cluster, Result = calc_rate(cluster_data, "Metal45_felt", 3)))
  detection_rates <- rbind(detection_rates, data.frame(Parameter = "Metal45_felt", Cluster = cluster, Result = calc_rate(cluster_data, "Metal45_felt", 2)))
  detection_rates <- rbind(detection_rates, data.frame(Parameter = "Metal22_pain", Cluster = cluster, Result = calc_rate(cluster_data, "Metal22_pain", 1)))
  detection_rates <- rbind(detection_rates, data.frame(Parameter = "Metal8_pain", Cluster = cluster, Result = calc_rate(cluster_data, "Metal8_pain", 1)))
  detection_rates <- rbind(detection_rates, data.frame(Parameter = "Metal37_painintensity", Cluster = cluster, Result = calc_rate(cluster_data, "Metal37_painintensity", 1)))
  detection_rates <- rbind(detection_rates, data.frame(Parameter = "Metal45_painintensity", Cluster = cluster, Result = calc_rate(cluster_data, "Metal45_painintensity", 1)))
  detection_rates <- rbind(detection_rates, data.frame(Parameter = "Qtip_intensity", Cluster = cluster, Result = "-"))
  detection_rates <- rbind(detection_rates, data.frame(Parameter = "Neurotip_intensity", Cluster = cluster, Result = calc_rate(cluster_data, "Neurotip_painful", 1)))
  detection_rates <- rbind(detection_rates, data.frame(Parameter = "Neurotip_interval", Cluster = cluster, Result = "-"))
  detection_rates <- rbind(detection_rates, data.frame(Parameter = "Qswipe_intensity", Cluster = cluster, Result = calc_rate(cluster_data, "Qswipe_painful", 1)))
  detection_rates <- rbind(detection_rates, data.frame(Parameter = "Pressure_intensity", Cluster = cluster, Result = calc_rate(cluster_data, "Pressure_painful", 1)))
  detection_rates <- rbind(detection_rates, data.frame(Parameter = "Pressure_threshold", Cluster = cluster, Result = "-"))
  detection_rates <- rbind(detection_rates, data.frame(Parameter = "Vibration", Cluster = cluster, Result = "-"))
}

final_table <- merge(final_table, detection_rates, by = c("Parameter", "Cluster"), all.x = TRUE)
colnames(final_table) <- c("Parameter", "Cluster", "Min", "Max", "Mean ± SD", "Detection rate (%)")
final_table[is.na(final_table)] <- "-"

row_name_mapping <- c(
  "Metal22_cold" = "Metal 22˚C Perception Intensity",
  "Metal22_felt" = "Metal 22˚C Paradoxical Heat Sensation",
  "Metal8_cold" = "Metal 8˚C Perception Intensity",
  "Metal8_felt" = "Metal 8˚C Paradoxical Heat Sensation",
  "Metal37_warm" = "Metal 37˚C Perception Intensity",
  "Metal37_felt" = "Metal 37˚C Paradoxical Heat Sensation",
  "Metal45_warm" = "Metal 45˚C Perception Intensity",
  "Metal45_felt" = "Metal 45˚C Paradoxical Heat Sensation",
  "Metal22_pain" = "Metal 22˚C Pain Intensity",
  "Metal8_pain" = "Metal 8˚C Pain Intensity",
  "Metal37_painintensity" = "Metal 37˚C Pain Intensity",
  "Metal45_painintensity" = "Metal 45˚C Pain Intensity",
  "Qtip_intensity" = "Q-tip Perception Intensity",
  "Neurotip_intensity" = "Neurotip Pain Intensity",
  "Neurotip_interval" = "Neurotip WUR Ratio Pain Intensity",
  "Qswipe_intensity" = "DMA Allodynia Pain Intensity",
  "Pressure_intensity" = "Pressure Algometer at 4-mL Pain Intensity",
  "Pressure_threshold" = "Pressure Algometer Pain Pressure Threshold",
  "Vibration" = "Vibration Detection Threshold"
)

final_table$Parameter <- ifelse(final_table$Parameter %in% names(row_name_mapping),
                                row_name_mapping[final_table$Parameter],
                                final_table$Parameter)

# Define section properties for Word export and create the coresponding flextable
library(officer)
sect_properties <- prop_section(
  page_size = page_size(orient = "portrait", width = 8.3, height = 11.7),
  type = "continuous",
  page_margins = page_mar()
)

ft <- flextable(final_table) %>%
  merge_v(j = "Parameter") %>%
  theme_vanilla() %>%
  autofit() %>%
  width(width = 1) %>%
  add_header_lines(values = "Table 1: Descriptive Analysis of Interval-Scaled Bedside-QST Parameters")
save_as_docx(ft, path = "table.docx", pr_section = sect_properties)

ft

# Plot stacked bar plots for Q-tip and Neurotip detection comparison across clusters
detection_qtip <- list()
detection_neuropen <- list()
for (cluster in clusters) {
  cluster_data <- subset(bedside_data_original, Clusters == cluster)
  detection_qtip <- rbind(detection_qtip, data.frame(Parameter = "Equal intensity\nacross arms", Cluster = cluster, Result = calc_rate(cluster_data, "Qtip_comparison", 1)))
  detection_qtip <- rbind(detection_qtip, data.frame(Parameter = "Intensity greater in\nthe control arm", Cluster = cluster, Result = calc_rate(cluster_data, "Qtip_comparison", 2)))
  detection_qtip <- rbind(detection_qtip, data.frame(Parameter = "Intensity greater in\nthe affected arm", Cluster = cluster, Result = calc_rate(cluster_data, "Qtip_comparison", 3)))
  detection_neuropen <- rbind(detection_neuropen, data.frame(Parameter = "Lack of sensation", Cluster = cluster, Result = calc_rate(cluster_data, "Neuropen_comparison", 1)))
  detection_neuropen <- rbind(detection_neuropen, data.frame(Parameter = "Blunt touch", Cluster = cluster, Result = calc_rate(cluster_data, "Neuropen_comparison", 2)))
  detection_neuropen <- rbind(detection_neuropen, data.frame(Parameter = "Pinprick", Cluster = cluster, Result = calc_rate(cluster_data, "Neuropen_comparison", 3)))
}

detection_qtip_matrix <- xtabs(Result ~ Cluster + Parameter, data = detection_qtip)
par(mar = c(5, 7, 4, 13))
detection_qtip_plot <- barplot(t(detection_qtip_matrix),
        horiz = TRUE,
        col = gray.colors(3, start = 0.3, end = 0.9),
        xlab = "Detection Rate (%)",
        ylab = "Cluster",
        main = "Variation in Q-tip Perception Across Clusters",
        cex.main = 1.4,
        font.main = 2,
        cex.lab = 1.2,
        font.lab = 2,
        cex.axis = 1,
        font.axis = 2,
        las = 2)
legend("topright", legend = colnames(detection_qtip_matrix), 
       fill = gray.colors(3, start = 0.3, end = 0.9),
       inset = c(-0.5, 0), xpd = TRUE, y.intersp = 2)

detection_neuropen_matrix <- xtabs(Result ~ Cluster + Parameter, data = detection_neuropen)
par(mar = c(5, 7, 4, 13))
detection_neuropen_plot <- barplot(t(detection_neuropen_matrix),
        horiz = TRUE,
        col = gray.colors(3, start = 0.3, end = 0.9),
        xlab = "Detection Rate (%)",
        ylab = "Cluster",
        main = "Variation in Neuropen Perception Across Clusters",
        cex.main = 1.4,
        font.main = 2,
        cex.lab = 1.2,
        font.lab = 2,
        cex.axis = 1,
        font.axis = 2,
        las = 2)
legend("topright", legend = colnames(detection_neuropen_matrix), 
       fill = gray.colors(3, start = 0.3, end = 0.9),

       inset = c(-0.5, 0), xpd = TRUE)
