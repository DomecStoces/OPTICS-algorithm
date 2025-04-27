library(tidyverse)
library(dbscan)
library(ggplot2)

# Step 1: Simulate presence/absence dataset
set.seed(42)
species <- paste0("Sp", 1:9)
intervals <- paste0("Y", seq(1925, 2025, 15))
presence_matrix <- matrix(0, nrow = 9, ncol = 7,
                          dimnames = list(species, intervals))

# Omnipresent species (Sp1–Sp3): present in all intervals
presence_matrix[1:3, ] <- 1

# Appearing species (Sp4–Sp6): appear later
presence_matrix[4, 1:5] <- 1
presence_matrix[5, 1:4] <- 1
presence_matrix[6, 1:3] <- 1

# Disappearing species (Sp7–Sp9): vanish over time
presence_matrix[7, 5:7] <- 1
presence_matrix[8, 4:7] <- 1
presence_matrix[9, 3:7] <- 1

# Step 2: Calculate temporal metrics
df_metrics <- as.data.frame(presence_matrix) %>%
  rownames_to_column("species") %>%
  rowwise() %>%
  mutate(
    ubiquity = mean(c_across(starts_with("Y"))),
    first = which(c_across(starts_with("Y")) == 1)[1],
    last = tail(which(c_across(starts_with("Y")) == 1), 1),
    slope = coef(lm(as.numeric(c_across(starts_with("Y"))) ~ 
                      seq_along(c_across(starts_with("Y")))))[2],
    turnover = sum(abs(diff(as.numeric(c_across(starts_with("Y")))))) / 2
  ) %>%
  ungroup()

# Step 3: Cluster species using OPTICS
metrics_scaled <- df_metrics %>%
  select(ubiquity, first, last, slope, turnover) %>%
  scale()

opt <- optics(metrics_scaled, minPts = 2)
plot(opt, main = "Reachability Plot")

clusters <- extractDBSCAN(opt, eps_cl = 1.5)
df_metrics$cluster <- as.character(clusters$cluster)

# Step 4: Assign dynamic labels based on average slope
cluster_behavior <- df_metrics %>%
  filter(!is.na(cluster)) %>%
  group_by(cluster) %>%
  summarise(avg_slope = mean(slope), .groups = "drop") %>%
  mutate(cluster_name = case_when(
    avg_slope > 0.1 ~ "Appearing",
    avg_slope < -0.1 ~ "Disappearing",
    TRUE ~ "Omnipresent"
  ))

df_metrics <- df_metrics %>%
  left_join(cluster_behavior, by = "cluster") %>%
  mutate(cluster_name = if_else(is.na(cluster_name), "Unclustered", cluster_name))

# Step 5: PCA visualization
pca <- prcomp(metrics_scaled)
df_pca <- as.data.frame(pca$x[, 1:2]) %>%
  mutate(species = df_metrics$species, cluster = df_metrics$cluster_name)

ggplot(df_pca, aes(PC1, PC2, color = cluster, label = species)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5) +
  theme_minimal() +
  labs(title = "OPTICS Clustering of Simulated Species")

# Step 6: Plot presence timeline grouped by cluster
presence_df <- as.data.frame(presence_matrix) %>%
  rownames_to_column("species") %>%
  pivot_longer(-species, names_to = "interval", values_to = "presence") %>%
  left_join(df_metrics %>% select(species, cluster_name), by = "species") %>%
  mutate(
    species = fct_reorder(species, as.numeric(factor(cluster_name))),
    interval = factor(interval, levels = intervals)
  )

ggplot(presence_df, aes(x = interval, y = species, fill = factor(presence))) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("0" = "white", "1" = "#2C7BB6"), name = "Presence") +
  facet_wrap(~ cluster_name, scales = "free_y", ncol = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold")) +
  labs(title = "Species Presence Timeline by Cluster Type",
       x = "Time Interval",
       y = "Species")
