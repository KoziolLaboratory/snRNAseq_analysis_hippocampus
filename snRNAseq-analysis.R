library(Seurat)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(dplyr)
library(tidyr)

setwd("D:/Koziol_lab/AD_Project")
getwd()
# Helper function to load one sample
load_sample <- function(barcode_path, features_path, matrix_path, sample_name) {
  counts <- ReadMtx(mtx = matrix_path,
                    features = features_path,
                    cells = barcode_path)
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name)
  seurat_obj$group <- sample_name
  return(seurat_obj)
}

# Load samples
A1 <- load_sample("D:/Koziol_lab/AD_Project/GSE175814/GSM5348374_A1_barcodes.tsv.gz",
                  "D:/Koziol_lab/AD_Project/GSE175814/GSM5348374_A1_features.tsv.gz",
                  "D:/Koziol_lab/AD_Project/GSE175814/GSM5348374_A1_matrix.mtx.gz", "AD")

A2 <- load_sample("D:/Koziol_lab/AD_Project/GSE175814/GSM5348375_A2_barcodes.tsv.gz",
                  "D:/Koziol_lab/AD_Project/GSE175814/GSM5348375_A2_features.tsv.gz",
                  "D:/Koziol_lab/AD_Project/GSE175814/GSM5348375_A2_matrix.mtx.gz", "Control")

A3 <- load_sample("D:/Koziol_lab/AD_Project/GSE175814/GSM5348376_A3_barcodes.tsv.gz",
                  "D:/Koziol_lab/AD_Project/GSE175814/GSM5348376_A3_features.tsv.gz",
                  "D:/Koziol_lab/AD_Project/GSE175814/GSM5348376_A3_matrix.mtx.gz", "AD")

A4 <- load_sample("D:/Koziol_lab/AD_Project/GSE175814/GSM5348377_A4_barcodes.tsv.gz",
                  "D:/Koziol_lab/AD_Project/GSE175814/GSM5348377_A4_features.tsv.gz",
                  "D:/Koziol_lab/AD_Project/GSE175814/GSM5348377_A4_matrix.mtx.gz", "Control")

# Merge and add cell identifiers
combined <- merge(A1, y = list(A2, A3, A4), add.cell.ids = c("A1", "A2", "A3", "A4"))
combined$condition <- ifelse(grepl("AD", combined$group), "AD", "Control")

Idents(combined) <- combined$condition


# Seurat processing
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- RunUMAP(combined, dims = 1:10)
combined <- FindNeighbors(combined, dims = 1:10)
combined <- FindClusters(combined, resolution = 0.5)

###################################################
# UMAP plot for AD and Control Feature
###################################################
# Get UMAP
umap <- Embeddings(combined, "umap") |> as.data.frame()
colnames(umap)[1:2] <- c("UMAP_1","UMAP_2")
umap$condition <- factor(combined$condition[rownames(umap)], levels = c("Control","AD"))

# Split so draw order is Control (below) then AD (above)
df_lo <- subset(umap, condition == "Control")
df_hi <- subset(umap, condition == "AD")

# Right-side label positions
xr <- range(umap$UMAP_1); yr <- range(umap$UMAP_2)
dx <- diff(xr); dy <- diff(yr)
x_lab <- xr[2] + 0.08*dx
y_mid <- mean(yr)

# 4) Plot
ad_ctl <- ggplot() +
  geom_point(data = df_lo, aes(UMAP_1, UMAP_2, color = condition),
             size = 0.5, alpha = 0.75, show.legend = TRUE) +
  geom_point(data = df_hi, aes(UMAP_1, UMAP_2, color = condition),
             size = 0.5, alpha = 0.9, show.legend = FALSE) +
  scale_color_manual(values = c(Control = "azure4", AD = "brown3"),
                     labels = c("Control","AD"), drop = TRUE) +
  coord_cartesian(xlim = c(xr[1], xr[2] + 0.12*dx), clip = "off") +
  coord_equal() +
  theme_classic(base_size = 12) +
  theme(
    plot.background  = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.background= element_rect(fill = "transparent", color = NA),
    legend.key       = element_rect(fill = "transparent", color = NA),
    plot.margin      = margin(5.5, 50, 5.5, 5.5)
  ) +
  guides(color = guide_legend(title = NULL, override.aes = list(size = 3, alpha = 1))) +
  labs(x = "UMAP_1", y = "UMAP_2")

# Save / show
dir.create("Figures2", showWarnings = FALSE, recursive = TRUE)
ggsave("Figures2/UMAP_condition_right_labels.png", ad_ctl,
       width = 5.8, height = 4, dpi = 600, bg = "transparent")
ad_ctl
#######################################################
#UMAP for ALKBH3 and PINK1
########################################################
alk_features <- FeaturePlot(combined,
                            features = c("ALKBH3"),
                            cols = c("lightgrey", "red"),
                            reduction = "umap",
                            pt.size = 0.5,
                            split.by = "condition",  # <- this adds condition panels
                            order = TRUE)
print(alk_features)
ggsave("Figures2/All_features_ALKBH3.png", alk_features, width = 8, height = 4, dpi = 600)

pink_features <- FeaturePlot(combined,
                             features = c("ALKBH3"),
                             cols = c("lightgrey", "red"),
                             reduction = "umap",
                             pt.size = 0.5,
                             split.by = "condition",  # <- this adds condition panels
                             order = TRUE)
ggsave("Figures2/All_features_PINK1.tif", pink_features, width = 8, height = 4, dpi = 600)
##################################################################
# Set High and low filter for ALKBH3
##################################################################
exp_data <- FetchData(combined, vars = c("ALKBH3","condition"))

# Thresholds (tweak as needed)
hi_thr  <- 1.0   # High in AD if ALKBH3 > 1
low_thr <- 1.0   # Low in Control if ALKBH3 < 1

grp <- rep("Other", nrow(exp_data))
grp[exp_data$condition == "AD"      & exp_data$ALKBH3 >=  hi_thr] <- "High_in_AD"
grp[exp_data$condition == "Control" & exp_data$ALKBH3 <  low_thr] <- "Low_in_Control"
combined$ALKBH3_expr_group <- factor(grp, levels = c("Other","Low_in_Control","High_in_AD"))

# --- Highlight overlay (blue under, red on top) ---
emb <- Embeddings(combined, "umap") %>% as.data.frame()
colnames(emb)[1:2] <- c("UMAP_1","UMAP_2")
emb$group <- combined$ALKBH3_expr_group[rownames(emb)]

df_hilo <- emb %>%
  filter(group %in% c("Low_in_Control","High_in_AD")) %>%
  mutate(group = factor(group, levels = c("Low_in_Control","High_in_AD")))
df_lo <- dplyr::filter(df_hilo, group == "Low_in_Control")
df_hi <- dplyr::filter(df_hilo, group == "High_in_AD")

# ranges for right-side labels
xr <- range(df_hilo$UMAP_1); yr <- range(df_hilo$UMAP_2)
dx <- diff(xr); dy <- diff(yr)
x_lab <- xr[2] + 0.08*dx
y_mid <- mean(yr)

p <- ggplot() +
  geom_point(data = df_lo, aes(UMAP_1, UMAP_2, color = group),
             size = 0.5, alpha = 0.75, show.legend = TRUE) +
  geom_point(data = df_hi, aes(UMAP_1, UMAP_2, color = group),
             size = 0.5, alpha = 0.9, show.legend = FALSE) +
  scale_color_manual(values = c("Low_in_Control" = "blue", "High_in_AD" = "brown3"),
                     labels = c("Low in Control", "High in AD"), drop = TRUE) +
  coord_cartesian(xlim = c(xr[1], xr[2] + 0.12*dx), clip = "off") +
  coord_equal() +
  theme_classic(base_size = 12) +
  theme(
    plot.background  = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.background= element_rect(fill = "transparent", color = NA),
    legend.key       = element_rect(fill = "transparent", color = NA),
    plot.margin      = margin(5.5, 50, 5.5, 5.5)
  ) +
  guides(color = guide_legend(title = NULL, override.aes = list(size = 3, alpha = 1))) +
  labs(x = "UMAP_1", y = "UMAP_2")

dir.create("Figures2", showWarnings = FALSE, recursive = TRUE)
ggsave("Figures2/alk_high_low_only.tif", p, width = 5.8, height = 4, dpi = 600)
print(p)
##################################################################
# Set High and low filter for PINK1
##################################################################
exp_data <- FetchData(combined, vars = c("PINK1","condition"))

# Thresholds (tweak as needed)
hi_thr  <- 1.0   # High in Control if PINK1 > 1
low_thr <- 1.0   # Low in AD if PINK1 <= 1

grp <- rep("Other", nrow(exp_data))
grp[exp_data$condition == "Control" & exp_data$PINK1 >=  hi_thr] <- "High_in_Control"
grp[exp_data$condition == "AD"      & exp_data$PINK1 < low_thr] <- "Low_in_AD"

combined$PINK1_expr_group <- factor(grp, levels = c("Other","High_in_Control","Low_in_AD"))

# --- Highlight overlay (green under, orange on top) ---
emb <- Embeddings(combined, "umap") %>% as.data.frame()
colnames(emb)[1:2] <- c("UMAP_1","UMAP_2")
emb$group <- combined$PINK1_expr_group[rownames(emb)]

df_hilo <- emb %>% filter(group %in% c("High_in_Control","Low_in_AD"))
df_hi   <- df_hilo %>% filter(group == "High_in_Control")
df_lo   <- df_hilo %>% filter(group == "Low_in_AD")        

# ranges for right-side labels
xr <- range(df_hilo$UMAP_1); yr <- range(df_hilo$UMAP_2)
dx <- diff(xr); dy <- diff(yr)
x_lab <- xr[2] + 0.08*dx
y_mid <- mean(yr)

p <- ggplot() +
  # points (green below, orange above)
  geom_point(data = df_hi, aes(UMAP_1, UMAP_2, color = group), size = 0.5, alpha = 0.75) +
  geom_point(data = df_lo, aes(UMAP_1, UMAP_2, color = group), size = 0.5, alpha = 0.9, show.legend = TRUE) +
  scale_color_manual(values = c("High_in_Control" = "green", "Low_in_AD" = "orange"),
                     labels = c("High in Control", "Low in AD"), drop = TRUE) +
  # extend x to make room for labels
  coord_cartesian(xlim = c(xr[1], xr[2] + 0.12*dx), clip = "off") +
  coord_equal() +
  theme_classic(base_size = 12) +
  theme(
    plot.background  = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.background= element_rect(fill = "transparent", color = NA),
    legend.key       = element_rect(fill = "transparent", color = NA),
    plot.margin      = margin(5.5, 50, 5.5, 5.5)
  ) +
  guides(color = guide_legend(title = NULL, override.aes = list(size = 3, alpha = 1))) +
  labs(x = "UMAP_1", y = "UMAP_2") 
dir.create("Figures2", showWarnings = FALSE, recursive = TRUE)
ggsave("Figures2/pink_high_low_only_clear.tif", p, width = 5.8, height = 4, dpi = 600, bg = "transparent")
print(p)
###############################################################
# Set filter for alkabh3 high and pink1 High
###############################################################
#Define thresholds for high/low expression
# Rebuild exp_data with gene columns
exp_data <- FetchData(combined, vars = c("ALKBH3","PINK1","condition"))

# ---- Define thresholds (adjust as you like)
alkbh3_high_thresh <- 1
alkbh3_low_thresh  <- 1
pink1_high_thresh  <- 1
pink1_low_thresh   <- 1

# ---- Label ONLY the two groups; leave others as NA so they won't be plotted
exp_data <- exp_data %>%
  mutate(
    expr_group = dplyr::case_when(
      .data[["ALKBH3"]] >  alkbh3_high_thresh & .data[["PINK1"]] <  pink1_low_thresh ~ "ALKBH3",
      .data[["ALKBH3"]] <  alkbh3_low_thresh  & .data[["PINK1"]] >  pink1_high_thresh ~ "PINK1",
      TRUE ~ NA_character_
    )
  )

# Write back and keep only the two groups
combined$expr_group <- factor(exp_data$expr_group, levels = c("ALKBH3","PINK1"))
cells_keep <- rownames(combined@meta.data)[!is.na(combined$expr_group)]

# ---- UMAP (no 'Other')
cell_group <- DimPlot(
  combined, cells = cells_keep, group.by = "expr_group",
  cols = c("ALKBH3"="red", "PINK1"="blue"),
  shuffle = TRUE, pt.size = 0.25
) + ggtitle(NULL)
ggsave("Figures2/cell_group_alk_pink.png", plot = cell_group, width = 5, height = 4, dpi = 600)

cell_group


# ---- Violins (no 'Other')
Idents(combined) <- "expr_group"
p_alk <- VlnPlot(combined, features = "ALKBH3",
                 idents = c("ALKBH3","PINK1"),
                 pt.size = 0.1) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.title.x = element_blank()) +
  labs(title= NULL,x = NULL, y = "ALKBH3 (log-normalized)")
ggsave("Figures2/ALKBH3_by_expr.tif", plot = p_alk, width = 5, height = 4, dpi = 600)

p_pink <- VlnPlot(combined, features = "PINK1",
                  idents = c("ALKBH3","PINK1"),
                  pt.size = 0.1) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.title.x = element_blank()) +
  labs(title= NULL,x = NULL, y = "PINK1 (log-normalized)")
ggsave("Figures2/PINK1_by_expr.tif", plot = p_pink, width = 5, height = 4, dpi = 600)

p_pink
p_alk
############################################################## 
# 'condition' metadata exists with "AD" and "Control"
###############################################################
combined$condition <- ifelse(grepl("AD", combined$group), "AD", "Control")
Idents(combined) <- combined$condition
# Extract metadata and gene expression
exp_data <- FetchData(combined, vars = c("ALKBH3", "PINK1", "condition", "seurat_clusters"))
# Create groupings for easier filtering
exp_data$group <- combined$condition
# Only cells where ALKBH3 or PINK1 is expressed
alkbh3_expr <- subset(exp_data, ALKBH3 > 0)
pink1_expr <- subset(exp_data, PINK1 > 0)
expr_long <- exp_data %>%
  filter(ALKBH3 > 0 | PINK1 > 0) %>%
  pivot_longer(cols = c(ALKBH3, PINK1), names_to = "Gene", values_to = "Expression")
# Violin with p-value
v_alk_pink<-ggviolin(expr_long, x = "group", y = "Expression", fill = "group", facet.by = "Gene",
                     add = "boxplot", add.params = list(fill = "white")) + #stat_compare_means(method = "wilcox.test") +
  scale_fill_manual(values = c("Control" = "azure4", "AD" = "brown3"))+
  theme_minimal()+labs(x = "", y = "Expression Level")
ggsave("Figures2/Violin_ALK_PINK1.tif", plot = v_alk_pink, width = 5, height = 4, dpi = 600)

v_alk_pink

