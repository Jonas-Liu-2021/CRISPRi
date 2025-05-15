# ---- 加载包 ----
library(ggrepel)
library(pheatmap)
library(dplyr)

# ---- 读入数据 ----
set1 <- read.csv("set1/set1.csv", row.names = 1,check.names = FALSE)
set2 <- read.csv("set2/set2.csv", row.names = 1,check.names = FALSE)

# ---- 筛选共同基因 ----
common_genes <- intersect(rownames(set1), rownames(set2))
set1 <- set1[common_genes, ]
set2 <- set2[common_genes, ]

# ---- 重新排序列，确保对照组+实验组顺序 ----
set1 <- set1[, c("1_R1_001", "2_R1_001", "3_R1_001", "4_R1_001", "5_R1_001", "6_R1_001")]
set2 <- set2[, c("7_R1_001", "8_R1_001", "9_R1_001", "10_R1_001", "11_R1_001", "12_R1_001")]

# ---- 合并数据 ----
expr_data <- cbind(set1, set2)

# ---- log2(count+1) 转换 ----
log_expr <- log2(expr_data + 1)

# ---- 定义分组信息 ----
group_labels <- c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment",
                  "Control", "Control", "Control", "Treatment", "Treatment", "Treatment")

# ---- PCA分析 ----
pca <- prcomp(t(log_expr), scale. = TRUE)
pca_data <- as.data.frame(pca$x)
pca_data$Group <- factor(group_labels)
pca_data$Sample <- colnames(log_expr)

# ---- 绘制 PCA ----
ggplot(pca_data, aes(PC1, PC2, color=Group, label=Sample)) +
  geom_point(size=2) +
  geom_text_repel(size=3, show.legend = FALSE) +  # 👈 自动避让标签
  labs(title="PCA of CRISPRi Data",
       x=paste0("PC1 (", round(summary(pca)$importance[2,1]*100,1), "%)"),
       y=paste0("PC2 (", round(summary(pca)$importance[2,2]*100,1), "%)")) +
  theme_minimal()+
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_blank()
  )

# ---- 热图（按基因标准化） ----
ann_colors <- list(Group = c(Control = "#1f77b4", Treatment = "#ff7f0e"))
annotation_col <- data.frame(Group = factor(group_labels))
rownames(annotation_col) <- colnames(log_expr)

pheatmap(log_expr,
         scale = "row",
         clustering_method = "complete",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_rownames = FALSE,
         fontsize_col = 10,
         main = "Hierarchical Clustering Heatmap")

# ---- 样本-样本距离热图 ----
# 使用欧氏距离
dist_matrix <- as.matrix(dist(t(log_expr), method = "euclidean"))

pheatmap(dist_matrix,
         clustering_method = "complete",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         display_numbers = FALSE,
         fontsize_number = 8,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Sample-Sample Distance Heatmap")
