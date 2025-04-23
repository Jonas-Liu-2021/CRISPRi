## 此代码是为了绘制有颜色分类的火山图
library(ggplot2)
library(tidyverse)
library(tidyr)

# 读取csv文件
df=read.csv("Data/Aeromonas_Integrated_Data.csv")
# 定义COG字母对应的类别含义
COG_meaning <- data.frame(
  COG_category = c("J", "A", "K", "L", "B", "D", "Y", "V", "T", "M", 
                   "N", "Z", "W", "U", "O", "C", "G", "E", "F", "H", "I", 
                   "P", "Q", "R", "S", "X"),
  Function = c("Translation, ribosomal structure and biogenesis",
               "RNA processing and modification",
               "Transcription",
               "Replication, recombination and repair",
               "Chromatin structure and dynamics",
               "Cell cycle control, cell division, chromosome partitioning",
               "Nuclear structure",
               "Defense mechanisms",
               "Signal transduction mechanisms",
               "Cell wall/membrane/envelope biogenesis",
               "Cell motility",
               "Cytoskeleton",
               "Extracellular structures",
               "Intracellular trafficking, secretion, and vesicular transport",
               "Posttranslational modification, protein turnover, chaperones",
               "Energy production and conversion",
               "Carbohydrate transport and metabolism",
               "Amino acid transport and metabolism",
               "Nucleotide transport and metabolism",
               "Coenzyme transport and metabolism",
               "Lipid transport and metabolism",
               "Inorganic ion transport and metabolism",
               "Secondary metabolites biosynthesis, transport and catabolism",
               "General function prediction only",
               "Function unknown",
               "Mobilome: prophages, transposons")
)

# 首先，创建一个映射向量以便后续快速替换
COG_map <- setNames(COG_meaning$Function, COG_meaning$COG_category)

# 将COG_category的有多个功能的基因拆分成几行
df=df %>% filter(!is.na(COG_category)) %>%
  separate_rows(COG_category,sep = "") %>%
  filter(COG_category != "")


# 对df进行转换，得到COG_anno新列
df <- df %>%
  mutate(
    COG_anno = COG_map[as.character(COG_category)]
  )

# 计算负对数调整p值，处理 padj=0 的情况（避免无限大）
df$neg_log10_padj <- -log10(df$padj + 1e-300)  # 添加极小值防止 log(0)

# 如果功能标签未转换为factor，转换一下以保证颜色自动分类
df$COG_anno <- as.factor(df$COG_anno)

# 去掉"logFC", "neg_log10_padj" 和"COG_anno"中有NA的行
df <- df %>% 
  filter(
    !is.na(logFC),
    !is.na(neg_log10_padj),
    !is.na(COG_anno)
  )

# 共有20种功能分类，给出一个20色自定义方案：
# 这个比较美观
my_colors <- c(
  "#1f78b4", "#33a02c", "#e31a1c", "#ff7f00",
  "#6a3d9a", "#b15928", "#a6cee3", "#b2df8a",
  "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99",
  "#8dd3c7", "#ffffb3", "#bebada", "#fb8072",
  "#80b1d3", "#fdb462", "#b3de69", "#fccde5"
)
#下面这个对比更加明显
# my_colors <- c(
#   "#e6194b", "#3cb44b", "#4363d8", "#f58231",
#   "#911eb4", "#46f0f0", "#f032e6", "#bcf60c",
#   "#fabebe", "#008080", "#e6beff", "#9a6324",
#   "#fffac8", "#800000", "#aaffc3", "#808000",
#   "#ffd8b1", "#000075", "#808080", "#000000"
# )
# 绘制火山图，按照COG功能分类设置颜色
p <- ggplot(df, aes(x = logFC, y = neg_log10_padj, color = COG_anno, shape = COG_anno)) +
  geom_point(position = position_jitter(width = 0, height = 0.1),alpha = 1, size = 2) +
  # 添加垂直虚线
  geom_vline(
    xintercept = -2,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.5
  ) +
  theme_minimal() +
  labs(
    title = "Volcano Plot by COG Functional Categories",
    x = expression(log[2]("Fold Change")),
    y = expression(-log[10]("Adjusted p-value")),
    color = "COG Category", # 图例名称
    shape = "COG Category"  # 形状图例名称
  ) +
  # 设置颜色主题
  scale_color_manual(values = my_colors )+
  scale_shape_manual(values = 1:length(unique(df$COG_anno))) # 为每个COG分类分配不同形状


# 显示图形
print(p)

# 保存为PDF
pdf("Figures/volcano_plot_COG_colored.pdf", width = 18, height = 10)
print(p)
dev.off()

