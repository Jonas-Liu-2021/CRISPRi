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
  "#1f78b4", "#3cb44b", "#e31a1c", "#ff7f00",
  "#6a3d9a", "#b15928", "#a6cee3", "#fffac8",
  "#fb9a99", "#008080", "#E10EF0", "#800000",
  "#aaffc3", "#EAFF03", "#09FF00", "#fb8072",
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
# 生成分组列（Non-essential genes + 各个 COG 功能）
df2 <- df %>%
  mutate(
    Group = if_else(
      essential_status_2 == "essential",
      as.character(COG_anno),
      "Non-essential genes"
    )
  ) %>%
  mutate(Group = factor(Group, levels = c("Non-essential genes", levels(COG_anno))))

df2 <- df2 %>%
  mutate(
    # 其他因子顺序不变，把 "Function unknown" 放到最后
    Group = fct_relevel(Group, "Function unknown", after = Inf)
  )

#################################################
# 初始版本，点有形状，图注在右边，无xy轴，有网格
#################################################

# 颜色映射：Non-essential 灰色 + 你的 my_colors
col_vals <- c(
  "Non-essential genes" = "grey80",
  setNames(my_colors, levels(df2$COG_anno))
)

# 形状映射：Non-essential 用 16，essential 的类别从 0–25 中取不含 16/19 的前 N 个
available_shapes <- setdiff(0:25, c(16, 18,19,20))
shape_choices   <- available_shapes[seq_len(length(levels(df2$COG_anno)))]
shape_vals <- c(
  "Non-essential genes" = 16,
  setNames(shape_choices, levels(df2$COG_anno))
)

# 绘图
p <- ggplot(df2, aes(x = logFC, y = neg_log10_padj, color = Group, shape = Group)) +
  geom_point(position = position_jitter(width = 0, height = 0.1),
             size = 2, alpha = 0.8) +
  geom_vline(xintercept = -2, linetype = "dashed", color = "grey50", size = 0.5) +
  scale_color_manual(
    name   = "Gene Type / COG Category",
    values = col_vals
  ) +
  scale_shape_manual(
    name   = "Gene Type / COG Category",
    values = shape_vals
  ) +
  labs(
    title = "Volcano Plot by COG Functional Categories",
    x     = expression(log[2]("Fold Change")),
    y     = expression(-log[10]("Adjusted p-value"))
  ) +
  theme_minimal(base_size = 12)

print(p)




# 保存为PDF
pdf("Figures/volcano_plot_COG_colored.pdf", width = 18, height = 10)
print(p)
dev.off()


#############################################################
# V2.0版本，点只有颜色，图注在下，有xy轴，无网格
#############################################################
# 自定义颜色向量（Non-essential + 各 COG 功能）
col_vals <- c(
  "Non-essential genes" = "grey80",
  setNames(my_colors, levels(df2$COG_anno))
)

# 开始绘图
p <- ggplot(df2, aes(x = logFC, y = neg_log10_padj, color = Group)) +
  # 实心圆，适当调大 size，可根据需要再微调
  geom_point(shape = 16,
             size = 3,
             alpha = 0.9,
             position = position_jitter(width = 0, height = 0.1)) +
  # 参考线
  geom_vline(xintercept = -2,
             linetype   = "dashed",
             color      = "grey50",
             linewidth       = 0.5) +
  # 水平线：padj = 0.05
  geom_hline(yintercept = -log10(0.05),
             linetype   = "dashed",
             color      = "grey50",
             linewidth       = 0.5) +

  # 颜色映射
  scale_color_manual(
    name   = "Gene Type / COG Category",
    values = col_vals,
    guide  = guide_legend(ncol = 2, byrow = TRUE)
  ) +
  labs(
    title = "Volcano Plot by COG Functional Categories",
    x     = expression(log[2]*FC),
    y     = expression(-log[10]*(p[adj]))
  ) +
  # 主题设置：透明背景 + 黑色文字 + 下方图例
  theme_minimal(base_size = 14) +
  theme(
    # 去掉网格
    panel.grid  = element_blank(),
    # 坐标轴、标题黑色
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(colour = "black"),
    plot.title       = element_text(color = "black", size = 16, face = "bold", hjust = 0.5),
    axis.title       = element_text(color = "black", size = 14),
    axis.text        = element_text(color = "black", size = 12),
    # 图例置底、透明键与背景
    legend.position   = "bottom",
    legend.title = element_blank(),
    legend.key       = element_rect(fill = NA, color = NA),
    legend.spacing.y   = unit(0.05, "cm"),   # 行与行之间仅 0.05cm
    legend.key.height   = unit(0.05, "cm"),   
    legend.text      = element_text(color = "black", size = 9)
  )

# 显示
print(p)

# 保存为透明背景的 PDF
pdf("Figures/volcano_plot_COG_colored2.pdf",
    width = 8, height = 11)
print(p)
dev.off()

