library(dplyr)
library(VennDiagram)
library(grid)
library(tidyr)
library(ggplot2)

# 读取数据
df1=read.csv("set1/Aeromonas_Integrated_Data.csv")
df2=read.csv("set2/Aeromonas_Integrated_Data.csv")

# 分别找出set1和set2独占的必须基因
# 提取必需基因子集
ess1 <- df1 %>% filter(essential_status_2 == "essential")
ess2 <- df2 %>% filter(essential_status_2 == "essential")

# 计算共同与各自独占
#    （a）共同的基因ID
common_ids <- intersect(ess1$gene_id, ess2$gene_id)
#    （b）只在 set1 出现的
unique1_ids <- setdiff(ess1$gene_id, ess2$gene_id)
#    （c）只在 set2 出现的
unique2_ids <- setdiff(ess2$gene_id, ess1$gene_id)

# 构造三个 data.frame
common_essential   <- ess1   %>% filter(gene_id %in% common_ids)
common_essential = common_essential[,c("gene_id","gene_name","Type","Function","COG_category")]
unique1_essential  <- ess1   %>% filter(gene_id %in% unique1_ids)
unique1_essential = unique1_essential[,c("gene_id","gene_name","Type","Function","COG_category")]
unique2_essential  <- ess2   %>% filter(gene_id %in% unique2_ids)
unique2_essential = unique2_essential[,c("gene_id","gene_name","Type","Function","COG_category")]

# 检查结果
length(common_ids)   # 共有多少
length(unique1_ids)  # set1 独占多少
length(unique2_ids)  # set2 独占多少

# 保存共同必须基因和各自独占的必需基因

# 绘制韦恩图
# 基本计数
n1 <- nrow(ess1)                # set1 中 essential 总数
n2 <- nrow(ess2)                # set2 中 essential 总数
n12 <- length(common_ids)       # 两者共有 essential 数

# 画图
grid.newpage()
venn.plot <- draw.pairwise.venn(
  area1      = n1,
  area2      = n2,
  cross.area = n12,
  category   = c("Set1", "Set2"),
  
  # 美化参数
  fill       = c("steelblue", "salmon"),
  alpha      = c(0.5, 0.5),
  lty        = "blank",      # 去掉边框线
  cex        = 2,            # 圆内数字大小
  cat.cex    = 1.5,          # 图例文字大小
  cat.pos    = c(180, 180),   # 图例文字位置
  cat.dist   = c(0.05, 0.05),  # 图例与圆的距离
  scaled     = TRUE,
  inverted   = TRUE,      # 翻转左右位置
)

pdf("Figures/essential_venn_plot.pdf", width = 12, height = 9)
# 使用 grid.draw 绘制 grid 图形对象
grid.draw(venn.plot)
dev.off()

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

# 将含多个COG字母的行拆分开，每个字母一行
common_essential <- common_essential %>%
  mutate(COG_category = strsplit(as.character(COG_category), "")) %>%
  unnest(COG_category)

unique1_essential <- unique1_essential %>%
  mutate(COG_category = strsplit(as.character(COG_category), "")) %>%
  unnest(COG_category)

unique2_essential <- unique2_essential %>%
  mutate(COG_category = strsplit(as.character(COG_category), "")) %>%
  unnest(COG_category)


# 统计每个子集在各 COG_category 下的基因数
common_cnt <- common_essential %>%
  count(COG_category) %>%
  rename(common = n)

set1_cnt <- unique1_essential %>%
  count(COG_category) %>%
  rename(set1 = n)

set2_cnt <- unique2_essential %>%
  count(COG_category) %>%
  rename(set2 = n)

# 准备长表
COG_count <- full_join(common_cnt, set1_cnt, by = "COG_category") %>%
  full_join(set2_cnt, by = "COG_category") %>%
  replace_na(list(common = 0, set1 = 0, set2 = 0)) %>%
  left_join(COG_meaning, by = "COG_category") %>%
  mutate(total = common + set1 + set2)

COG_count_long <- COG_count %>%
  pivot_longer(
    cols      = c(common, set1, set2),
    names_to  = "subset",
    values_to = "Count"
  ) %>%
  # 强制指定 stacking 顺序
  mutate(subset = factor(subset, levels = c("common", "set1", "set2"))) 

# 绘图 
p = ggplot(COG_count_long, aes(
  x    = reorder(Function, total),
  y    = Count,
  fill = subset
)) +
  geom_bar(stat = "identity", width = 0.7) +
  # 只给 Count>0 的那部分画文字 
  geom_text(
    data     = filter(COG_count_long, Count > 0),
    aes(label = Count),
    position = position_stack(vjust = 0.5),
    size     = 3,
    color    = "white"
  ) +
  coord_flip() +
  labs(
    x     = "COG Functional Category",
    y     = "Number of Genes",
    fill  = "Subset",
    title = "Gene Count per COG Functional Category"
  ) +
  scale_fill_manual(
    values = c(common = "steelblue", set1 = "skyblue", set2 = "salmon"),
    labels = c(common = "Common", set1 = "Only Set1", set2 = "Only Set2")
  ) +
  theme_minimal() +
  theme(

    axis.text.y   = element_text(size = 12, margin = margin(r = -20)),
    axis.text.x   = element_text(size = 11),
    axis.title    = element_text(size = 12, face = "bold"),
    plot.title    = element_text(size = 14, face = "bold", hjust = 0),
    legend.title  = element_text(size = 12),
    legend.text   = element_text(size = 11)
  )

p
# 保存图像信息
pdf("Figures/COG_category_difference.pdf", width = 16, height = 9)
print(p)
dev.off()
