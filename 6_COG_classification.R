library(tidyverse)
library(dplyr)

# 读取数据
df=read.csv("Data/Aeromonas_Integrated_Data.csv")

# 筛选必须基因
df=df[df$essential_status_2=="essential", c("gene_id","essential_status_2","COG_category")]

# 去掉df中有NA的行
df=df[!is.na(df$COG_category),]  # 留下的都是CDS

# 将含多个COG字母的行拆分开，每个字母一行
df_split <- df %>%
  mutate(COG_category = strsplit(as.character(COG_category), "")) %>%
  unnest(COG_category)

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

# 统计每个COG类别数量
COG_count <- df_split %>%
  group_by(COG_category) %>%
  summarise(Count = n()) %>%
  left_join(COG_meaning, by = "COG_category") %>%
  arrange(desc(Count))

# 绘制条形图(有数字)
p=ggplot(COG_count, aes(x = reorder(Function, Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = Count),hjust = -0.4,color = "black",size = 4) +
  coord_flip() +
  labs(x = "COG Functional Category", y = "Number of Genes", title = "Gene Count per COG Functional Category") +
  theme_minimal()+
  theme(
    axis.text.y = element_text(size = 12, margin = margin(r = -20)),  # 标签靠近柱子
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0)
  )

p


# 保存图像和统计信息
pdf("COG/COG_category.pdf", width = 12, height = 9)
print(p)
dev.off()

# 绘制条形图(无数字)
p=ggplot(COG_count, aes(x = reorder(Function, Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(x = "COG Functional Category", y = "Number of Genes", title = "Gene Count per COG Functional Category") +
  theme_minimal()+
  theme(
    axis.text.y = element_text(size = 12, margin = margin(r = -20)),  # 标签靠近柱子
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0)
  )

p


# 保存图像和统计信息
pdf("COG/COG_category(no_count).pdf", width = 12, height = 9)
print(p)
dev.off()

write.csv(COG_count,"COG/COG_category.csv")
