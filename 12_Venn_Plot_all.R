library(GenomicRanges)
library(org.EcK12.eg.db) # 专门用于 E. coli K‐12 MG1655 基因组注释
library(dplyr)
library(VennDiagram)
library(grid)

# 读取数据
hg=read.csv("Homologues/Homologues_genes_SonicParanoid2.csv")
df1=read.csv("set1/Aeromonas_Integrated_Data.csv")
df2=read.csv("set2/Aeromonas_Integrated_Data.csv")
ess=read.csv("set1/ecoli_essential_genes_Keio.csv")

# 去除Ecoli蛋白编号的版本号（如 NP_418415.1 -> NP_418415）
hg$Ecoli <- sub("\\.[0-9]+$", "", hg$Ecoli)

# 构建映射
mapping_info <- AnnotationDbi::select(org.EcK12.eg.db,
                                      keys =hg$Ecoli,
                                      keytype = "REFSEQ",
                                      columns = "ALIAS")
# 筛选出以 "ECK" 开头的别名（符合系统编号格式）
mapping_info_filtered <- mapping_info %>%
  filter(grepl("^(ECK)[0-9]+", ALIAS))

# 在同源数据中添加大肠杆菌的ECK编号
idx=match(hg$Ecoli,mapping_info_filtered$REFSEQ)
hg$Ecoli_ECK=mapping_info_filtered$ALIAS[idx]

# 在同源数据中添加set1的Aeromonas的必须性信息
idx1=match(hg$Aeromonas,df1$gene_id)
hg$Aero1_essential_status=df1$essential_status_2[idx1]

# 在同源数据中添加E.coli的必须性信息
idx2=match(hg$Ecoli_ECK,ess$ECK_number)
hg$Ecoli_essential_status=ifelse(is.na(idx2),"uncertain","essential")

# 在同源数据中添加set2的Aeromonas的必须性信息
idx1=match(hg$Aeromonas,df2$gene_id)
hg$Aero2_essential_status=df2$essential_status_2[idx1]

# 保存数据
write.csv(hg,"Homologues/Homologues_Essential_Status_all.csv",row.names = FALSE)



# 1. 提取三组“essential”基因的行号
set1_idx <- which(hg$Aero1_essential_status == "essential")
set2_idx <- which(hg$Aero2_essential_status == "essential")
ecoli_idx <- which(hg$Ecoli_essential_status   == "essential")

# 2. 计算各个区域大小
a1   <- length(set1_idx)
a2   <- length(set2_idx)
a3   <- length(ecoli_idx)
n12  <- length(intersect(set1_idx, set2_idx))
n23  <- length(intersect(set2_idx, ecoli_idx))
n13  <- length(intersect(set1_idx, ecoli_idx))
n123 <- length(Reduce(intersect, list(set1_idx, set2_idx, ecoli_idx)))

# 3. 绘制三集韦恩图
grid.newpage()
venn_plot=draw.triple.venn(
  area1    = a1,
  area2    = a2,
  area3    = a3,
  n12      = n12,
  n23      = n23,
  n13      = n13,
  n123     = n123,
  category = c("Aero Set1", "Aero Set2", "E. coli"),
  fill     = c("steelblue", "salmon", "orange"),
  alpha    = c(0.5, 0.5, 0.5),
  cex      = 2,            # 圆内数字大小
  cat.cex  = 1.5,          # 图例文字大小
  cat.pos  = c(-20, 20, 180),  # 三个类别标签的位置
  cat.dist = c(0.05, 0.05, 0.05),# 三个类别标签与圆的距离
  scaled   = TRUE,
  lty = "blank",
)

pdf("Figures/venn_plot.pdf", width = 12, height = 9)
# 使用 grid.draw 绘制 grid 图形对象
grid.draw(venn_plot)
dev.off()
