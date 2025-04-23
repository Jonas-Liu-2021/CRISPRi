library(GenomicRanges)
library(org.EcK12.eg.db) # 专门用于 E. coli K‐12 MG1655 基因组注释
library(dplyr)
library(VennDiagram)
library(grid)

# 读取数据
hg=read.csv("Homologues/Homologues_genes_SonicParanoid2.csv")
df=read.csv("Data/Aeromonas_Integrated_Data.csv")
ess=read.csv("Data/ecoli_essential_genes_Keio.csv")

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

# 在同源数据中添加Aeromonas的必须性信息
idx1=match(hg$Aeromonas,df$gene_id)
hg$Aero_essential_status=df$essential_status_2[idx1]

# 在同源数据中添加E.coli的必须性信息
idx2=match(hg$Ecoli_ECK,ess$ECK_number)
hg$Ecoli_essential_status=ifelse(is.na(idx2),"uncertain","essential")

# 保存数据
write.csv(hg,"Homologues/Homologues_Essential_Status.csv",row.names = FALSE)

### 绘制共同必须基因的韦恩图
# 计算共同必需基因（两侧均为 "essential"）
common <- sum(hg$Aero_essential_status == "essential" & hg$Ecoli_essential_status == "essential", na.rm = TRUE)

# 计算 Aeromonas 独有必需：仅 Aero 为 essential，而 Ecoli 不是 essential
aero_only <- sum(hg$Aero_essential_status == "essential" & hg$Ecoli_essential_status != "essential", na.rm = TRUE)

# 计算 Ecoli 独有必需：仅 Ecoli 为 essential，而 Aero 不是 essential
ecoli_only <- sum(hg$Aero_essential_status != "essential" & hg$Ecoli_essential_status == "essential", na.rm = TRUE)

# 总数分别为：
aero_total <- common + aero_only  # Aeromonas 总必需基因数
ecoli_total <- common + ecoli_only  # E. coli 总必需基因数

# 输出计算结果
cat("Aeromonas必需基因总数:", aero_total, "\n")
cat("E. coli必需基因总数:", ecoli_total, "\n")
cat("共同必需基因数:", common, "\n")
cat("Aeromonas独有必需基因数:", aero_only, "\n")
cat("E. coli独有必需基因数:", ecoli_only, "\n")

# 绘制韦恩图
# 参数说明：
# - area1、area2 分别为两侧的必需基因总数
# - cross.area 为共同必需基因数
# - category 设置为空字符串，且 cat.cex=0 隐藏类别文字
venn.plot <- draw.pairwise.venn(
  area1 = aero_total,
  area2 = ecoli_total,
  cross.area = common,
  category = c("", ""),
  print.mode = "raw",   # 只显示数字
  fill = c("skyblue", "pink"),
  lty = "blank",
  cex = 2,
  cat.cex = 0
)


pdf("Figures/venn_plot.pdf", width = 12, height = 9)
# 使用 grid.draw 绘制 grid 图形对象
grid.draw(venn.plot)
dev.off()


# 绘制外圈表示所有同源基因（共 2122 个）的集合
# 这里假设韦恩图绘制在 npc 坐标系中，中心在 (0.5, 0.5)
# 可适当调整半径以使外圈包住两个小圈
# grid.circle(x = 0.5, y = 0.5, r = 0.55, 
#             gp = gpar(lwd = 2, col = "black", fill = NA))
# 
# # 在外圈上方添加数字 "2122"
# grid.text("2122", x = 0.5, y = 0.95, gp = gpar(cex = 2))
