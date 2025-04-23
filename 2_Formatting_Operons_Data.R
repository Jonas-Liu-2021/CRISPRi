library(tidyverse)
library(GenomicRanges)
library(readr)
library(data.table)

## 格式化Operon Mapper获取的操纵子文件
# 读取操纵子列表文件，填充缺失列并保留NA字符串为缺失值
operons_df <- read.delim("Data/list_of_operons_2286949.txt", header = TRUE, 
                         sep = "\t", fill = TRUE, stringsAsFactors = FALSE, 
                         na.strings = "NA")


# 使用fill函数将Operon列的缺失值向下填充，确保每行都有正确的操纵子编号
operons_df <- fill(operons_df, Operon)

# 移除没有基因信息的行（这些行仅包含Operon编号且IdGene等列为空）
operons_df <- operons_df[ !is.na(operons_df$IdGene) & operons_df$IdGene != "", ]


# 将位置列转换为数值型（整数），操纵子编号也转为整数类型
operons_df$Operon   <- as.integer(operons_df$Operon)
operons_df$PosLeft  <- as.integer(operons_df$PosLeft)
operons_df$postRight <- as.integer(operons_df$postRight)


# 生成所有可能的基因ID（按顺序应该是连续的）
all_gene_ids <- sprintf("PPGIJJGB_%05d", 1:4597)

# 找出缺失的基因ID
missing_gene <- setdiff(all_gene_ids, operons_df$IdGene)

# 打印缺失的基因
print(missing_gene) # 少了一个PPGIJJGB_03094，这个是tmRNA，Operon-Mapper不处理tmRNA，验证后发现该基因也不与相邻基因共用一个操纵子

### 补充缺少的PPGIJJGB_03094的那一行
new_row = data.frame(Operon=2546,IdGene="PPGIJJGB_03094",Type="tmRNA",COGgene=NA,PosLeft=3401328,
                     postRight=3401687, Strand="+",Function=NA)
new_row$Operon   <- as.integer(new_row$Operon)    # 转为int
new_row$COGgene   <- as.character(new_row$COGgene)    # 转为字符
new_row$PosLeft  <- as.integer(new_row$PosLeft)     # 转为整数
new_row$postRight <- as.integer(new_row$postRight)    # 转为整数
new_row$Function <- as.character(new_row$Function)  # 转为字符

# 定义插入的位置
insert_pos <- 3094

# 使用 rbind 将数据框分成两部分，并在中间插入新行
operons_df_new <- rbind(operons_df[1:(insert_pos - 1), ], new_row, operons_df[insert_pos:nrow(operons_df), ])


### 将操纵子数据与DESeq2的结果合并
# 读取DESeq2数据
deg_path <- "Data/DESeq2_results_filtered.csv"
deg <- fread(deg_path, header = TRUE, sep = ",")


# 确保 logFC 是数值型
deg[, log2FoldChange := as.numeric(log2FoldChange)]

# 创建 logFC 映射表
logFC_dict <- setNames(deg$log2FoldChange, deg$V1)

# 创建 GRanges 对象,并添加必要信息
gr <- GRanges(
  seqnames = "chromosome",  # 假设所有基因来自同一条染色体，若有多个染色体，需要修改此列
  ranges = IRanges(start = operons_df_new$PosLeft, end = operons_df_new$postRight),
  gene_id = operons_df_new$IdGene,
  logFC=ifelse(operons_df_new$IdGene %in% names(logFC_dict), logFC_dict[operons_df_new$IdGene],NA)
)

seqinfo(gr) <- Seqinfo(
  seqnames = "chromosome",  # 必须明确指定染色体名字
  seqlengths = c("chromosome" = 4962974),  # 长度信息
  isCircular = c("chromosome" = TRUE)  # 必须按 `seqnames` 长度匹配
)

# 添加 Operon,Type和Function 三列到 gr 的 metadata 中
gr$Type=operons_df_new$Type
gr$Function=operons_df_new$Function
gr$Operon=operons_df_new$Operon

# 添加 Strand 到 gr 的核心数据 strand中
strand(gr) <- Rle(operons_df_new$Strand)

# 添加padj信息到gr
idx=match(gr$gene_id,deg$V1)
gr$padj=deg$padj[idx]


# 在 GRanges 的 metadata 中添加logFC < -2的列
gr$fc_lt_2 <- ifelse(is.na(gr$logFC),FALSE,gr$logFC< -2)

### 从EMAPPER得到的信息添加GO，KEGG, COG编号和基因名信息
# 读取注释文件（假设文件名为 annotation.csv）
annot_data <- read_csv("Data/out.emapper.annotations.csv") %>%
  dplyr::select(gene = query, GOs, KEGG_Pathway,COG_category,Preferred_name)

# 把注释数据的GO和KEGG号加入到gr中
idx=match(gr$gene_id,annot_data$gene)
gr$GOs = annot_data$GOs[idx]
gr$GOs[gr$GOs=="-"]=NA
gr$KEGG_Pathway=annot_data$KEGG_Pathway[idx]
gr$KEGG_Pathway[gr$KEGG_Pathway=="-"]=NA
gr$COG_category=annot_data$COG_category[idx]
gr$COG_category[gr$COG_category=="-"]=NA
gr$gene_name=annot_data$Preferred_name[idx]
gr$gene_name[gr$gene_name=="-"]=NA

# 打印 GRanges 对象查看
gr

# 保存该 GRanges 对象以供后续分析使用
saveRDS(gr,"Output_RDS/GRanges_log2FC_Operon.rds")


## 提取多顺反子操纵子 （多基因操纵子）
# 统计每个 Operon 的频次
operon_counts <- table(gr$Operon)

# 筛选出出现多次的 Operon 编号
multi_operon_ids <- names(operon_counts)[operon_counts > 1]

# 从 gr 中提取这些基因
gr_multi <- gr[gr$Operon %in% multi_operon_ids]

# 将 gr_multi 转换为 data.frame
gr_multi_df <- as.data.frame(gr_multi)
gr_multi_df <- gr_multi_df[, c("seqnames", 
                               "gene_id",
                               "gene_name",
                               "Type",
                               "Operon",
                               "start", 
                               "end", 
                               "width", 
                               "strand", 
                               "logFC", 
                               "fc_lt_2")] # 调整列的顺序

# 写出 CSV 文件到 "Polycistron" 文件夹中
write.csv(gr_multi_df, file = "Polycistron/Polycistron.csv", row.names = FALSE)


