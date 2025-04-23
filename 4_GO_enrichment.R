library(clusterProfiler)
library(GenomicRanges)
library(GO.db)
library(KEGGREST)
library(AnnotationDbi)



# 读取RDS文件的整合数据
gr=readRDS("Output_RDS/GRanges_Integrated_Data_2.rds")


# --------------------------
# 构建GO映射表
# --------------------------
# 提取基因ID和对应的GO注释
gene_ids <- gr$gene_id      # 基因ID向量
go_terms <- gr$GOs         # 对应的逗号分隔的GO字符串

# 对每个基因拆分GO编号，生成映射列表
# 使用空字符向量替代 NA
gene2go_list <- lapply(go_terms, function(x) {
  if (is.na(x)) {
    character(0)
  } else {
    unlist(strsplit(x, split = ","))
  }
})

###--------------------------
# 删除废旧过时的GO条目
###--------------------------
# 所有 GO ID
all_go_ids <- keys(GO.db)

# # 查看是否是过时条目的关键方法
# is_obsolete <- function(go_id) {
#   term <- Term(go_id)
#   is.na(term) || grepl("obsolete", term, ignore.case = TRUE)
# }
# 
# # 用 lapply 检查所有 GO ID 是否过时（建议小批量做）
# obsolete_flags <- sapply(all_go_ids, is_obsolete)
# 
# # 得到所有过时的 GO ID
# obsolete_go_ids <- all_go_ids[obsolete_flags] # 结果发现GO.db似乎并没有保留过时的GO条目

# 清理 emapper 提供的 GO 注释 （取emapper提供的GO注释与GO.db中所有的有效GO条目取交集）
gene2go_list_clean <- lapply(gene2go_list, function(go_vec) {
  go_vec[go_vec %in% all_go_ids]
})

# 可选：只保留含有至少一个有效 GO 的基因
# gene2go_list_clean <- gene2go_list_clean[lengths(gene2go_list_clean) > 0]



# 构建映射表，注意只保留有GO注释的基因
gene_counts <- sapply(gene2go_list_clean, length)
gene2go_df <- data.frame(
  GO = unlist(gene2go_list_clean),
  gene = rep(gene_ids, gene_counts)
)

# -----------------------------
# 构建TERM2NAME映射表：GO编号 -> GO描述
# -----------------------------
unique_go_ids <- unique(gene2go_df$GO)
# Term() 函数从GO.db中查询每个GO编号对应的简短描述
go_desc <- Term(unique_go_ids)
go2desc_df <- data.frame(term = unique_go_ids, 
                         name = go_desc, 
                         stringsAsFactors = FALSE)

# # 手动注释因过时等缺失的GO name
# go2desc_df$name[go2desc_df$term == "GO:0034645"] <- "macromolecule biosynthetic process"
# go2desc_df$name[go2desc_df$term == "GO:0051186"] <- "cofactor metabolism"



# -----------------------------
# 定义目标基因集和背景基因集
# -----------------------------
# 目标基因集：essential_status 为 "essential" 的基因
essential_genes <- gr$gene_id[gr$essential_status == "essential"]
# 背景基因集：所有基因
background_genes <- gr$gene_id

# -----------------------------
# GO 富集分析
# -----------------------------
ego <- enricher(gene          = essential_genes,
                universe      = background_genes,
                TERM2GENE     = gene2go_df,
                TERM2NAME     = go2desc_df,
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2)

# 查看富集结果
head(ego)

p_dot=dotplot(ego, showCategory = 20) 
p_bar=barplot(ego, showCategory = 20)

p_dot
p_bar
pdf("GO/dotplot_ego.pdf", width = 12, height = 9)
print(p_dot)
dev.off()

pdf("GO/barplot_ego.pdf", width = 12, height = 9)
print(p_bar)
dev.off()
# cnetplot(ego, showCategory = 10)

# goplot(ego)
# 保存csv结果
ego_df=as.data.frame(ego)
write.csv(ego_df,"GO/ego.csv",row.names = FALSE)



###--------------------------------------
# 将数据分为BP,CC,MF三大类分别进行GO富集
###--------------------------------------
# 1. 获取每个 GO term 的 ontology 信息
go_ont <- AnnotationDbi::select(GO.db, 
                                keys = unique(gene2go_df$GO), 
                                columns = "ONTOLOGY", 
                                keytype = "GOID")

# 根据 ontology 信息划分 GO ID
BP_ids <- go_ont$GOID[go_ont$ONTOLOGY == "BP"]
CC_ids <- go_ont$GOID[go_ont$ONTOLOGY == "CC"]
MF_ids <- go_ont$GOID[go_ont$ONTOLOGY == "MF"]

# 2. 根据 GO 类型过滤 gene2go_df 与 go2desc_df
gene2go_df_BP <- gene2go_df[gene2go_df$GO %in% BP_ids, ]
gene2go_df_CC <- gene2go_df[gene2go_df$GO %in% CC_ids, ]
gene2go_df_MF <- gene2go_df[gene2go_df$GO %in% MF_ids, ]

go2desc_df_BP <- go2desc_df[go2desc_df$term %in% BP_ids, ]
go2desc_df_CC <- go2desc_df[go2desc_df$term %in% CC_ids, ]
go2desc_df_MF <- go2desc_df[go2desc_df$term %in% MF_ids, ]

# 3. 分别进行 GO 富集分析
ego_BP <- enricher(gene          = essential_genes,
                   universe      = background_genes,
                   TERM2GENE     = gene2go_df_BP,
                   TERM2NAME     = go2desc_df_BP,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2)

ego_CC <- enricher(gene          = essential_genes,
                   universe      = background_genes,
                   TERM2GENE     = gene2go_df_CC,
                   TERM2NAME     = go2desc_df_CC,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2)

ego_MF <- enricher(gene          = essential_genes,
                   universe      = background_genes,
                   TERM2GENE     = gene2go_df_MF,
                   TERM2NAME     = go2desc_df_MF,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2)

# 查看各分类的富集结果（可选）
head(ego_BP)
head(ego_CC)
head(ego_MF)

# 4. 绘制富集结果图：dotplot 与 barplot
# BP 类
p_dot_BP <- dotplot(ego_BP, showCategory = 20)
p_bar_BP <- barplot(ego_BP, showCategory = 20)

# CC 类
p_dot_CC <- dotplot(ego_CC, showCategory = 20)
p_bar_CC <- barplot(ego_CC, showCategory = 20)

# MF 类
p_dot_MF <- dotplot(ego_MF, showCategory = 20)
p_bar_MF <- barplot(ego_MF, showCategory = 20)

# 在 R 中显示图形
p_dot_BP
p_bar_BP
p_dot_CC
p_bar_CC
p_dot_MF
p_bar_MF

# 5. 保存图形为 PDF 文件
pdf("GO/Classify_before_GO/dotplot_ego_BP.pdf", width = 12, height = 9)
print(p_dot_BP)
dev.off()

pdf("GO/Classify_before_GO/barplot_ego_BP.pdf", width = 12, height = 9)
print(p_bar_BP)
dev.off()

pdf("GO/Classify_before_GO/dotplot_ego_CC.pdf", width = 12, height = 9)
print(p_dot_CC)
dev.off()

pdf("GO/Classify_before_GO/barplot_ego_CC.pdf", width = 12, height = 9)
print(p_bar_CC)
dev.off()

pdf("GO/Classify_before_GO/dotplot_ego_MF.pdf", width = 12, height = 9)
print(p_dot_MF)
dev.off()

pdf("GO/Classify_before_GO/barplot_ego_MF.pdf", width = 12, height = 9)
print(p_bar_MF)
dev.off()

# 6. 保存富集分析结果为 CSV 文件
ego_BP_df <- as.data.frame(ego_BP)
write.csv(ego_BP_df, "GO/Classify_before_GO/ego_BP.csv", row.names = FALSE)

ego_CC_df <- as.data.frame(ego_CC)
write.csv(ego_CC_df, "GO/Classify_before_GO/ego_CC.csv", row.names = FALSE)

ego_MF_df <- as.data.frame(ego_MF)
write.csv(ego_MF_df, "GO/Classify_before_GO/ego_MF.csv", row.names = FALSE)


###-----------------------------------
# 把整体富集分析的结果分类展示
###-----------------------------------


# 1. 将整体富集分析结果转换为数据框
ego_df <- as.data.frame(ego)

# 2. 利用 GO.db 获取每个 GO term 的 ONTOLOGY 信息
go_ont <- AnnotationDbi::select(GO.db,
                                keys = ego_df$ID,
                                columns = "ONTOLOGY",
                                keytype = "GOID")
# 注意：返回的 go_ont 可能包含重复项或顺序不一，下面按ID合并到结果中

# 3. 将ontology信息合并到整体富集结果中
ego_df_annotated <- merge(ego_df, go_ont, by.x = "ID", by.y = "GOID", all.x = TRUE)

# 4. 根据ONTOLOGY对整体结果进行拆分
ego_BP <- subset(ego_df_annotated, ONTOLOGY == "BP")
ego_CC <- subset(ego_df_annotated, ONTOLOGY == "CC")
ego_MF <- subset(ego_df_annotated, ONTOLOGY == "MF")

# （可选）将拆分后的结果保存为 CSV 文件
write.csv(ego_BP, "GO/Classify_after_GO/ego_overall_BP.csv", row.names = FALSE)
write.csv(ego_CC, "GO/Classify_after_GO/ego_overall_CC.csv", row.names = FALSE)
write.csv(ego_MF, "GO/Classify_after_GO/ego_overall_MF.csv", row.names = FALSE)

# 5. 可视化：一种方法是利用 clusterProfiler 的 dotplot/barplot 后用 facet_wrap 分面展示
# 将ontology信息添加回原来的 enrichResult 对象，便于 dotplot/barplot 绘图时使用
ego@result <- ego_df_annotated

# 使用 dotplot 进行展示，并根据ONTOLOGY分面
p_dot_overall <- dotplot(ego, showCategory = 20) +
  ggplot2::facet_wrap(~ONTOLOGY, scales = "free_y")
p_bar_overall <- barplot(ego, showCategory = 20) +
  ggplot2::facet_wrap(~ONTOLOGY, scales = "free_y")

# 显示图形
p_dot_overall
p_bar_overall

# 6. 保存图形为 PDF 文件
pdf("GO/Classify_after_GO/dotplot_ego_overall.pdf", width = 12, height = 9)
print(p_dot_overall)
dev.off()

pdf("GO/Classify_after_GO/barplot_ego_overall.pdf", width = 12, height = 9)
print(p_bar_overall)
dev.off()
