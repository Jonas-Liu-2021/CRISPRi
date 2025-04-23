library(clusterProfiler)
library(GenomicRanges)
library(KEGGREST)
library(AnnotationDbi)


# 读取RDS文件的整合数据
gr=readRDS("Output_RDS/GRanges_Integrated_Data_2.rds")


# --------------------------
# 构建 KEGG 映射表
# --------------------------
# 提取基因ID和对应的 KEGG 注释
gene_ids <- gr$gene_id      # 基因ID向量
kegg_terms <- gr$KEGG_Pathway     # KEGG注释（可能为逗号分隔的字符串）

# 对每个基因拆分 KEGG 编号，生成映射列表
gene2kegg_list <- lapply(kegg_terms, function(x) {
  if (is.na(x)) {
    character(0)
  } else {
    unlist(strsplit(x, split = ","))
  }
})
gene_counts_kegg <- sapply(gene2kegg_list, length)

# 构建映射表，只保留有KEGG注释的基因
gene2kegg_df <- data.frame(
  KEGG = unlist(gene2kegg_list),
  gene = rep(gene_ids, gene_counts_kegg)
)

# -----------------------------
# 构建 TERM2NAME 映射表：KEGG编号 -> KEGG描述
# -----------------------------
unique_kegg_ids <- unique(gene2kegg_df$KEGG)

# -----------------------------
# 使用parallel进行并行计算
# -----------------------------
# 定义查询函数：返回名称或 NA
# query_kegg <- function(id) {
#   tryCatch({
#     res <- KEGGREST::keggGet(id)
#     if (length(res) > 0 && !is.null(res[[1]]$NAME)) {
#       return(res[[1]]$NAME)
#     } else {
#       return(NA)
#     }
#   }, error = function(e) {
#     return(NA)
#   })
# }
# 
# # 初始化结果向量，所有编号初始设为 NA
# results <- setNames(rep(NA, length(unique_kegg_ids)), unique_kegg_ids)
# remaining_ids <- unique_kegg_ids
# 
# max_attempts <- 7  # 最大尝试次数
# attempt <- 1
# 
# # 创建并行集群（使用可用核心数减 1）
# numCores <- detectCores() - 2
# cl <- makeCluster(numCores)
# 
# while(length(remaining_ids) > 0 && attempt <= max_attempts) {
#   cat(sprintf("尝试次数 %d: 还有 %d 个编号未成功查询。\n", attempt, length(remaining_ids)))
#   
#   # 使用 pblapply 带进度条进行并行查询
#   temp_results <- pblapply(remaining_ids, query_kegg, cl = cl)
#   names(temp_results) <- remaining_ids
#   
#   # 更新 results 向量
#   for(id in remaining_ids) {
#     results[id] <- temp_results[[id]]
#   }
#   
#   # 筛选仍然未获取到结果的编号
#   remaining_ids <- names(results)[is.na(results)]
#   attempt <- attempt + 1
#   
#   if(length(remaining_ids) > 0) {
#     cat("等待 10 秒后重新尝试...\n")
#     Sys.sleep(10)
#   }
# }
# 
# # 关闭并行集群
# stopCluster(cl)
# 
# # 保存结果到文件
# saveRDS(results, file = "Output_RDS/kegg_TERMS2NAME.rds")
# 
# # 输出未能获取到结果的 KEGG 编号
# failed_ids <- names(results)[is.na(results)]
# if(length(failed_ids) > 0) {
#   cat("以下编号未能获取到结果：\n")
#   print(failed_ids)
# } else {
#   cat("所有编号均已成功获取到结果！\n")
# }

# results=readRDS("Output_RDS/kegg_TERMS2NAME.rds")



# 获取 KEGG 通路信息
map_info <- keggList("pathway")            # 获取参考通路 (mapXXXX)
ko_info  <- keggList("pathway", "ko")       # 获取 KO 通路 (koXXXX)

# 去除前缀 "path:"，构建映射表
map2name <- setNames(map_info, sub("path:", "", names(map_info)))
ko2name  <- setNames(ko_info,  sub("path:", "", names(ko_info)))

# 定义函数，根据 ID 前缀返回对应名称
getTerm2Name <- function(ids) {
  sapply(ids, function(id) {
    if (startsWith(id, "map")) {
      map2name[id]
    } else if (startsWith(id, "ko")) {
      ko2name[id]
    } else {
      NA
    }
  })
}

# 对 unique_kegg_ids 进行映射
term2name <- getTerm2Name(unique_kegg_ids)
term2name_df <- data.frame(TERM = unique_kegg_ids, 
                           NAME = term2name, 
                           stringsAsFactors = FALSE)
term2name_df




# -----------------------------
# 定义目标基因集和背景基因集（与GO分析相同）
# -----------------------------
essential_genes <- gr$gene_id[gr$essential_status == "essential"]
background_genes <- gr$gene_id

# -----------------------------
# KEGG 富集分析
# -----------------------------
ekegg <- enricher(gene          = essential_genes,
                  universe      = background_genes,
                  TERM2GENE     = gene2kegg_df,
                  TERM2NAME     = term2name_df,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.9,
                  qvalueCutoff  = 0.9)

# 查看富集结果
head(ekegg)

# 绘图展示富集结果
dotplot(ekegg, showCategory = 20)
barplot(ekegg, showCategory = 20)

p_dot=dotplot(ekegg, showCategory = 20) 
p_bar=barplot(ekegg, showCategory = 20)

pdf("KEGG/dotplot_ekegg.pdf", width = 8, height = 6)
print(p_dot)
dev.off()

pdf("KEGG/barplot_ekegg.pdf", width = 8, height = 6)
print(p_bar)
dev.off()

# 保存为CSV文件
ekegg_df=as.data.frame(ekegg)
write.csv(ekegg_df,"KEGG/ekegg.csv")
