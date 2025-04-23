library(dplyr)
library(readr)
library(GenomicRanges)

# 读取RDS数据（GRanges数据）
gr=readRDS("Output_RDS/GRanges_log2FC_Operon.rds")

# 添加fitness_defect列，条件为log2FC<-2且padj<0.05
gr$fitness_defect <- gr$logFC<(-2) & gr$padj <0.05

# 按照Operon分组处理，标记essential、uncertain、nonessential和nonessential/internal promoter
# 逻辑顺序说明：
# 从操纵子的最后一个基因往前逐一检查：
# 1. 如果基因fitness_defect为FALSE，标记为 "nonessential"。
# 2. 如果基因fitness_defect为TRUE，第一个遇到的TRUE标记为 "essential"；
#    继续往前，如果仍为TRUE，则标记为 "uncertain"，直到遇到第一个FALSE；
#    遇到的第一个FALSE标记为 "nonessential/internal promoter"；
#    该位置前的所有基因，不论fitness_defect状态如何，均标记为 "uncertain/others"，
#    表示操纵子内部可能存在复杂的未知调控机制。
# *. 无论何种情况，一旦遇到NA，该多顺反子内当前基因和其前面的基因全部标注为uncertain
operon_df=as.data.frame(gr)

operon_df <- operon_df %>%
  group_by(Operon) %>%
  mutate(
    essential_status = {
      n <- n()
      status <- rep("", n)
      i <- n
      while(i >= 1) {
        # 如果当前基因的数据为 NA，则将当前位置及之前所有基因标记为 uncertain，并退出循环
        if (is.na(fitness_defect[i])) {
          status[1:i] <- "uncertain"
          break
        }
        if (!fitness_defect[i]) {
          status[i] <- "nonessential"
          i <- i - 1
        } else {
          status[i] <- "essential"
          i <- i - 1
          while(i >= 1) {
            # 同样检查 NA
            if (is.na(fitness_defect[i])) {
              status[1:i] <- "uncertain"
              i <- 0  # 退出内部循环
              break
            }
            if (fitness_defect[i]) {
              status[i] <- "uncertain"
              i <- i - 1
            } else {
              break
            }
          }
          if(i >= 1) {
            # 此处保证 fitness_defect[i] 为 FALSE（且非 NA，因为前面已检查）
            status[i] <- "nonessential/internal promoter"
            i <- i - 1
            while(i >= 1) {
              if (is.na(fitness_defect[i])) {
                status[1:i] <- "uncertain"
                i <- 0
                break
              }
              status[i] <- "uncertain/others"
              i <- i - 1
            }
          }
        }
      }
      status
    }
  )
# 多顺反子基因信息更新到GRanges对象
match_idx <- match(mcols(gr)$gene_id, operon_df$gene_id)
mcols(gr)$essential_status[!is.na(match_idx)] <- operon_df$essential_status[match_idx[!is.na(match_idx)]]


# 修改essential_status_2列
gr$essential_status_2=gr$essential_status
gr$essential_status_2[gr$essential_status_2=="uncertain/others"] = "uncertain"

# 保存csv数据
gr_df=as.data.frame(gr)
gr_df <- gr_df[, c("seqnames", 
                   "gene_id",
                   "gene_name",
                   "Operon",
                   "start", 
                   "end", 
                   "width", 
                   "strand",
                   "Type",
                   "logFC", 
                   "padj",
                   "fc_lt_2",
                   "fitness_defect",
                   "essential_status",
                   "essential_status_2",
                   "Function",
                   "GOs",
                   "KEGG_Pathway",
                   "COG_category")] # 调整列的顺序

write.csv(gr_df, file = "Data/Aeromonas_Integrated_Data.csv", row.names = FALSE)

# 保存RDS数据
saveRDS(gr, "Output_RDS/GRanges_Integrated_Data_2.rds")

