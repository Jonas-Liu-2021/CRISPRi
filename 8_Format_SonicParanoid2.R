library(dplyr)

# 读取CSV文件，根据实际情况调整分隔符（例如sep="\t"如果是制表符）
df <- read.csv("Data/AeromonasSSU-EcoliMG1655.csv", stringsAsFactors = FALSE)

# 只取最后两列（假设列名为 OrthoA 和 OrthoB）
df_subset <- df %>% select(OrthoA, OrthoB)

# 定义函数：从单元格文本中提取出得分为1的基因
# 解析思路：单元格内数据形如 "GeneID score GeneID score ..."（各项以空格分隔）
# 如果仅有一个基因得分为1，则返回该基因，否则返回 NA
extract_single_gene <- function(cell_text) {
  if (is.na(cell_text) || cell_text == "") return(NA)
  # 按空格分割
  tokens <- unlist(strsplit(cell_text, "\\s+"))
  # 至少需要两个 token 才能组成一对
  if (length(tokens) < 2) return(NA)
  
  # 每两个 token 为一组：第一个为基因ID，第二个为得分
  gene_ids <- tokens[seq(1, length(tokens), by = 2)]
  scores   <- tokens[seq(2, length(tokens), by = 2)]
  
  # 筛选得分为 "1" 的基因ID
  valid_genes <- gene_ids[scores == "1"]
  
  # 仅当正好只有1个符合条件时返回该基因
  if (length(valid_genes) == 1) {
    return(valid_genes)
  } else {
    return(NA)
  }
}

# 对 OrthoA 和 OrthoB 分别提取得分为1的基因
df_subset$Aeromonas <- sapply(df_subset$OrthoA, extract_single_gene)
df_subset$Ecoli     <- sapply(df_subset$OrthoB, extract_single_gene)

# 筛选出两侧都提取出基因的行，即1对1匹配的基因对
result <- df_subset %>% 
  filter(!is.na(Aeromonas) & !is.na(Ecoli)) %>% 
  select(Aeromonas, Ecoli)

# 查看结果
print(result)

# 如需保存结果到文件
write.csv(result, "Homologues/Homologues_genes_SonicParanoid2.csv", row.names = FALSE)
