library(DESeq2)

# 读取数据
countData <- read.csv("Data/set1.csv", row.names = 1)
colData <- data.frame(
  row.names = colnames(countData),
  condition = factor(c("treatment","treatment", "treatment","control", "control", "control"))
) # 根据实验设计选出实验组和对照组

# 建立DESeq对象
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData   = colData,
  design    = ~ condition
)

dds <- DESeq(dds)

# 输出结果
res <- results(dds, contrast = c("condition", "treatment", "control"))
res_filtered <- res[!is.na(res$log2FoldChange), ] #去掉带NA的行

# 保存数据
write.csv(as.data.frame(res), file = "Data/DESeq2_results.csv")
write.csv(as.data.frame(res_filtered), file = "Data/DESeq2_results_filtered.csv")
