# 加载必要的包
library(GenomicRanges)
library(ggbio)

# 读取数据
gr = readRDS("Output_RDS/GRanges_Integrated_Data_2.rds")

# 筛选 CDS 基因
gr_CDS <- gr[gr$Type == "CDS"]
gr_CDS <- GenomicRanges::reduce(gr_CDS)

# 根据链信息进一步分为正链和负链
gr_CDS_pos <- gr_CDS[strand(gr_CDS) == "+"]
gr_CDS_pos <- GenomicRanges::reduce(gr_CDS_pos)
gr_CDS_neg <- gr_CDS[strand(gr_CDS) == "-"]
gr_CDS_neg <- GenomicRanges::reduce(gr_CDS_neg)

## 绘制环形图(有正负链)
# 先筛选并合并essential genes，防止基因堆叠
gr_fc1 <- gr[gr$fc_lt_1]   # 864个essential gene（logFC < -1）
gr_fc1_merged <- GenomicRanges::reduce(gr_fc1)

gr_fc2 <- gr[gr$fc_lt_2]   # 630个essential gene (logFC < -2)
gr_fc2_merged <- GenomicRanges::reduce(gr_fc2)

# 根据 logFC 的绝对值范围构造对称的 ylim
ylim_range <- c(-max(abs(gr$logFC), na.rm = TRUE), max(abs(gr$logFC), na.rm = TRUE))

p <- ggbio() +
  circle(gr[!is.na(gr$logFC)], 
         geom = "bar",
         aes(y = logFC), 
         color = "gray70",
         radius = 11,      # 基线半径，0 对应半径 10
         trackWidth = 15,  # 环带的厚度
         space.skip = 0,
         ylim = ylim_range) +
  circle(gr_fc2_merged,
         geom = "rect",
         color = "green", 
         space.skip = 0,
         radius = 26,
         trackWidth = 2,
         fill = "green",
         rect.inter.n = 1800) +
  circle(gr_CDS_pos,        # 正链 CDS
         geom = "rect",
         color = "steelblue",
         space.skip = 0,
         radius = 29,       # 根据需要调整 radius，使正负链层不重叠
         trackWidth = 2,
         fill = "steelblue",
         rect.inter.n = 1800) +
  circle(gr_CDS_neg,        # 负链 CDS，可选用另一种颜色
         geom = "rect",
         color = "red",
         space.skip = 0,
         radius = 32,       # 注意选择合适的 radius 使两圈不重叠
         trackWidth = 2,
         fill = "red",
         rect.inter.n = 1800) +
  circle(gr,
         geom = "ideo",
         fill = "gray",
         radius = 35.05,
         trackWidth = 0.005,
         space.skip = 0,
         rect.inter.n = 1800) +
  circle(gr,
         geom = "scale",
         scale.n = 5,
         size = 5,
         radius = 35,
         trackWidth = 2,
         space.skip = 0) + # 必须用 space.skip = 0 才能让首尾闭合
  annotate("text",
           x = 0,
           y = 0,
           label = "Aeromonas \ndhakensis \ngenome",
           size = 9,
           fontface = "bold")
p

# 保存图片
pdf("Figures/cicular_pos_neg.pdf", width = 10, height = 9.5)
print(p)  # 显示图形
dev.off()

# 绘制环形图（考虑CDS但不考虑正负链）
p <- ggbio() +
  circle(gr[!is.na(gr$logFC)], 
         geom = "bar",
         aes(y = logFC), 
         color = "gray70",
         radius = 11,      # 基线半径，0 对应半径 10
         trackWidth = 15,  # 环带的厚度
         space.skip = 0,
         ylim = ylim_range) +
  circle(gr_fc2_merged,
         geom = "rect",
         color = "red", 
         space.skip = 0,
         radius = 26,
         trackWidth = 2,
         fill = "red",
         rect.inter.n = 1800) +
  circle(gr_CDS,        # CDS
         geom = "rect",
         color = "steelblue",
         space.skip = 0,
         radius = 29,       
         trackWidth = 2,
         fill = "steelblue",
         rect.inter.n = 1800) +
  
  circle(gr,
         geom = "ideo",
         fill = "gray",
         radius = 32.05,
         trackWidth = 0.005,
         space.skip = 0,
         rect.inter.n = 1800) +
  circle(gr,
         geom = "scale",
         scale.n = 5,
         size = 5,
         radius = 32,
         trackWidth = 2,
         space.skip = 0) + # 必须用 space.skip = 0 才能让首尾闭合
  annotate("text",
           x = 0,
           y = 0,
           label = "Aeromonas \ndhakensis \ngenome",
           size = 9,
           fontface = "bold")
p

# 保存图片
pdf("Figures/cicular_CDS.pdf", width = 10, height = 9.5)
print(p)  # 显示图形
dev.off()


# 绘制环形图（包含essential）
# 筛选 essential 基因
gr_ess <- gr[gr$essential_status_2 == "essential"]
gr_ess <- GenomicRanges::reduce(gr_ess)


p <- ggbio() +
  circle(gr[!is.na(gr$logFC)], 
         geom = "bar",
         aes(y = logFC), 
         color = "gray70",
         radius = 11,      # 基线半径，0 对应半径 10
         trackWidth = 15,  # 环带的厚度
         space.skip = 0,
         ylim = ylim_range) +
  circle(gr_ess,
         geom = "rect",
         color = "orange", 
         space.skip = 0,
         radius = 26,
         trackWidth = 2,
         fill = "orange",
         rect.inter.n = 1800) +
  circle(gr_fc2_merged,
         geom = "rect",
         color = "green", 
         space.skip = 0,
         radius = 29,
         trackWidth = 2,
         fill = "green",
         rect.inter.n = 1800) +
  circle(gr_CDS_pos,        # 正链 CDS
         geom = "rect",
         color = "steelblue",
         space.skip = 0,
         radius = 32,       # 根据需要调整 radius，使正负链层不重叠
         trackWidth = 2,
         fill = "steelblue",
         rect.inter.n = 1800) +
  circle(gr_CDS_neg,        # 负链 CDS，可选用另一种颜色
         geom = "rect",
         color = "red",
         space.skip = 0,
         radius = 35,       # 注意选择合适的 radius 使两圈不重叠
         trackWidth = 2,
         fill = "red",
         rect.inter.n = 1800) +
  circle(gr,
         geom = "ideo",
         fill = "gray",
         radius = 38.05,
         trackWidth = 0.005,
         space.skip = 0,
         rect.inter.n = 1800) +
  circle(gr,
         geom = "scale",
         scale.n = 5,
         size = 5,
         radius = 38,
         trackWidth = 2,
         space.skip = 0) + # 必须用 space.skip = 0 才能让首尾闭合
  annotate("text",
           x = 0,
           y = 0,
           label = "Aeromonas \ndhakensis \ngenome",
           size = 9,
           fontface = "bold")
p

# 保存图片
pdf("Figures/cicular_pos_neg_ess.pdf", width = 10, height = 9.5)
print(p)  # 显示图形
dev.off()


# 绘制环形图（包含essential） V 2.0
p <- ggbio() +
  circle(gr[!is.na(gr$logFC)], 
         geom = "bar",
         aes(y = logFC), 
         color = "mediumpurple",
         radius = 11,      # 基线半径，0 对应半径 10
         trackWidth = 15,  # 环带的厚度
         space.skip = 0,
         ylim = ylim_range) +
  circle(gr_ess,
         geom = "rect",
         color = "orange", 
         space.skip = 0,
         radius = 26,
         trackWidth = 2,
         fill = "orange",
         rect.inter.n = 1800) +
  circle(gr_fc2_merged,
         geom = "rect",
         color = "green", 
         space.skip = 0,
         radius = 29,
         trackWidth = 2,
         fill = "green",
         rect.inter.n = 1800) +
  circle(gr_CDS_pos,        # 正链 CDS
         geom = "rect",
         color = "steelblue",
         space.skip = 0,
         radius = 32,       # 根据需要调整 radius，使正负链层不重叠
         trackWidth = 2,
         fill = "steelblue",
         rect.inter.n = 1800) +
  circle(gr_CDS_neg,        # 负链 CDS，可选用另一种颜色
         geom = "rect",
         color = "red",
         space.skip = 0,
         radius = 35,       # 注意选择合适的 radius 使两圈不重叠
         trackWidth = 2,
         fill = "red",
         rect.inter.n = 1800) +
  circle(gr,
         geom = "ideo",
         fill = "black",
         radius = 38.05,
         trackWidth = 0.005,
         space.skip = 0,
         rect.inter.n = 1800) +
  circle(gr,
         geom = "scale",
         scale.n = 5,
         size = 500,
         radius = 38,
         trackWidth = 2,
         space.skip = 0) + # 必须用 space.skip = 0 才能让首尾闭合
  annotate("text",
           x = 0,
           y = 0,
           label = "Aeromonas \ndhakensis \nSSU genome",
           size = 900,
           fontface = "bold")
p
pdf("Figures/cicular_pos_neg_ess3.pdf", width = 1000, height = 950)
print(p)  # 显示图形
dev.off()

# 绘制环形图（包含essential） V 3.0
p <- ggbio() +
  circle(gr[!is.na(gr$logFC)], 
         geom = "bar",
         aes(y = logFC), 
         color = "gray70",
         radius = 11,      # 基线半径，0 对应半径 10
         trackWidth = 15,  # 环带的厚度
         space.skip = 0,
         ylim = ylim_range) +
  circle(gr_ess,
         geom = "rect",
         color = "steelblue", 
         space.skip = 0,
         radius = 26,
         trackWidth = 2,
         fill = "steelblue",
         rect.inter.n = 1800) +
  circle(gr_fc2_merged,
         geom = "rect",
         color = "red", 
         space.skip = 0,
         radius = 29,
         trackWidth = 2,
         fill = "red",
         rect.inter.n = 1800) +
  circle(gr_CDS_pos,        # 正链 CDS
         geom = "rect",
         color = "gray70",
         space.skip = 0,
         radius = 32,       # 根据需要调整 radius，使正负链层不重叠
         trackWidth = 2,
         fill = "gray70",
         rect.inter.n = 1800) +
  circle(gr_CDS_neg,        # 负链 CDS，可选用另一种颜色
         geom = "rect",
         color = "gray70",
         space.skip = 0,
         radius = 35,       # 注意选择合适的 radius 使两圈不重叠
         trackWidth = 2,
         fill = "gray70",
         rect.inter.n = 1800) +
  circle(gr,
         geom = "ideo",
         fill = "black",
         radius = 38.05,
         trackWidth = 0.005,
         space.skip = 0,
         rect.inter.n = 1800) +
  circle(gr,
         geom = "scale",
         scale.n = 5,
         size = 500,
         radius = 38,
         trackWidth = 2,
         space.skip = 0) + # 必须用 space.skip = 0 才能让首尾闭合
  annotate("text",
           x = 0,
           y = 0,
           label = "Aeromonas \ndhakensis \nSSU genome",
           size = 900,
           fontface = "bold")
p
pdf("Figures/cicular_pos_neg_ess4.pdf", width = 1000, height = 950)
print(p)  # 显示图形
dev.off()

