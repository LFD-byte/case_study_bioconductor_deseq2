# Run from the project folder: Rscript scripts/deseq2_demo.r
library(DESeq2)
library(airway)

# 加载示例数据
data(airway)
counts  <- assay(airway)                  # 原始 count 矩阵：行是基因，列是样本
coldata <- as.data.frame(colData(airway)) # 样本信息

# 设置对照组，否则默认按字母顺序选择
coldata$dex <- relevel(factor(coldata$dex), ref = "untrt")
coldata$cell <- factor(coldata$cell)

# 构建对象：cell 用于校正细胞系差异，dex 是要检验的因素（放在最后）
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ cell + dex
)

# 过滤低表达基因：至少 4 个样本的 count >= 10（4 为最小组的样本数）
keep <- rowSums(counts(dds) >= 10) >= 4
dds  <- dds[keep, ]

# 运行差异分析
dds <- DESeq(dds)

# 提取结果：处理组 vs 对照组
res <- results(dds, contrast = c("dex", "trt", "untrt"), alpha = 0.05)
summary(res)

# 筛选显著差异基因并保存
res_df <- as.data.frame(res[order(res$padj), ])
sig <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
nrow(sig)
head(sig)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
write.csv(res_df, "results/tables/deseq2_all_results.csv")
write.csv(sig, "results/tables/deseq2_significant.csv")

# 所有图保存到同一个 PDF
pdf("results/figures/deseq2_plots.pdf")

# MA 图
plotMA(res, ylim = c(-5, 5))

# PCA 图：检查样本分组和批次
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = c("dex", "cell"))

# 火山图
plot_df <- res_df[!is.na(res_df$padj), ]
sig_gene <- plot_df$padj < 0.05 & abs(plot_df$log2FoldChange) > 1
plot(plot_df$log2FoldChange, -log10(plot_df$padj),
     pch = 20, cex = 0.5,
     col = ifelse(sig_gene, "red", "grey"),
     xlab = "log2 fold change", ylab = "-log10 adjusted p-value")
abline(v = c(-1, 1), h = -log10(0.05), lty = 2)

dev.off()
