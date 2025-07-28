#Author: Ziyan Lin
#r/4.0.3
###############################
#load package
library(devtools)
library(ashr)
library(dplyr)
library(openxlsx)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(glue)
library(magrittr)
library(RColorBrewer)
library(plotly)
library(pheatmap)
library(stringr)

###############################
#input data
sampleTable <- read.csv("sampleTable.csv",header=T)
#sampleTable
#  SampleName  fileName group
#1  C_2_1 C_2_1.txt        C2
#2  C_2_2 C_2_2.txt        C2

Cond1 <- "STR-M-DSVA" #control
Cond2 <- "STR-M-XDP"
sampleTable$condition <- sampleTable$group
sampleTable$condition <- factor(sampleTable$condition,levels=c(Cond1,Cond2)) #2/1 #control group as the first condition
sampleTable

dds <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory="./rawcounts/",design= ~condition)
dds <- DESeq(dds)

###############################
#VST
vsd = varianceStabilizingTransformation(dds, blind = TRUE)

###############################
#PCA
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
saveRDS(pcaData,"pcaData.rds")

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=9) +
  geom_text(aes(label=name), vjust=0, nudge_y=1, size=2) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic() + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add border around the plot area
  )


###############################

###############################
#distance plot
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, names(vsd$sizeFactor), sep="-")
colnames(sampleDistMatrix) <- paste(vsd$condition, names(vsd$sizeFactor), sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

###############################
#Normalized counts
norm_counts_table = counts(dds, normalized = TRUE) %>% round(2) %>% as_tibble(rownames = "gene")
write.csv(norm_counts_table, "cts_norm.csv")

###############################
#get DEGs
resultsNames(dds)

###############################
group_name <- paste0(gsub("-",".",Cond2),"_vs_",gsub("-",".",Cond1))
group_name

###############################
res_unshrunk <- results(dds, name=paste0("condition_",group_name),cooksCutoff = F)

###############################
#get shirnk log fold change
res <- lfcShrink(dds,res_unshrunk,coef=paste0("condition_",group_name),type = "ashr")
res = res[order(res$padj, res$pvalue, -res$baseMean), ]# sort results so most significant are first

###############################
#save unmodified results object
saveRDS(res, "desep2-res.rds")
#save the unmodified results table as csv
res_tbl = as_tibble(res, rownames = "gene") %>% dplyr::arrange(padj, pvalue, -baseMean)
#write.csv(res_tbl,paste0(group_name,".csv"))

###############################
#add unshrunk fold change to results
res_unshrunk_tbl = as_tibble(res_unshrunk, rownames = "gene")
res_unshrunk_tbl = dplyr::select(res_unshrunk_tbl, gene, log2FCunshrunk = log2FoldChange)
res_tbl = left_join(res_tbl, res_unshrunk_tbl, by = "gene")

res_sum <- full_join(norm_counts_table,res_tbl,by="gene")
head(res_sum)
colnames(res_sum)
###############################
group_name <- gsub("\\.","-",group_name) #format file name

#format results for excel export
res_tbl =
res_sum %>%
dplyr::mutate(
baseMean       = round(baseMean, 1),
log2FC         = round(log2FoldChange, 2),
log2FCunshrunk = round(log2FCunshrunk, 2),
pvalue         = if_else(pvalue < 0.00001, pvalue, round(pvalue, 5)),
padj           = if_else(padj < 0.00001, padj, round(padj, 5))
) %>%
dplyr::select(gene,paste0(Cond1,"1"),paste0(Cond1,"2"),paste0(Cond1,"3"),paste0(Cond2,"1"),paste0(Cond2,"2"),paste0(Cond2,"3"),paste0(Cond2,"4"),baseMean, log2FC, log2FCunshrunk, pvalue, padj)

write.csv(res_tbl,paste0(group_name,".csv"))

###############################
#save DEGs
write.xlsx(res_tbl,paste0(group_name,".xlsx"))
write.xlsx(subset(res_tbl,padj<0.05),paste0(group_name,"-q005.xlsx"))
#write.xlsx(subset(res_tbl,padj<0.1),paste0(group_name,"-q01.xlsx"))
write.xlsx(subset(res_tbl,abs(log2FoldChange)>1 & padj<0.05),paste0(group_name,"-q005-2fold.xlsx"))

###############################
#volcano plot
##P-value, or adjusted p-value, be careful
res_tbl <- res_tbl %>%
  mutate(color = ifelse(log2FoldChange > 0 & padj < 0.05, "Treated",
                        ifelse(log2FoldChange < 0 & padj < 0.05, "Untreated", "none")))

dot.col = c("#F76F61", "#4D8AC9", "#636363")

my.VPlot <- ggplot(res_tbl, aes(x = log2FoldChange, y = -log10(padj), text = gene)) +
  geom_point(aes(color = factor(color)), size = 3, alpha = 0.5, na.rm = TRUE) +  # Add gene points
  geom_text_repel(
    data = head(res_tbl[order(res_tbl$padj),], 20),
    aes(label = gene),
    color = "black",
    size = 6,
    point.padding = 0.1,
    max.overlaps = Inf
  ) +
  theme_bw(base_size = 16) +  # Clean up theme
  theme(legend.position = "none") +  # Remove legend
  theme(
    legend.position = "none",  # Remove legend
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.ticks = element_line(color = "black"),  # Add axis ticks
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove panel background
    plot.background = element_blank()  # Remove plot background
  ) +
  ggtitle("Str XDP vs Dsva") +  # Add title
  xlab("log2FoldChange") +  # x-axis label
  ylab("-log10 (adjusted p-value)") +  # y-axis label
  geom_vline(xintercept = 0, colour = "black") +  # Add line at 0
  geom_hline(yintercept = 0, colour = "black") +  # p(0.05) = 1.3
  scale_color_manual(values = c("Treated" = dot.col[1],
                               "Untreated" = dot.col[2],
                                "none" = dot.col[3])) +
  scale_y_continuous(trans = "log1p")

print(my.VPlot)

# An interactive plot using plotly and display it
interactive_plot <- ggplotly(my.VPlot)
print(interactive_plot)

# If needed, save the interactive plot to an HTML file
OUT.PUT <- "volcano.html"
htmlwidgets::saveWidget(interactive_plot, OUT.PUT)

###############################
#Heatmaps
# VST - Perform variance stabilizing transformation
vsd = assay(varianceStabilizingTransformation(dds, blind = TRUE))

mat = vsd
samples_all = colnames(dds)

#row_subset <- readLines("genes.txt")

##################
#top 50/100 genes
#row_subset <- readLines("genes.txt")
n_top = 50
tbl <- res_tbl
tbl_sub <- tbl[order(tbl$padj),] %>% head(n_top)
row_subset <- tbl_sub$gene
##################

row_subset <- unique(row_subset)

row_subset = intersect(row_subset, rownames(mat))
col_subset  = samples_all

mat2 = mat[row_subset, col_subset]

cell_colors <- colorRampPalette(c("#0033A0", "#4D8AC9", "#F7F7F7", "#F76F61", "#C70039"))(50)

pheatmap(mat2, 
         color = cell_colors, 
         border_color = NA, 
         scale = "row", 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         main = "Heatmap", 
         fontsize_row = 9, 
         fontsize_col = 12, 
         show_rownames = TRUE,
         cellwidth=12,
         cellheight=10)
