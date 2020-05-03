library(Seurat)
library(dplyr)
install.packages('patchwork')
library(patchwork)
library(reticulate)
Install.packages('umap')
library(umap)
library(virtualenv)
reticulate::py_install(packages ='umap-learn')

val <- load("/projectnb/bf528/users/group6/project_4/output/processed_panc_with_clusters.rda")
#cells <- readRDS("/projectnb/bf528/project_4_scrnaseq/GSM2230760_seurat.rda")

pbmc.markers <- FindAllMarkers(panc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10 = pbmc.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
#write.csv(top10,"DE_genes.csv")


database <- read.table('PanglaoDB_markers_27_Mar_2020.tsv', sep = '\t', header=TRUE, stringsAsFactors = FALSE)
pbmc.markers$cell_type <- NA
pbmc.markers$cannonical <- NA
pbmc.markers$ubiquitous_index <- NA
pbmc.markers$specificity <- NA
match <- which(database$official.gene.symbol == gene)
for (row in 1:nrow(pbmc.markers)) {
  gene <- pbmc.markers[row, "gene"]
  if (gene %in% database$official.gene.symbol)
    pbmc.markers$cell_type[row] <- database[match[1], 3]
}


new_cluster=c('Acinar', 'Adipocyte','Acinar','Mast', 'Muller cells', 'Schwann', 'Acinar','Germ cells','Germ cells','Beta', 'Not found','Delta','Dendritic', 'Mast','Alpha cells','Not found')
names(new_cluster)=levels(panc)
panc=RenameIdents(panc,new_cluster)
plot_umap <- RunUMAP(panc, dims = 1:10)
plot_umap
DimPlot(plot_umap, reduction="umap", pt.size=0.5)
                                                         

DoHeatmap(panc, features = top10$gene) 



