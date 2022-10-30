# Data integration
library(dplyr)
library(Seurat)
library(harmony)
library(viridis)
sample <- c("CTR1","CTR2","CTR3","CTR4","PBC1","PBC2","PBC3","PBC4","PBC5")
group <- c("CTR","CTR","CTR","CTR","PBC","PBC","PBC","PBC","PBC")
#Input: exp_list
for(i in 1:length(exp_list)){
    if(i == 1){
        input <- Read10X(data.dir = exp_list[i])
        data <- CreateSeuratObject(counts = input, project = sample[i], min.cells = 5)
        data <- RenameCells(data,add.cell.id = sample[i])
        data@meta.data$sample <- sample[i]
        data@meta.data$group <- group[i]
        data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
        ngenemax <- unname(quantile(data$nFeature_RNA,0.95))
        umimax <- unname(quantile(data$nCount_RNA,0.95))
        data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < ngenemax & percent.mt < 50 & nCount_RNA > 0 &  nCount_RNA  < umimax)
        data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
    }else{
        input <- Read10X(data.dir = exp_list[i])
        input <- Read10X(data.dir = exp_list[i])
        tmp <- CreateSeuratObject(counts = input, project = sample[i], min.cells = 5)
        tmp <- RenameCells(tmp,add.cell.id = sample[i])
        tmp@meta.data$sample <- sample[i]
        tmp@meta.data$group <- group[i]
        tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^MT-")
        ngenemax <- unname(quantile(tmp$nFeature_RNA,0.95))
        umimax <- unname(quantile(tmp$nCount_RNA,0.95))
        tmp <- subset(tmp, subset = nFeature_RNA > 200 & nFeature_RNA < ngenemax & percent.mt < 50 & nCount_RNA > 0 &  nCount_RNA  < umimax)
        tmp <- NormalizeData(tmp, normalization.method = "LogNormalize", scale.factor = 10000)
        data <- merge(data,tmp)
    }
    data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(data)
    data <- ScaleData(data, features = all.genes)
    data <- RunPCA(data, features = VariableFeatures(object = data))
    DimPlot(data, reduction = "pca")
    data <- JackStraw(data, num.replicate = 100)
    data <- ScoreJackStraw(data, dims = 1:20)
    plot1 <- JackStrawPlot(data, dims = 1:20)
    plot2<-ElbowPlot(data)
    CombinePlots(plots = list(plot1, plot2))
    data <- RunHarmony(data,"sample", plot_convergence = FALSE)
    data <- FindNeighbors(data, reduction = "harmony",dims = 1:20)
    data <- FindClusters(data, resolution = 1.2)
    data <- RunUMAP(data, reduction = "harmony",dims = 1:20)
    data <- RunTSNE(data, reduction = "harmony",dims = 1:20)
    pdf('./umap.cluster.pdf')
    DimPlot(data, reduction = "umap",label = F)
    dev.off()
	saveRDS(data,"./integrated_cluster.rds")
    for (l in levels (data)) {
        cluster1.markers <- FindMarkers(data, ident.1 = l, min.pct = 0.01,logfc.threshold = 0.25)
        cluster1.markers <- data.frame (gene = rownames (cluster1.markers),cluster1.markers)
        write.table (cluster1.markers,paste0("./cluster.",l,".diffgenes.xls"),row.names = F,col.names = T, sep='\t',quote = F)
    }
}

# Annotation:
#Input: celltypes
new.cluster.ids <- celltypes
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
saveRDS(data,"./integrated_add_celllabel.rds")
pdf('./umap.label.pdf')
DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

## Visualization
#Input: features
FeaturePlot(data, features = features,cols = c("grey","red"))
DoHeatmap(data, features = features,label = F) + scale_fill_viridis()
DotPlot(data, features = features,cols = c("blue","red"))

## Cell-cell interaction

#Input: inputdata
library(dplyr)
library(Seurat)
rds  <- readRDS (inputdata)
write.table(as.matrix(rds@assays$RNA@data),'cellphonedb_count.txt',sep='\t',quote = F)
rds@meta.data$celltype <- rds@active.ident
meta.data <- data.frame (rownames(rds@meta.data),rds@meta.data[,'celltype',drop = F]) %>% as.matrix()
write.table(meta.data,"cellphonedb_meta.txt",sep = "\t",quote = F,row.names = F)
### Shell
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --counts-data=gene_name
cellphonedb plot dot_plot
cellphonedb plot heatmap_plot cellphonedb_meta.txt

## Enrichment analysis.
library(DOSE)
library(org.Hs.eg.db)
library(dplyr)
library(clusterProfiler)

#Input: diffgene,GO_term,KEGG_term

genes <- read.table (diffgene,header = T)
data = bitr(genes$gene,
    fromType="SYMBOL",
    toType="ENTREZID",
    OrgDb="org.Hs.eg.db")


ggo <- groupGO(gene = data$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC",level = 3,readable = TRUE)
ego_ALL <- enrichGO(gene = data$ENTREZID, 
                OrgDb = org.Hs.eg.db,
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 1,
                qvalueCutoff = 1)

write.table(ego_ALL@result,"Enrichment_GO.xls",sep = "\t",col.names = T, row.names = F,quote = F)
ego_ALL <- filter(ego_ALL, ego_ALL@result$Description %in% GO_term)
pdf("Enrichment_GO.pdf")
dotplot(ego_ALL,title="Enrichment_GO_dot") + facet_grid(ONTOLOGY~.,scale='free', space='free_y')
dev.off()

kk <- enrichKEGG(gene = data$ENTREZID,
                 organism = 'hsa',
                 pvalueCutoff = 1)
write.table(kk@result,"Enrichment_KEGG.xls",sep = "\t",col.names = T, row.names = F,quote = F)
kk <- filter(kk, kk@result$Description %in% KEGG_term)
pdf("Enrichment_KEGG.pdf")
dotplot(kk,title="Enrichment_KEGG_dot")
dev.off()


