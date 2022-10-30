#! /usr/bin/env Rscript
library(Seurat)
library(dplyr)
library(ggplot2)
#prepare information
#clusterrds: Integrated Seruat object rds before annotation.
#marker: Marker gene list with the first colum was celltype and the second cloumn was the marker genes corresponding to celltype.
#cluster: A table whit the first column was celltype, and the second column was cluster number corresponding to celltype.
#scale.max: 1
#scale.min: 0
#outdir: Output path of result.
#prefix: Prefix of output result,

print("rds loading...")
print(clusterrds)
PRO.cluster <- readRDS(clusterrds)
print("rds loaded.")
marker_list <- read.table(marker,header = T,sep = "\t",quote = "")
marker_tmp <- as.character(marker_list[,2])
names(marker_tmp) <- marker_list[,1]
marker_list <- marker_tmp
print("marker gene list: ")
print(marker_list)
cluster_list <- read.table(cluster,header = T,sep = "\t",quote = "")
cluster_tmp <- cluster_list[,2]
cluster_to_subset <- cluster_tmp
names(cluster_tmp) <- cluster_list[,1]
cluster_list <- cluster_tmp
print("relationship between cluster and cell types:")
print(cluster_list)
cluster_to_subset <- paste(cluster_to_subset,collapse = ",")
cluster_to_subset <- as.numeric(unique(unlist(strsplit(cluster_to_subset,split = ","))))
print("cluster_to_plot:")
print(cluster_to_subset)

#prepaer information use in following precess
print("start processing data.")
PRO.cluster <- subset(PRO.cluster,idents = cluster_to_subset)
data <- PRO.cluster@assays$RNA@data
cluster_info <- data.frame(PRO.cluster@active.ident)
names(cluster_info) <- "cluster"
print("cluster:")
print(head(cluster_info))
#define function
percentage <- function(x){
	whole_number <- length(x)
	number <- length(x[x>0])
	percent <- number/whole_number
	return(percent)
}
#set list to store signature_score of different cell types
signature_score <- list()
#####strat calculating scores and plotting
ct <- names(marker_list)
dir.create(paste(outdir,"/oridata_celltype",sep = ""))
for(i in 1:length(ct)){
	gene_tmp <- unlist(strsplit(marker_list[ct[i]],split = ","))
	print(gene_tmp)
	data_tmp <- data[gene_tmp,]
	data_tmp_sig_score <- data.frame(apply(data_tmp,2,function(x) prod(x)^(1/length(x))))
	names(data_tmp_sig_score) <- "sig_score_no_scale"
	if(max(data_tmp_sig_score$sig_score_no_scale) != min(data_tmp_sig_score$sig_score_no_scale)){
		k <- (scale.max-scale.min)/(max(data_tmp_sig_score$sig_score_no_scale)-min(data_tmp_sig_score$sig_score_no_scale))
	  }else if(max(data_tmp_sig_score$sig_score_no_scale) == min(data_tmp_sig_score$sig_score_no_scale)&max(data_tmp_sig_score$sig_score_no_scale) == 0){
		k <- 1
	  }
	print(max(data_tmp_sig_score$sig_score_no_scale))
	print(min(data_tmp_sig_score$sig_score_no_scale))
	print("k:")
	print(k)
	data_tmp_sig_score$sig_score_scale <- scale.min+k*(data_tmp_sig_score$sig_score_no_scale-min(data_tmp_sig_score$sig_score_no_scale))
	data_tmp_sig_score <- merge(data_tmp_sig_score,cluster_info,by = 0)
	print("data:")
	print(head(data_tmp_sig_score))
	write.table(data_tmp_sig_score,paste(outdir,"/oridata_celltype/",prefix,"_",ct[i],"_cell_data.xls",sep = ""),sep = "\t",quote = F,row.names = F)
	store <- data_tmp_sig_score %>% group_by(cluster) %>% summarise(mean=mean(sig_score_scale),percent = percentage(sig_score_scale))
	store <- data.frame(store)
	print(ct[i])
	print(head(store))
	store$celltype <- ct[i]
	signature_score[[ct[i]]] <- store
}
data_plot <- signature_score[[ct[1]]]
print(class(data_plot))
for(i in 2:length(names(signature_score))){
	print(names(signature_score)[i])
	data_plot <- rbind(data_plot,signature_score[[ct[i]]])
}
names(data_plot)[1] <- "cluster"
print("data_plot: ")
print(head(data_plot))
cluster_levels <- c()
celltype_levels <- c()
for(i in names(cluster_list)){
	print(i)
	celltype_tmp <- i
	celltype_levels <- c(celltype_levels,celltype_tmp)
	print(cluster_list[i])
	print(class(cluster_list[i]))
	cluster_tmp <- unlist(strsplit(as.character(cluster_list[i]),split = ","))
	cluster_levels <- c(cluster_levels,cluster_tmp)
}
print("cluster:")
print(cluster_levels)
cluster_levels <- unique(cluster_levels)
print("cluster without repeat:")
print(cluster_levels)
print("celltype:")
print(celltype_levels)
data_plot$cluster <- factor(data_plot$cluster,levels = cluster_levels)
data_plot$celltype <- factor(data_plot$celltype,levels = celltype_levels)
print("data_to_plot:")
print(names(data_plot))
p <- ggplot(data_plot,aes(x = cluster,y = celltype))+geom_point(aes(size = percent,color = mean))+scale_color_gradient(low="blue",high="red")+theme(panel.background = element_rect(color = "white",fill = "white"),axis.line = element_line(size = 0.5))+xlab("Cluster")+ylab("Lineage signature")+theme(legend.key = element_blank())+labs(size = "Fraction of cells",color = "Sig-score")
path <- paste(outdir,"/",prefix,"_sig_score_dotplot.pdf",sep = "")
pdf(path,width = 10)
print(p)
dev.off()
write.table(data_plot,paste(outdir,"/",prefix,"_data.xls",sep = ""),sep = "\t",quote = F)

	


