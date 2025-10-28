
library(pheatmap)
library(RColorBrewer)


rm(list=ls())
gc()

setwd("/home/luchun/scRNA/Bone_neu/03_age_ovx2/pyscenic/")

gene_len <- read.csv("Age_len_names.csv",header=TRUE)


regulonActivity_byGroup <- read.csv("Age_activity_average.csv",header=TRUE,row.names=1)

rownames(regulonActivity_byGroup) <- gene_len$new_name[match(rownames(regulonActivity_byGroup),gene_len$tfs)]



TFs_1 <- read.table("./Age_GLM/Neu_VS_Non_neu_TFs.txt",header=TRUE)
TFs <- gene_len$new_name[match(TFs_1$TFs,gene_len$tfs)]

TFs_1 <- read.table("./Age_GLM/Neu_Cluster_6_TFs.txt",header=TRUE)
TFs <- gene_len$new_name[match(TFs_1$TFs,gene_len$tfs)]

TFs_1 <- read.table("./Age_GLM/Neu_Cluster_TFs.txt",header=TRUE)
TFs <- gene_len$new_name[match(TFs_1$TFs,gene_len$tfs)]


regulons.use <- TFs



colnames(regulonActivity_byGroup)

regulonActivity_byGroup2 <- regulonActivity_byGroup[,c(2,3,6,4,7,1,5)]


result <- regulonActivity_byGroup2[regulons.use,]
head(result)

colnames(result) <- c('Elane+Prtn3+ progenitor', 'H2afz+Hmgb2+ proliferating', 'Ngp+Lcn2+ immature','Ifitm6+Ltf+ immature', 'Retnlg+Mmp8+ mature', 'Ccl6+Sell+ mature','Il1b+Srgn+ mature')




#画热图————————————————————————————————————
typelist = c('Elane+Prtn3+ progenitor', 'H2afz+Hmgb2+ proliferating', 'Ngp+Lcn2+ immature','Ifitm6+Ltf+ immature', 'Retnlg+Mmp8+ mature', 'Ccl6+Sell+ mature','Il1b+Srgn+ mature')

col_ann <- data.frame(preAnno=typelist)
rownames(col_ann) <- typelist



ann_colors <- list(
  preAnno = c(`Elane+Prtn3+ progenitor` = "#76afda",`H2afz+Hmgb2+ proliferating` = "#8264CC",`Ngp+Lcn2+ immature` = "#009E73",`Ifitm6+Ltf+ immature` = "#F0E442",`Retnlg+Mmp8+ mature`="#f06152",`Ccl6+Sell+ mature`="#0072B2",`Il1b+Srgn+ mature`="#fc9272"))

 
cols <- colorRampPalette(rev(brewer.pal(n=11,name="RdBu")[c(1:3,5:7,9:11)]))(100)



jpeg("02_Age_AS_data_heatmap.jpg",height=700)
pheatmap(result,scale = "row",
		cluster_cols=FALSE,cluster_rows=TRUE,
		show_colnames=FALSE,show_rownames=TRUE,
		annotation_col=col_ann,
		annotation_colors=ann_colors,
		fontsize_row=8,
		fontsize_col=10,
		cellheight=10,cellwidth=20,
		color=cols,
		treeheight_row=20,treeheight_col=20,
		border_color="white", 
		main="Regulon activity\n(neutrophil versus non-neutrophil)\nTop 20 sorted by t-value",
		angle_col ="45")
dev.off()



pdf("02_Age_AS_data_heatmap.pdf",height=10.21)
pheatmap(result,scale = "row",
		cluster_cols=FALSE,cluster_rows=TRUE,
		show_colnames=FALSE,show_rownames=TRUE,
		annotation_col=col_ann,
		annotation_colors=ann_colors,
		fontsize_row=8,
		fontsize_col=10,
		cellheight=10,cellwidth=20,
		color=cols,
		treeheight_row=20,treeheight_col=20,
		border_color="white", 
		main="Regulon activity\n(neutrophil versus non-neutrophil)\nTop 20 sorted by t-value",
		angle_col ="45")
dev.off()





jpeg("03_Age_AS_data_heatmap.jpg",width=550)
pheatmap(result,scale = "row",
		cluster_cols=FALSE,cluster_rows=TRUE,
		show_colnames=FALSE,show_rownames=TRUE,
		annotation_col=col_ann,
		annotation_colors=ann_colors,
		fontsize_row=8,
		fontsize_col=10,
		cellheight=10,cellwidth=20,
		color=cols,
		treeheight_row=20,treeheight_col=20,
		border_color="white", 
		main="Regulon activity\n(Ccl6+Sell+ mature versus other clusters)\nTop 20 sorted by t-value",
		angle_col ="45")
dev.off()



pdf("03_Age_AS_data_heatmap.pdf",width=8)
pheatmap(result,scale = "row",
		cluster_cols=FALSE,cluster_rows=TRUE,
		show_colnames=FALSE,show_rownames=TRUE,
		annotation_col=col_ann,
		annotation_colors=ann_colors,
		fontsize_row=8,
		fontsize_col=10,
		cellheight=10,cellwidth=20,
		color=cols,
		treeheight_row=20,treeheight_col=20,
		border_color="white", 
		main="Regulon activity\n(Ccl6+Sell+ mature versus other clusters)\nTop 20 sorted by t-value",
		angle_col ="45")
dev.off()




jpeg("04_Age_AS_data_heatmap.jpg",width=550,height=800)
pheatmap(result,scale = "row",
		cluster_cols=FALSE,cluster_rows=TRUE,
		show_colnames=FALSE,show_rownames=TRUE,
		annotation_col=col_ann,
		annotation_colors=ann_colors,
		fontsize_row=8,
		fontsize_col=10,
		cellheight=10,cellwidth=20,
		color=cols,
		treeheight_row=20,treeheight_col=20,
		border_color="white", 
		main="Regulon activity\n(Each neu cluster versus other neu clusters)\nTop 20 sorted by t-value",
		angle_col ="45")
dev.off()



pdf("04_Age_AS_data_heatmap.pdf",width=8,height=11.67)
pheatmap(result,scale = "row",
		cluster_cols=FALSE,cluster_rows=TRUE,
		show_colnames=FALSE,show_rownames=TRUE,
		annotation_col=col_ann,
		annotation_colors=ann_colors,
		fontsize_row=8,
		fontsize_col=10,
		cellheight=10,cellwidth=20,
		color=cols,
		treeheight_row=20,treeheight_col=20,
		border_color="white", 
		main="Regulon activity\n(Each neu cluster versus other neu clusters)\nTop 20 sorted by t-value",
		angle_col ="45")
dev.off()






#—————————————————————————————————— network  —————————————————————————
sce_SCENIC <- open_loom("./data/Age_total_SCENIC.loom")

regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
class(regulons)
ls(regulons)


net_dat <- data.frame()

for(i in c("OTX1(+)","NHLH1(+)",'EGR3(+)','NFIL3(+)','NFE2L2(+)')){
	to_lis <- regulons[[i]]	
	from_lis <- rep(substr(i, 1, nchar(i) - 3),length(to_lis))
	
	tmp <- data.frame(from=from_lis,to=to_lis,cols=from_lis)	
	net_dat <- rbind(net_dat,tmp)
}

names(which(table(net_dat$to) != 1))

index <- which(net_dat$to %in% names(which(table(net_dat$to) != 1)))


net_dat$cols[index] <- "Shared"

write.csv(net_dat[,c(1,2)],"net_dat.csv",quote=FALSE,row.names=FALSE)


meta <- net_dat[,c(2,3)]
colnames(meta) <- c("Gene","Group")
index <- which(meta$Gene %in% c("OTX1","NHLH1",'EGR3','NFIL3','NFE2L2'))
meta2 <- meta[-index,]

meta3 <- data.frame(Gene=c("OTX1","NHLH1",'EGR3','NFIL3','NFE2L2'),Group=rep("TF",5))

meta_dat <- rbind(meta2,meta3)

write.csv(meta_dat,"meta_dat.csv",quote=FALSE,row.names=FALSE)





#—————————————————————————————————— network  —————————————————————————
library(ggraph)

rm(list=ls())
gc()


grn <- read.table("/home/luchun/scRNA/Bone_neu/03_age_ovx2/pyscenic/data/Age_adjacencies.tsv",sep='\t',header=T,stringsAsFactors=F)
head(grn)



###OTX1————————————————————————————
ntop2=50

inregulons1=c("Nfib")
c1 <- which(grn$TF %in% inregulons1)
grn2 <- grn[c1,]

head(grn2,20)

for (tf in inregulons1) {
	
	tf_col <- "#d73027"
	target_col <- "#2873B3"
	motif_target_col <- "#e8743c"
	edge_col <- "#737373"
	
	tmp <- subset(grn2,TF==tf)
	
	if (dim(tmp)[1] > ntop2) {
		tmp <- tmp[order(tmp$importance,decreasing=T),]
		
		index <- which(tmp$target %in% c('Ccl6', 'Ly6c1', 'Rnf220','Ptprb','Col1a2'))
		index2 <- unique(c(1:(50-length(index)),index))
		
		tmp <- tmp[index2,]
	}

	node2 <- data.frame(tmp$target)
	node2$node.size=tmp$importance
	node2$node.colour <- target_col
	node2$node.colour[which(node2$tmp.target %in% c('Ccl6', 'Ly6c1', 'Rnf220','Ptprb','Col1a2'))] <- motif_target_col
	
	colnames(node2) <- c('node','node.size','node.colour')
	

	df1 <- data.frame(node=tf,node.size=20,node.colour=tf_col)
	node2 <- rbind(df1,node2)
	
	node2$node.colour <- factor(node2$node.colour,levels=c(tf_col,target_col,motif_target_col))
	

	edge2 <- tmp
	colnames(edge2) <- c('from','to','edge.width')
	edge2$edge.colour <- edge_col
	
	torange=c(1,3)
	edge2$edge.width <- scales::rescale(edge2$edge.width,to=torange)


	graph_data <- tidygraph::tbl_graph(nodes = node2, edges = edge2, directed = T)
	
	p1 <- ggraph(graph = graph_data, layout = "stress", circular = TRUE) + 
		geom_edge_arc(aes(edge_colour = edge.colour, edge_width = edge.width)) +
		scale_edge_width_continuous(range = c(1,0.4)) +
		geom_node_point(aes(colour = node.colour, size = node.size))+
		scale_size_continuous(range = c(1, 8)) +		
		theme_void() +
		geom_node_label(aes(label = node,colour = node.colour),size = 3, repel = TRUE)

	p1 <- p1 + 
		scale_color_manual(values=c(tf_col,target_col,motif_target_col))+
		scale_edge_color_manual(values=c(edge_col))
	
}

pdf(paste("05_Age_",inregulons1,"_regulon_netplot.pdf",sep=""),width=8)
p1
dev.off()
	
jpeg(paste("05_Age_",inregulons1,"_regulon_netplot.jpg",sep=""),width=550)
p1
dev.off()




###STAT5B————————————————————————————
ntop2=50

inregulons1=c("STAT5B")
c1 <- which(grn$TF %in% inregulons1)
grn2 <- grn[c1,]


for (tf in inregulons1) {
	
	tf_col <- "#d73027"
	target_col <- "#2873B3"
	motif_target_col <- "#e8743c"
	edge_col <- "#737373"
	
	tmp <- subset(grn2,TF==tf)
	
	if (dim(tmp)[1] > ntop2) {
		tmp <- tmp[order(tmp$importance,decreasing=T),]
		tmp <- tmp[1:ntop2,]
	}

	node2 <- data.frame(tmp$target)
	node2$node.size=tmp$importance
	node2$node.colour <- target_col
	node2$node.colour[which(node2$tmp.target %in% c('KDELC1', 'FAM110A','LINC00937','ZNF267','RNF183','PRB3', 'MCOLN2', 'NPPA','DLGAP5', 'PIAS4','CXCR1','CCL22','PKP3','DHRS3', 'FGF17', 'HAMP', 'LIN28B-AS1', 'ZNF710','TMEM45B','TNFRSF10C'))] <- motif_target_col
	
	colnames(node2) <- c('node','node.size','node.colour')
	

	df1 <- data.frame(node=tf,node.size=20,node.colour=tf_col)
	node2 <- rbind(df1,node2)
	
	node2$node.colour <- factor(node2$node.colour,levels=c(tf_col,target_col,motif_target_col))
	

	edge2 <- tmp
	colnames(edge2) <- c('from','to','edge.width')
	edge2$edge.colour <- edge_col
	
	torange=c(1,3)
	edge2$edge.width <- scales::rescale(edge2$edge.width,to=torange)


	graph_data <- tidygraph::tbl_graph(nodes = node2, edges = edge2, directed = T)
	
	p1 <- ggraph(graph = graph_data, layout = "stress", circular = TRUE) + 
		geom_edge_arc(aes(edge_colour = edge.colour, edge_width = edge.width)) +
		scale_edge_width_continuous(range = c(1,0.4)) +
		geom_node_point(aes(colour = node.colour, size = node.size))+
		scale_size_continuous(range = c(1, 8)) +		
		theme_void() +
		geom_node_label(aes(label = node,colour = node.colour),size = 3, repel = TRUE)

	p1 <- p1 + 
		scale_color_manual(values=c(tf_col,target_col,motif_target_col))+
		scale_edge_color_manual(values=c(edge_col))
	
}

pdf(paste("05_",inregulons1,"_regulon_netplot.pdf",sep=""),width=8)
p1
dev.off()
	
jpeg(paste("05_",inregulons1,"_regulon_netplot.jpg",sep=""),width=550)
p1
dev.off()





