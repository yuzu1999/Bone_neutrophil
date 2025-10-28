


library(SCopeLoomR) 
#library(devtools)
#devtools::install_local("/home/luchun/Install/SCopeLoomR-master.zip")


library(AUCell)
library(SCENIC) 
#devtools::install_local("/home/luchun/Install/SCENIC-master.zip")


library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(ggplot2)
library(pheatmap)
library(Seurat)
library(dior)
library(SCENIC)


rm(list=ls())
gc()

setwd("/home/luchun/scRNA/Bone_neu/03_age_ovx2/pyscenic/")



sce_SCENIC <- open_loom("./data/OVX_total_SCENIC.loom")

regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
class(regulons)
ls(regulons)




len=c()
for(i in ls(regulons)){
	tmp <- length(regulons[[i]])
	len <- c(len,tmp)
}

len_names <- data.frame(tfs=names(regulons),len=paste(" (",len," genes)",sep=""))

len_names$new_name <- paste(substr(len_names$tfs, 1, nchar(len_names$tfs) - 3),len_names$len,sep="")

write.csv(len_names,"OVX_len_names.csv",quote=FALSE,row.names=FALSE)






regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC')


pydata <- dior::read_h5('OVX_total_adata.h5')
table(pydata$preAnno)


pydata$preAnno2 <- "Neu"
pydata$preAnno2[which(pydata$preAnno == "non-neutrophil")] <- "Non_neu"

table(pydata$preAnno2)




pydata$preAnno3 <- "Non_neu"
pydata$preAnno3[which(pydata$preAnno == 'Elane+Prtn3+ progenitor')] <- "Neu_1"
pydata$preAnno3[which(pydata$preAnno == 'H2afz+Hmgb2+ proliferating')] <- "Neu_2"
pydata$preAnno3[which(pydata$preAnno == 'Ngp+Lcn2+ immature')] <- "Neu_3"
pydata$preAnno3[which(pydata$preAnno == 'Ifitm6+Ltf+ immature')] <- "Neu_4"
pydata$preAnno3[which(pydata$preAnno == 'Retnlg+Mmp8+ mature')] <- "Neu_5"
pydata$preAnno3[which(pydata$preAnno == 'Ccl6+Sell+ mature')] <- "Neu_6"
pydata$preAnno3[which(pydata$preAnno == 'Il1b+Srgn+ mature')] <- "Neu_7"
table(pydata$preAnno3)



cellinfo <- pydata@meta.data[,c('preAnno','preAnno2','preAnno3')] 
colnames(cellinfo) <- c('preAnno','preAnno2','preAnno3')

cellinfo$cell <- rownames(cellinfo)
dim(cellinfo)
#[1] 40605     4



head(cellinfo)

write.csv(cellinfo,"OVX_cellinfo.csv",quote=FALSE)




#——————————————————计算细胞特意性的转录因子————————————————————————

#通过点图展示————————————————————————

cellTypes <-  as.data.frame(subset(cellinfo,select = 'preAnno'))
selectedResolution <- "preAnno"
sub_regulonAUC <- regulonAUC


identical(rownames(cellinfo),colnames(getAUC(sub_regulonAUC)))

AUC=getAUC(sub_regulonAUC)
dim(AUC)
#[1]   355 40605

write.csv(AUC,"OVX_activity.csv",quote=FALSE)





#计算活性均值————————————————————————————
cellTypes <- data.frame(row.names = colnames(pydata), 
                        celltype = pydata$preAnno)


selectedResolution <- "celltype" # select resolution
cellsPerGroup <- split(rownames(cellTypes), cellTypes[,selectedResolution])
					   
					   
# 去除extened regulons
sub_regulonAUC <- regulonAUC

sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
dim(sub_regulonAUC)
#[1]   355 40605


regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))
		
dim(regulonActivity_byGroup)
#[1] 355   8


write.csv(regulonActivity_byGroup,"OVX_activity_average.csv",quote=FALSE)



'''
regulonActivity_byGroup2 <- regulonActivity_byGroup[,c(2,3,6,4,7,1,5)]


result <- regulonActivity_byGroup2[regulons.use,]
head(result)





#画热图————————————————————————————————————
typelist = c("AZU1+KCNQ5+ progenitor","PDE4D+MCU+ proliferating","LTF+CAMP+ immature","S100A12+MMP9+ immature","S100A6+NCF1+ mature","NAMPT+IFITM2+ mature","NEAT1+FTH1+ mature")

col_ann <- data.frame(preAnno=typelist)
rownames(col_ann) <- typelist



ann_colors <- list(
  preAnno = c(`AZU1+KCNQ5+ progenitor` = "#76afda",`PDE4D+MCU+ proliferating` = "#8264CC",`LTF+CAMP+ immature` = "#009E73",`S100A12+MMP9+ immature` = "#F0E442",`S100A6+NCF1+ mature`="#f06152",`NAMPT+IFITM2+ mature`="#0072B2",`NEAT1+FTH1+ mature`="#fc9272"))

 
cols <- colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100)



jpeg("02_AS_data_heatmap.jpg")
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
		main="Activity Score",
		angle_col ="45")
dev.off()



pdf("02_AS_data_heatmap.pdf")
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
		main="Activity Score",
		angle_col ="45")
dev.off()








#利用热图展示细胞类型特异性转录因子————————————————————————
library(pheatmap)
library(reshape2)

rss_data <- rssPlot$plot$data
rss_data 

rss_data <- dcast(rss_data,Topic~rss_data$preAnno,value.var = 'Z')
rss_data 

rownames(rss_data) <- rss_data[,1]
rss_data <- rss_data[,-1]


typelist = c("AZU1+KCNQ5+ progenitor","PDE4D+MCU+ proliferating","LTF+CAMP+ immature","S100A12+MMP9+ immature","S100A6+NCF1+ mature","NAMPT+IFITM2+ mature","NEAT1+FTH1+ mature")

col_ann <- data.frame(preAnno=typelist)
rownames(col_ann) <- typelist



ann_colors <- list(
  preAnno = c(`AZU1+KCNQ5+ progenitor` = "#76afda",`PDE4D+MCU+ proliferating` = "#8264CC",`LTF+CAMP+ immature` = "#009E73",`S100A12+MMP9+ immature` = "#F0E442",`S100A6+NCF1+ mature`="#f06152",`NAMPT+IFITM2+ mature`="#0072B2",`NEAT1+FTH1+ mature`="#fc9272"))

 
cols <- colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100)



jpeg("03_rssPlot_data_heatmap.jpg",height=800,width=600)
pheatmap(rss_data,scale = "row",
		cluster_cols=FALSE,cluster_rows=TRUE,
		show_colnames=TRUE,show_rownames=TRUE,
		annotation_col=col_ann,
		annotation_colors=ann_colors,
		fontsize_row=10,
		fontsize_col=10,
		cellheight=15,cellwidth=20,
		color=cols,
		treeheight_row=20,treeheight_col=20,
		border_color="white", 
		angle_col ="45")
dev.off()






#————————————————————————独立展示 rank后 转录因子——————————————————

# regulonActivity_byGroup
randdf = rss2
pplist = list()


for( i in  1:dim(randdf)[2]){
    groupcol = colnames(randdf)[i]
    rankdf = randdf[,groupcol] %>% na.omit()
    randf = data.frame(tf = names(rankdf),scorez = rankdf )
	
    randf = randf[order(randf$scorez,decreasing=TRUE),]
    randf$xx = seq(1,dim(randf)[1])
    randf$colo = rep("#76afda",dim(randf)[1])
    randf$colo[1:10] = "#d73027"
    randf$tf[11:dim(randf)[1]] = NA
	
	
    library(ggrepel)
    pp = ggplot(data=randf,aes(x=xx,y=scorez))+geom_point(aes(color=colo),size=3)+
	scale_color_manual(values = c("#76afda","#d73027")) +
    geom_text_repel(data=randf,aes(x=xx,label=tf),max.overlaps = 15,nudge_x = 50-randf$xx,direction = "y",hjust = 0)+
      labs(y="Specificity Score",x="Regulons",title=groupcol)+
      theme_bw()+
	  theme(panel.grid=element_blank())+
      scale_x_continuous(breaks = NULL)+
      theme(legend.position = "none",
      plot.title = element_text( hjust=0.5, size =16))
	  
    pplist[[i]] = pp
}

gglist = patchwork::wrap_plots(plots = pplist,ncol = 4)
ggsave(file="04_rank_plot_list.pdf",gglist,width=16,height=8,limitsize=F)




