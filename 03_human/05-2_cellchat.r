

cd /d G:\模式所\Li lab\bioinformatics\project\Bone_neu\program2\02_human
scp adata_cellchat.h5ad guest@172.16.4.11:/public/home/guest/scRNA/Bone_neu/



#remotes::install_github("sqjin/CellChat")
#devtools::install_local("/home/luchun/Install/CellChat-master.zip")

library(rhdf5)
library(Matrix)
library(data.table)
library(CellChat)
library(tidyr)
library(Seurat)
library(dplyr)
library(patchwork)
library(gplots)
library(Matrix)
set.seed(101)
library(RColorBrewer)


rm(list=ls())
gc()

setwd("/home/luchun/scRNA/Bone_neu/04_cell/cellchat2/")


mydata <- h5read("mat.h5","mat")

mat <- mydata$block0_values

rownames(mat) <- mydata$axis0
colnames(mat) <- mydata$axis1

mat <- Matrix(mat, sparse = TRUE)
mat[1:4,1:4]
dim(mat)
#[1] 27437 26129



meta <- read.table('metadata.tsv',sep="\t",header=TRUE,row.names=1)
head(meta)

table(meta$preAnno)


identical(rownames(meta),colnames(mat))
identical(match(rownames(meta),colnames(mat)),1:dim(mat)[2])


obj <- CreateSeuratObject(mat,meta.data=meta)
obj
#An object of class Seurat
#27437 features across 26129 samples within 1 assay
#Active assay: RNA (27437 features, 0 variable features)



#————————————————————————cellchat————————————————————————————————
color.vector <- c(brewer.pal(8, "Set1"), brewer.pal(12, "Set3"), brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"))


data.input = obj@assays$RNA@data
meta.data =  obj@meta.data


unique(meta.data$preAnno)


#设置因子水平
meta.data$preAnno = factor(meta.data$preAnno,levels = c("Neu_AZU1+","Neu_LTF+","Neu_SELL+","Neu_IL1B+",'Fibro-MSC','THY1+ MSC','Adipo-MSC','APOD+ MSC','Osteo-MSC','Osteoblast'))



#——————————————————————Step1. 构建cellchat对象——————————————————
### 1.2 Create a CellChat object
cellchat <- createCellChat(object = data.input, 
                           meta = meta.data, 
                           group.by = "preAnno")


### 1.3 可在cellchat对象的meta插槽中添加表型信息
cellchat <- addMeta(cellchat, meta = meta.data)

# 设置默认的labels
levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents)) 




#———————————————————Step2. 加载CellChatDB数据库———————————————————
### 1.4 加载CellChat受配体数据库
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

jpeg("01_showDatabaseCategory.jpg",width=600,height=400)
showDatabaseCategory(CellChatDB)
dev.off()


# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)


CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use



#———————————————————Step3. 对表达数据进行预处理———————————————————
cellchat <- subsetData(cellchat) 
#future::plan("multiprocess", workers = 4)# do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#cellchat <- projectData(cellchat, PPI.human) 
#找到配体受体关系后，projectData将配体受体对的表达值投射到PPI上




#——————————————Step4. 计算通讯概率，推断细胞通讯网络——————————————
#population. size = FALSE：如果分析分选富集的单细胞，以消除潜在的人为因素细胞群的大小。
#population. size = TRUE，如果分析非分选富集的单细胞转录组，原因是丰富的细胞群往往比稀少的细胞群发送集体更强的信号。

cellchat <- computeCommunProb(cellchat, type = "truncatedMean",trim = 0.1,population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)


saveRDS(cellchat,"cellchat.rds")


#cellchat <- readRDS("cellchat.rds")



#—————————————————————— cellchat取子集————————————————————————————————
subsetCellChat <- function (object, cells.use = NULL, idents.use = NULL, group.by = NULL,invert = FALSE, thresh = 0.05){
    if (!is.null(idents.use)) {
        if (is.null(group.by)) {
            labels <- object@idents
            if (object@options$mode == "merged") {
                message("Use the joint cell labels from the merged CellChat object")
                labels <- object@idents$joint
            }
        }
        else {
            labels <- object@meta[[group.by]]
        }
        if (!is.factor(labels)) {
            labels <- factor(labels)
        }
        level.use0 <- levels(labels)
        level.use <- levels(labels)[levels(labels) %in% unique(labels)]
        if (invert) {
            level.use <- level.use[!(level.use %in% idents.use)]
        }
        else {
            level.use <- level.use[level.use %in% idents.use]
        }
        cells.use.index <- which(as.character(labels) %in% level.use)
        cells.use <- names(labels)[cells.use.index]
    }
    else if (!is.null(cells.use)) {
        labels <- object@idents
        if (object@options$mode == "merged") {
            message("Use the joint cell labels from the merged CellChat object")
            labels <- object@idents$joint
        }
        level.use0 <- levels(labels)
        level.use <- levels(labels)[levels(labels) %in% unique(as.character(labels[cells.use]))]
        cells.use.index <- which(names(labels) %in% cells.use)
    }
    else {
        stop("USER should define either `cells.use` or `idents.use`!")
    }
    cat("The subset of cell groups used for CellChat analysis are ",
        level.use, "\n")
    if (nrow(object@data) > 0) {
        data.subset <- object@data[, cells.use.index]
    }
    else {
        data.subset <- matrix(0, nrow = 0, ncol = 0)
    }
    if (nrow(object@data.project) > 0) {
        data.project.subset <- object@data.project[, cells.use.index]
    }
    else {
        data.project.subset <- matrix(0, nrow = 0, ncol = 0)
    }
    data.signaling.subset <- object@data.signaling[, cells.use.index]
    meta.subset <- object@meta[cells.use.index, , drop = FALSE]
    if (object@options$mode == "merged") {
        idents <- object@idents[1:(length(object@idents) - 1)]
        group.existing <- level.use0[level.use0 %in% level.use]
        group.existing.index <- which(level.use0 %in% level.use)
        net.subset <- vector("list", length = length(object@net))
        netP.subset <- vector("list", length = length(object@netP))
        idents.subset <- vector("list", length = length(idents))
        names(net.subset) <- names(object@net)
        names(netP.subset) <- names(object@netP)
        names(idents.subset) <- names(object@idents[1:(length(object@idents) -
            1)])
        images.subset <- vector("list", length = length(idents))
        names(images.subset) <- names(object@idents[1:(length(object@idents) -
            1)])
        for (i in 1:length(idents)) {
            cat("Update slots object@images, object@net, object@netP, object@idents in dataset ",
                names(object@idents)[i], "\n")
            images <- object@images[[i]]
            for (images.j in names(images)) {
                values <- images[[images.j]]
                if (images.j %in% c("coordinates")) {
                  values.new <- values[cells.use.index, ]
                  images[[images.j]] <- values.new
                }
                if (images.j %in% c("distance")) {
                  values.new <- values[group.existing.index,
                    group.existing.index, drop = FALSE]
                  images[[images.j]] <- values.new
                }
            }
            images.subset[[i]] <- images
            net <- object@net[[i]]
            for (net.j in names(net)) {
                values <- net[[net.j]]
                if (net.j %in% c("prob", "pval")) {
                  values.new <- values[group.existing.index,
                    group.existing.index, ]
                  net[[net.j]] <- values.new
                }
                if (net.j %in% c("count", "sum", "weight")) {
                  values.new <- values[group.existing.index,
                    group.existing.index]
                  net[[net.j]] <- values.new
                }
            }
            net.subset[[i]] <- net
            netP = computeCommunProbPathway(net = net.subset[[i]],
                pairLR.use = object@LR[[i]]$LRsig, thresh = thresh)
            netP$centr = netAnalysis_computeCentrality(net = net.subset[[i]]$prob)
            netP.subset[[i]] <- netP
            idents.subset[[i]] <- idents[[i]][names(idents[[i]]) %in%
                cells.use]
            idents.subset[[i]] <- factor(idents.subset[[i]],
                levels = levels(idents[[i]])[levels(idents[[i]]) %in%
                  level.use])
        }
        idents.subset$joint <- factor(object@idents$joint[cells.use.index],
            levels = level.use)
    }
    else {
        cat("Update slots object@images, object@net, object@netP in a single dataset...",
            "\n")
        group.existing <- level.use0[level.use0 %in% level.use]
        group.existing.index <- which(level.use0 %in% level.use)
        images <- object@images
        for (images.j in names(images)) {
            values <- images[[images.j]]
            if (images.j %in% c("coordinates")) {
                values.new <- values[cells.use.index, ]
                images[[images.j]] <- values.new
            }
            if (images.j %in% c("distance")) {
                values.new <- values[group.existing.index, group.existing.index,
                  drop = FALSE]
                images[[images.j]] <- values.new
            }
        }
        images.subset <- images
        net <- object@net
        for (net.j in names(net)) {
            values <- net[[net.j]]
            if (net.j %in% c("prob", "pval")) {
                values.new <- values[group.existing.index, group.existing.index, 1:dim(values)[3], drop = FALSE]
                net[[net.j]] <- values.new
            }
            if (net.j %in% c("count", "sum", "weight")) {
                values.new <- values[group.existing.index, group.existing.index, drop = FALSE]
                net[[net.j]] <- values.new
            }
        }
        net.subset <- net
        netP = computeCommunProbPathway(net = net.subset, pairLR.use = object@LR$LRsig,
            thresh = thresh)
        netP$centr = netAnalysis_computeCentrality(net = net.subset$prob)
        netP.subset <- netP
        idents.subset <- object@idents[cells.use.index]
        idents.subset <- factor(idents.subset, levels = level.use)
    }
    object.subset <- methods::new(Class = "CellChat", data = data.subset,
        data.signaling = data.signaling.subset, data.project = data.project.subset,
        images = images.subset, net = net.subset, netP = netP.subset,
        meta = meta.subset, idents = idents.subset, var.features = object@var.features,
        LR = object@LR, DB = object@DB, options = object@options)
    return(object.subset)
}



cellchat.subset <- cellchat
#cellchat.subset <- subsetCellChat(object=cellchat, idents.use=c("Neu_LTF+","Neu_SELL+","Neu_IL1B+",'Fibro-MSC','THY1+ MSC','Adipo-MSC','APOD+ MSC','Osteo-MSC','Osteoblast'))




##——————————————Step6. 在信号通路水平推断细胞通讯————————————

#CellChat通过汇总与每个信号通路相关的所有配体-受体相互作用的通讯概率，计算信号通路层面的通讯概率。
cellchat.subset <- computeCommunProbPathway(cellchat.subset)
#head(cellchat.subset@net)
#head(cellchat.subset@netP)



#——————————————Step7. 计算加和的cell-cell通讯网络————————————
cellchat.subset <- aggregateNet(cellchat.subset)

groupSize <- as.numeric(table(cellchat.subset@idents))



levels(cellchat.subset@idents)
# [1] "Neu_AZU1+"  "Neu_LTF+"   "Neu_SELL+"  "Neu_IL1B+"  "Fibro-MSC"
# [6] "THY1+ MSC"  "Adipo-MSC"  "APOD+ MSC"  "Osteo-MSC"  "Osteoblast"


sources <- 1:4
targets <- 5:10


#cols <- c("#009E73", "#F0E442", "#f06152","#0072B2","#fc9272","#428AC9","#129392","#FFCC4F","#F37E78","#883A96","#cca69c")
cols <- c("#5066a1","#7fb961","#ffc556","#f06152",
"#428AC9","#129392","#FFCC4F","#F37E78","#883A96","#cca69c")



pdf("02-1_circle_weights.pdf")
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat.subset@net$weight,color.use = cols,sources.use = sources, targets.use = targets, vertex.weight = groupSize,edge.width.max = 12,weight.scale = T, label.edge= F,title.name = "Outgoing interaction strength of neutrophil",arrow.size = 0.6)
dev.off()


#反过来，source变成targets
pdf("02-2_circle_weights.pdf")
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat.subset@net$weight,color.use = cols,sources.use = targets, targets.use = sources, vertex.weight = groupSize,edge.width.max = 12,weight.scale = T, label.edge= F,title.name = "Incoming interaction strength of neutrophil",arrow.size = 0.6)
dev.off()





dat1 <- cellchat.subset@net$weight[sources,targets]

dat <- dat1

my_palette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(n = 100)
mat_rounded <- round(dat, digits = 2)

# 创建热图并标注数值
pdf("02-3_heatmap_out_weight.pdf",width=9,height=9)
par(cex.main=1)
heatmap.2(dat,   # 你的数值矩阵
         dendrogram = "both",   # 在行和列上显示聚类图
          trace = "none",        # 禁用行列边缘的注释
          col = my_palette,      # 指定颜色谱
          cellnote = mat_rounded,        # 标注每个单元格的数值
          notecol = "black",     # 标注文本的颜色
          density.info = "none", # 不显示密度信息
          key = TRUE,            # 显示颜色键
          keysize = 1.0,         # 颜色键大小
          margins = c(5, 5),     # 设置边距
          cexCol = 0.8,          # 列标签字体大小
          cexRow = 0.8,          # 行标签字体大小
          xlab = "Target",       # x轴标签
          ylab = "Source",       # y轴标签
          main = "Outgoing interaction strength of neutrophil"
)
dev.off()





dat2 <- cellchat.subset@net$weight[targets,sources]
dat <- t(dat2)

my_palette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(n = 100)
mat_rounded <- round(dat, digits = 2)

# 创建热图并标注数值
pdf("02-4_heatmap_in_weight.pdf",width=9,height=9)
par(cex.main=1)
heatmap.2(dat,   # 你的数值矩阵
         dendrogram = "both",   # 在行和列上显示聚类图
          trace = "none",        # 禁用行列边缘的注释
          col = my_palette,      # 指定颜色谱
          cellnote = mat_rounded,        # 标注每个单元格的数值
          notecol = "black",     # 标注文本的颜色
          density.info = "none", # 不显示密度信息
          key = TRUE,            # 显示颜色键
          keysize = 1.0,         # 颜色键大小
          margins = c(5, 5),     # 设置边距
          cexCol = 0.8,          # 列标签字体大小
          cexRow = 0.8,          # 行标签字体大小
          xlab = "Source",       # x轴标签
          ylab = "Target",       # y轴标签
          main = "Incoming interaction strength of neutrophil"
)
dev.off()





dat1 <- cellchat.subset@net$weight[sources,targets]
dat2 <- cellchat.subset@net$weight[targets,sources]

outgoing_dat <- rowSums(dat1)
incoming_dat <- colSums(dat2)

dat <- data.frame(
  outgoing = outgoing_dat,      
  incoming = incoming_dat,
  size=groupSize[sources]  
)

dat$preAnno <- rownames(dat)
dat$preAnno <- factor(dat$preAnno,levels=c("Neu_AZU1+","Neu_LTF+","Neu_SELL+","Neu_IL1B+"))

dat$outgoing <- dat$outgoing*1000
dat$incoming <- dat$incoming*1000

dat

# 使用ggplot绘制散点图
cols <- c("#5066a1","#7fb961","#ffc556","#f06152")
 

pdf("02-5_dot_in_out.pdf",height=9,width=11)
ggplot(dat, aes(x = outgoing, y = incoming, size=size,color = preAnno)) +
  geom_point() +
  scale_size_continuous(range = c(8, 20)) +
  scale_x_continuous(limits = c(0, 9),breaks=seq(0,9,2))+
  scale_y_continuous(limits = c(0, 50),breaks=seq(0,50,10))+
  geom_text(aes(label = preAnno), hjust = 0.5, vjust = 5,size=5)+
  labs(x = "Outgoing interaction strength", y = "Incoming interaction strength") +
  scale_color_manual(values=cols)+
  theme_bw() +
  theme(aspect.ratio = 1/1,panel.grid=element_blank(),
  axis.title.x = element_text(size = 20),  # 设置x轴标签字体大小
    axis.title.y = element_text(size = 20),  # 设置y轴标签字体大小
    axis.text.x = element_text(size = 15),    # 设置x轴刻度字体大小
    axis.text.y = element_text(size = 15),     # 设置y轴刻度字体大小
	legend.title = element_text(size = 15),    # 设置图例标题字体大小
    legend.text = element_text(size = 12))
dev.off() 
  
  




##——————————————————————Step11. 观察多种配体受体或信号通路介导的细胞-细胞通讯——————————————————

### Bubble plot
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')

pdf("02-8_dotplot.pdf",width=10,height=9)
netVisual_bubble(cellchat.subset,sources.use = sources, targets.use = targets, font.size = 13,angle.x = 45,remove.isolate = FALSE)
dev.off()

pdf("02-9_dotplot.pdf",width=10,height=12)
netVisual_bubble(cellchat.subset,sources.use = targets, targets.use = sources, font.size = 13,angle.x = 45,remove.isolate = FALSE)
dev.off()					 
					 


'''
#——————————————————————————Step13. 识别细胞亚群的信号作用（例如主要的发送者，接收者）以及主要的贡献信号————————————————————————————

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
cellchat@netP$pathways


pdf("04_signaling.pdf")
netAnalysis_signalingRole_scatter(cellchat)+ggtitle("All pathway")
dev.off()


ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",height = 12,font.size = 5)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",height = 12,font.size = 5)


pdf("04_signaling.pdf",width=14,height=14)
ht1 + ht2
dev.off()


pdf("04_signaling.pdf",width=11.67,height=11.67)
ht1 + ht2
dev.off()


jpeg("05_circle_plot.jpg",width=600)
netVisual_chord_gene(cellchat, sources.use = c(19,25), targets.use = c(9:12), signaling = c("COLLAGEN","LAMININ","APP","FN1"), lab.cex = 0.5,legend.pos.x = 8)
dev.off()



pdf("05_circle_plot.pdf",width=8.75)
netVisual_chord_gene(cellchat, sources.use = c(19,25), targets.use = c(9:12), signaling = c("COLLAGEN","LAMININ","APP","FN1"), lab.cex = 0.5,legend.pos.x = 8)
dev.off()



jpeg("06_circle_plot.jpg",width=600)
netVisual_chord_gene(cellchat, sources.use = c(13), targets.use = c(1:12), signaling = c("ANNEXIN","SELL"), legend.pos.x = 8)
dev.off()

#,"SELPLG","RESISTIN","TGFb","CDH1"


'''




#———————————————————————————————————————————————————————————————————
cellchat <- readRDS("cellchat.rds")


df.net <- subsetCommunication(cellchat)
head(df.net)

#以通路为单位提取通讯信息
df.pathway = subsetCommunication(cellchat,slot.name = "netP")

levels(cellchat@idents)
# [1] "AZU1+KCNQ5+ progenitor"   "PDE4D+MCU+ proliferating"
# [3] "LTF+CAMP+ immature"       "S100A12+MMP9+ immature"
# [5] "S100A6+NCF1+ mature"      "NAMPT+IFITM2+ mature"
# [7] "NEAT1+FTH1+ mature"       "Fibro-MSC"
# [9] "THY1+ MSC"                "Adipo-MSC"
#[11] "APOD+ MSC"                "Osteo-MSC"
#[13] "Osteoblast"



sources <- 1:7
targets <- 8:13


df.net <- subsetCommunication(cellchat, sources.use = sources, targets.use = targets)
head(df.net)


df.net <- subsetCommunication(cellchat, signaling = c("MIF", "TGFb"))
head(df.net)




##——————————————Step6. 在信号通路水平推断细胞通讯————————————

#CellChat通过汇总与每个信号通路相关的所有配体-受体相互作用的通讯概率，计算信号通路层面的通讯概率。
cellchat <- computeCommunProbPathway(cellchat)
#head(cellchat@net)
#head(cellchat@netP)



#——————————————Step7. 计算加和的cell-cell通讯网络————————————
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

cols <- c("#76afda","#8264CC","#009E73", "#F0E442", "#f06152","#0072B2","#fc9272","#428AC9","#129392","#FFCC4F","#F37E78","#883A96","#cca69c")



pdf("03-1_circle_Number.pdf")
netVisual_circle(cellchat@net$count,color.use = cols,sources.use = sources, targets.use = targets, vertex.weight = groupSize,edge.width.max = 12,weight.scale = T, label.edge= F,title.name = "Number of interactions",arrow.size = 0.6)
dev.off()

#反过来，source变成targets
pdf("03-2_circle_Number.pdf")
netVisual_circle(cellchat@net$count,color.use = cols,sources.use = targets, targets.use = sources, vertex.weight = groupSize,edge.width.max = 12,weight.scale = T, label.edge= F,title.name = "Number of interactions",arrow.size = 0.6)
dev.off()



pdf("03-3_circle_weights.pdf")
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$weight,color.use = cols,sources.use = sources, targets.use = targets, vertex.weight = groupSize,edge.width.max = 12,weight.scale = T, label.edge= F,title.name = "Outgoing interaction strength of neutrophil",arrow.size = 0.6)
dev.off()


#反过来，source变成targets
pdf("03-4_circle_weights.pdf")
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$weight,color.use = cols,sources.use = targets, targets.use = sources, vertex.weight = groupSize,edge.width.max = 12,weight.scale = T, label.edge= F,title.name = "Incoming interaction strength of neutrophil",arrow.size = 0.6)
dev.off()





dat1 <- cellchat@net$weight[sources,targets]

dat <- dat1

my_palette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(n = 100)
mat_rounded <- round(dat, digits = 2)

# 创建热图并标注数值
pdf("03-5_heatmap_out_weight.pdf",width=9,height=9)
par(cex.main=1)
heatmap.2(dat,   # 你的数值矩阵
         dendrogram = "both",   # 在行和列上显示聚类图
          trace = "none",        # 禁用行列边缘的注释
          col = my_palette,      # 指定颜色谱
          cellnote = mat_rounded,        # 标注每个单元格的数值
          notecol = "black",     # 标注文本的颜色
          density.info = "none", # 不显示密度信息
          key = TRUE,            # 显示颜色键
          keysize = 1.0,         # 颜色键大小
          margins = c(5, 5),     # 设置边距
          cexCol = 0.8,          # 列标签字体大小
          cexRow = 0.8,          # 行标签字体大小
          xlab = "Target",       # x轴标签
          ylab = "Source",       # y轴标签
          main = "Outgoing interaction strength of neutrophil"
)
dev.off()





dat2 <- cellchat@net$weight[targets,sources]
dat <- t(dat2)

my_palette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(n = 100)
mat_rounded <- round(dat, digits = 2)

# 创建热图并标注数值
pdf("03-6_heatmap_in_weight.pdf",width=9,height=9)
par(cex.main=1)
heatmap.2(dat,   # 你的数值矩阵
         dendrogram = "both",   # 在行和列上显示聚类图
          trace = "none",        # 禁用行列边缘的注释
          col = my_palette,      # 指定颜色谱
          cellnote = mat_rounded,        # 标注每个单元格的数值
          notecol = "black",     # 标注文本的颜色
          density.info = "none", # 不显示密度信息
          key = TRUE,            # 显示颜色键
          keysize = 1.0,         # 颜色键大小
          margins = c(5, 5),     # 设置边距
          cexCol = 0.8,          # 列标签字体大小
          cexRow = 0.8,          # 行标签字体大小
          xlab = "Source",       # x轴标签
          ylab = "Target",       # y轴标签
          main = "Incoming interaction strength of neutrophil"
)
dev.off()





dat1 <- cellchat@net$weight[sources,targets]
dat2 <- cellchat@net$weight[targets,sources]

outgoing_dat <- rowSums(dat1)
incoming_dat <- colSums(dat2)

dat <- data.frame(
  outgoing = outgoing_dat,      
  incoming = incoming_dat,
  size=groupSize[sources]  
)

dat$preAnno <- rownames(dat)
dat$preAnno <- factor(dat$preAnno,levels=c("AZU1+KCNQ5+ progenitor","PDE4D+MCU+ proliferating","LTF+CAMP+ immature", "S100A12+MMP9+ immature",'S100A6+NCF1+ mature',"NAMPT+IFITM2+ mature","NEAT1+FTH1+ mature"))

dat$outgoing <- dat$outgoing*1000
dat$incoming <- dat$incoming*1000


# 使用ggplot绘制散点图
cols <- c("#76afda","#8264CC","#009E73", "#F0E442", "#f06152","#0072B2","#fc9272")


pdf("03-7_dot_in_out.pdf",height=9,width=11)
ggplot(dat, aes(x = outgoing, y = incoming, size=size,color = preAnno)) +
  geom_point() +
  scale_size_continuous(range = c(20, 30)) +
  scale_x_continuous(limits = c(0, 9),breaks=seq(0,8,2))+
  scale_y_continuous(limits = c(0, 26),breaks=seq(0,25,5))+
  geom_text(aes(label = preAnno), hjust = 0.5, vjust = 5,size=5)+
  labs(x = "Outgoing interaction strength", y = "Incoming interaction strength") +
  scale_color_manual(values=cols)+
  theme_bw()+
  theme(panel.grid=element_blank())
dev.off() 
  
  




##——————————————————————Step11. 观察多种配体受体或信号通路介导的细胞-细胞通讯——————————————————

### Bubble plot
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')

pairLR.use <- extractEnrichedLR(cellchat, signaling = c("TGFb"))
pairLR.use <- data.frame(interaction_name=c('HBEGF_ERBB4', 'HBEGF_EGFR', 'HBEGF_ERBB2_ERBB4', 'CTSG_PARD3', 'AREG_EGFR', 'AREG_EGFR_ERBB2','TGFB1_TGFBR1_TGFBR2', 'TGFB1_ACVR1C_TGFBR2', 'TGFB1_ACVR1B_TGFBR2', 'TGFB1_ACVR1_TGFBR1', 'NAMPT_INSR', 'NAMPT_ITGA5_ITGB1', 'BMP6_BMPR1B_BMPR2', 'BMP6_BMPR1B_ACVR2A', 'BMP6_BMPR1A_BMPR2', 'BMP6_BMPR1A_ACVR2B', 'BMP6_BMPR1A_ACVR2A', 'BMP6_ACVR1_BMPR2', 'BMP6_ACVR1_ACVR2B', 'BMP6_ACVR1_ACVR2A'))

pdf("03-8_dotplot.pdf",width=10,height=5)
netVisual_bubble(cellchat,sources.use = sources, targets.use = targets,pairLR.use=pairLR.use, font.size = 10,angle.x = 45,remove.isolate = FALSE)
dev.off()

pdf("03-9_dotplot.pdf",width=10,height=18)
netVisual_bubble(cellchat,sources.use = targets, targets.use = sources, font.size = 13,angle.x = 45,remove.isolate = FALSE)
dev.off()					 
					 



'''


#——————————————————————————Step13. 识别细胞亚群的信号作用（例如主要的发送者，接收者）以及主要的贡献信号————————————————————————————

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
cellchat@netP$pathways


pdf("04_signaling.pdf")
netAnalysis_signalingRole_scatter(cellchat)+ggtitle("All pathway")
dev.off()


ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",height = 12,font.size = 5)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",height = 12,font.size = 5)


pdf("04_signaling.pdf",width=14,height=14)
ht1 + ht2
dev.off()


pdf("04_signaling.pdf",width=11.67,height=11.67)
ht1 + ht2
dev.off()


jpeg("05_circle_plot.jpg",width=600)
netVisual_chord_gene(cellchat, sources.use = c(19,25), targets.use = c(9:12), signaling = c("COLLAGEN","LAMININ","APP","FN1"), lab.cex = 0.5,legend.pos.x = 8)
dev.off()



pdf("05_circle_plot.pdf",width=8.75)
netVisual_chord_gene(cellchat, sources.use = c(19,25), targets.use = c(9:12), signaling = c("COLLAGEN","LAMININ","APP","FN1"), lab.cex = 0.5,legend.pos.x = 8)
dev.off()



jpeg("06_circle_plot.jpg",width=600)
netVisual_chord_gene(cellchat, sources.use = c(13), targets.use = c(1:12), signaling = c("ANNEXIN","SELL"), legend.pos.x = 8)
dev.off()

#,"SELPLG","RESISTIN","TGFb","CDH1"
