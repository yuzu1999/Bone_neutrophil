


rm(list=ls())
gc()

setwd("/home/luchun/scRNA/Bone_neu/02_ovx2/")

gene_mapping <- read.csv("/home/luchun/scRNA/Bone_neu/mouse_human_orthologs.csv",header=TRUE)


convert_mouse_to_human <- function(mouse_genes) {
  # 创建一个空向量保存人类同源基因
  human_genes <- character(length(mouse_genes))
  
  # 遍历每个鼠基因,查找对应的人类同源基因
  for (i in seq_along(mouse_genes)) {
    mouse_gene <- mouse_genes[i]
    matching_rows <- which(gene_mapping$mouse_gene == mouse_gene)
    
    # 如果找到一个一对一的对应关系,保存人类基因名
    if (length(matching_rows) == 1) {
      human_genes[i] <- gene_mapping$human_gene[matching_rows]
    } else {
      human_genes[i] <- NA # 否则标记为 NA
    }
  }
  
  return(human_genes)
}




Mus_deg <- read.csv("OVX_Neu_preAnno4_deg.csv",header=TRUE)
Mus_deg_sig <- subset(Mus_deg,logfoldchanges > 0 & pvals_adj < 0.05)



library(dplyr)
Mus_top100_2 <- list()

for(i in unique(Mus_deg_sig$group)){
	genes <- Mus_deg_sig$names[Mus_deg_sig$group == i][1:300]
	
	my_vector <- convert_mouse_to_human(genes)
	
	cleaned_vector <- subset(my_vector, !is.na(my_vector))
	
	name <- paste("Mice_",i,sep="")
	Mus_top100_2[[name]] <- cleaned_vector

}
lengths(Mus_top100_2)
#head(Mus_top100_2[["Mice_Neu_1"]])



Homo_deg <- read.csv("/home/luchun/scRNA/Bone_neu/04_cell/neu/Human_Neu_preAnno4_deg.csv",header=TRUE)
Homo_deg_sig <- subset(Homo_deg,logfoldchanges > 0 & pvals_adj < 0.05)
head(Homo_deg_sig)


Homo_top100_2 <- list()

for(i in unique(Homo_deg_sig$group)){
	genes <- Homo_deg_sig$names[Homo_deg_sig$group == i][1:300]

	name <- paste("Human_",i,sep="")
	Homo_top100_2[[name]] <- genes

}
lengths(Homo_top100_2)
#head(Homo_top100_2[["Human_Neu_1"]])


shared_genes <- lapply(Homo_top100_2, function(x) {
   lapply(Mus_top100_2, function(y) {
     length(intersect(x, y))
   })
})


shared_genes


tmp <- c()
for(i in 1:length(Homo_top100_2)){	
	for(j in 1:length(Mus_top100_2)){
		inter <- intersect(Homo_top100_2[[i]],Mus_top100_2[[j]])
        tmp <- c(tmp,length(inter))
	}
}

tmp

#这里的数字要更改，列数是人的neu类别数
tmp <- matrix(tmp,ncol=4)

colnames(tmp) <- names(Homo_top100_2)
rownames(tmp) <- names(Mus_top100_2)

#colnames(tmp) <- paste("Human_",1:9,sep="")
#rownames(tmp) <- paste("Mice_",1:7,sep="")

tmp



# 假设这是你的除数向量


divisors <- lengths(Mus_top100_2)

# 对每一行进行除法操作
for (i in 1:nrow(tmp)) {
  tmp[i, ] <- tmp[i, ] / divisors[i]*100
}


tmp




library(pheatmap)
library(RColorBrewer)
mycol <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(100)

pdf("UP_similarity.pdf",width=9)
pheatmap(tmp,
		cluster_rows = FALSE,
		cluster_cols = FALSE,
		border_color = "white",
		display_numbers = TRUE,
		fontsize_number = 12, 
		color=mycol,
		cellwidth =40, cellheight=40,
		number_format = "%.2f",
		number_color = "black",
		angle_col=45,
		main="Up-regulated genes similarity")
dev.off()





#——————————————————————————————  下调  —————————————————————————————

rm(list=ls())
gc()

setwd("/home/luchun/scRNA/Bone_neu/01_age/")

gene_mapping <- read.csv("/home/luchun/scRNA/Bone_neu/mouse_human_orthologs.csv",header=TRUE)


convert_mouse_to_human <- function(mouse_genes) {
  # 创建一个空向量保存人类同源基因
  human_genes <- character(length(mouse_genes))
  
  # 遍历每个鼠基因,查找对应的人类同源基因
  for (i in seq_along(mouse_genes)) {
    mouse_gene <- mouse_genes[i]
    matching_rows <- which(gene_mapping$mouse_gene == mouse_gene)
    
    # 如果找到一个一对一的对应关系,保存人类基因名
    if (length(matching_rows) == 1) {
      human_genes[i] <- gene_mapping$human_gene[matching_rows]
    } else {
      human_genes[i] <- NA # 否则标记为 NA
    }
  }
  
  return(human_genes)
}




Mus_deg <- read.csv("Age_Neu_preAnno_deg.csv",header=TRUE,row.names=1)
Mus_deg_sig <- subset(Mus_deg,logfoldchanges < 0 & pvals_adj < 0.05)



library(dplyr)
Mus_top100_2 <- list()

for(i in unique(Mus_deg_sig$group)){
	genes <- tail(Mus_deg_sig$names[Mus_deg_sig$group == i],300)
	
	my_vector <- convert_mouse_to_human(genes)
	
	cleaned_vector <- subset(my_vector, !is.na(my_vector))
	
	name <- paste("Mice_",i,sep="")
	Mus_top100_2[[name]] <- cleaned_vector

}
lengths(Mus_top100_2)
#head(Mus_top100_2[["Mice_Neu_1"]])



Homo_deg <- read.csv("/home/luchun/scRNA/Bone_neu/04_cell/neu/Neu_preAnno_deg.csv",header=TRUE,row.names=1)
Homo_deg_sig <- subset(Homo_deg,logfoldchanges < 0 & pvals_adj < 0.05)
head(Homo_deg_sig)


Homo_top100_2 <- list()

for(i in unique(Homo_deg_sig$group)){
	genes <- tail(Homo_deg_sig$names[Homo_deg_sig$group == i],300)

	name <- paste("Human_",i,sep="")
	Homo_top100_2[[name]] <- genes

}
lengths(Homo_top100_2)
#head(Homo_top100_2[["Human_Neu_1"]])


shared_genes <- lapply(Homo_top100_2, function(x) {
   lapply(Mus_top100_2, function(y) {
     length(intersect(x, y))
   })
})


shared_genes


tmp <- c()
for(i in 1:length(Homo_top100_2)){	
	for(j in 1:length(Mus_top100_2)){
		inter <- intersect(Homo_top100_2[[i]],Mus_top100_2[[j]])
        tmp <- c(tmp,length(inter))
	}
}


#这里的数字要更改，列数是人的neu类别数
tmp <- matrix(tmp,ncol=7)

colnames(tmp) <- names(Homo_top100_2)
rownames(tmp) <- names(Mus_top100_2)

#colnames(tmp) <- paste("Human_",1:9,sep="")
#rownames(tmp) <- paste("Mice_",1:7,sep="")

tmp



# 假设这是你的除数向量


divisors <- lengths(Mus_top100_2)


# 对每一行进行除法操作
for (i in 1:nrow(tmp)) {
  tmp[i, ] <- tmp[i, ] / divisors[i]*100
}


tmp


'''

for (i in 1:nrow(tmp)) {

	for (j in 1:ncol(tmp)) {
		tmp[i,j] <- tmp[i,j ] / (lengths(Mus_top100_2)[i]*0.5+lengths(Homo_top100_2)[j]*0.5)*100
	}
}


tmp
'''


library(pheatmap)
library(RColorBrewer)
mycol <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(100)

pdf("DOWN_similarity.pdf",width=9)
pheatmap(tmp,
		cluster_rows = FALSE,
		cluster_cols = FALSE,
		border_color = "white",
		display_numbers = TRUE,
		fontsize_number = 12, 
		color=mycol,
		cellwidth =40, cellheight=40,
		number_format = "%.2f",
		number_color = "black",
		angle_col=45,
		main="Down regulated genes similarity")
dev.off()








