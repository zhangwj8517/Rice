 
library(pheatmap)   #加载pheatmap 包；
data1=read.table("nomalize_counts_k4_deng_and_mine.txt",header=T,row.names=1,sep="\t")
data2=read.table("nomalize_counts_k27_deng_and_mine.txt",header=T,row.names=1,sep="\t")

matrix=cor(data1[,1:8])   #计算相关系数；
####write.table(matrix,"coefficient_matrix_overlap_nomalize_k27_ROOT_reads_all_peak_counts.txt",sep="\t")             #将相关系数计算结果输出存储到你的电脑里，存储为1个txt文件；
pdf("nomalize_counts_k4_deng_and_mine.pdf")
pheatmap(matrix,cluster_rows=F,cluster_cols=F,display_numbers=T) # 行和列都不聚类，并且在热图中显示数值；
dev.off()



pdf("aaa.pdf")
pheatmap(matrix,cluster_rows=F,cluster_cols=F,display_numbers=T,fontsize_number = 18,   face_number="bold") 
dev.off()


pdf("aaa.pdf")
pheatmap(matrix,cluster_rows=F,cluster_cols=F,display_numbers=T,fontsize_number = 18) 
dev.off()
