rm(list=ls())
WORK_SPACE="E:/Zhenling_论文数据计算/2016_Rice/20171022_rice_data/1.原始数据文件reads_counts/"
setwd(WORK_SPACE)
library(gsubfn)
library(sqldf)
library(pheatmap)
library(tcltk)

HOMEHETER=c("home","heter")
N99N=c("N9","9N")
K4K27=c("k4","k27")
ROOTLEAF=c("root","leaf")
WTJ=c("WTJ7","WTJ12","WTJ14","WTJ24","WTJ34","WTJ37")

SOURCE1="1.原始reads的snpsplit1016_counts分布/"
SOURCE2="2.marker_edgeR/"
STEP1_PATH="STEP1/"
STEP2_PATH="STEP2/"


SOURCE_K27_LEAF_9N_DATA =read.table("snpslit_uniq_macs14_k27_leaf1022_9N.txt",head=T,fill=T,sep = "\t" )
SOURCE_K27_LEAF_N9_DATA =read.table("snpslit_uniq_macs14_k27_leaf1022_N9.txt",head=T,fill=T,sep = "\t" )
SOURCE_K27_ROOT_9N_DATA =read.table("snpslit_uniq_macs14_k27_root1022_9N.txt",head=T,fill=T,sep = "\t" )
SOURCE_K27_ROOT_N9_DATA =read.table("snpslit_uniq_macs14_k27_root1022_N9.txt",head=T,fill=T,sep = "\t" )
SOURCE_K4_LEAF_9N_DATA  =read.table("snpslit_uniq_macs14_k4_leaf1022_9N.txt",head=T,fill=T,sep = "\t" )
SOURCE_K4_LEAF_N9_DATA  =read.table("snpslit_uniq_macs14_k4_leaf1022_N9.txt",head=T,fill=T,sep = "\t" )
SOURCE_K4_ROOT_9N_DATA  =read.table("snpslit_uniq_macs14_k4_root1022_9N.txt",head=T,fill=T,sep = "\t" )
SOURCE_K4_ROOT_N9_DATA  =read.table("snpslit_uniq_macs14_k4_root1022_N9.txt",head=T,fill=T,sep = "\t" )



uniq_macs_k27_leaf_DATA =read.table("counts_merged_sorted_uniq_macs14_k27_leaf20171022.txt",head=T,fill=T,sep = "\t" )
uniq_macs_k27_root_DATA =read.table("counts_merged_sorted_uniq_macs14_k27_root20171022.txt",head=T,fill=T,sep = "\t" )
uniq_macs_k4_leaf_DATA =read.table("counts_merged_sorted_uniq_macs14_k4_leaf20171022.txt",head=T,fill=T,sep = "\t" )
uniq_macs_k4_root_DATA =read.table("counts_merged_sorted_uniq_macs14_k4_root20171022.txt",head=T,fill=T,sep = "\t" )



SOURCE_FILTER_K27_LEAF_9N_DATA=sqldf(" select * from SOURCE_K27_LEAF_9N_DATA where mix_nip + mix_9311 + x9n_nip + x9n_9311   >=6 ")
SOURCE_FILTER_K27_LEAF_N9_DATA=sqldf(" select * from SOURCE_K27_LEAF_N9_DATA where mix_nip + mix_9311 + n9_nip  + n9_9311   >=6 ")
SOURCE_FILTER_K27_ROOT_9N_DATA=sqldf(" select * from SOURCE_K27_ROOT_9N_DATA where mix_nip + mix_9311 + x9n_nip + x9n_9311   >=6 ")
SOURCE_FILTER_K27_ROOT_N9_DATA=sqldf(" select * from SOURCE_K27_ROOT_N9_DATA where mix_nip + mix_9311 + n9_nip  + n9_9311   >=6 ")
SOURCE_FILTER_K4_LEAF_9N_DATA =sqldf(" select * from SOURCE_K4_LEAF_9N_DATA  where mix_nip + mix_9311 + x9n_nip + x9n_9311   >=6 ")
SOURCE_FILTER_K4_LEAF_N9_DATA =sqldf(" select * from SOURCE_K4_LEAF_N9_DATA  where mix_nip + mix_9311 + n9_nip  + n9_9311   >=6 ")
SOURCE_FILTER_K4_ROOT_9N_DATA =sqldf(" select * from SOURCE_K4_ROOT_9N_DATA  where mix_nip + mix_9311 + x9n_nip + x9n_9311   >=6 ")
SOURCE_FILTER_K4_ROOT_N9_DATA =sqldf(" select * from SOURCE_K4_ROOT_N9_DATA  where mix_nip + mix_9311 + n9_nip  + n9_9311   >=6 ")

nrow(SOURCE_FILTER_K27_LEAF_9N_DATA);
nrow(SOURCE_FILTER_K27_LEAF_N9_DATA);
nrow(SOURCE_FILTER_K27_ROOT_9N_DATA);
nrow(SOURCE_FILTER_K27_ROOT_N9_DATA);
nrow(SOURCE_FILTER_K4_LEAF_9N_DATA );
nrow(SOURCE_FILTER_K4_LEAF_N9_DATA );
nrow(SOURCE_FILTER_K4_ROOT_9N_DATA );
nrow(SOURCE_FILTER_K4_ROOT_N9_DATA );

###　uniq_macs_*_DATA join (SOURCE_FILTER_*_9N_DATA union SOURCE_FILTER_*_N9_DATA)

filter_uniq_macs_k27_leaf_DATA=sqldf(" select A.*  from uniq_macs_k27_leaf_DATA A,  (select peak from SOURCE_FILTER_K27_LEAF_9N_DATA union select peak from SOURCE_FILTER_K27_LEAF_N9_DATA) B where a.peak = b.peak" )
filter_uniq_macs_k27_root_DATA=sqldf(" select A.*  from uniq_macs_k27_root_DATA A,  (select peak from SOURCE_FILTER_K27_ROOT_9N_DATA union select peak from SOURCE_FILTER_K27_ROOT_N9_DATA) B where a.peak = b.peak" )
filter_uniq_macs_k4_leaf_DATA =sqldf(" select A.*  from uniq_macs_k4_leaf_DATA  A,  (select peak from SOURCE_FILTER_K4_LEAF_9N_DATA  union select peak from SOURCE_FILTER_K4_LEAF_N9_DATA)  B where a.peak = b.peak" )
filter_uniq_macs_k4_root_DATA =sqldf(" select A.*  from uniq_macs_k4_root_DATA  A,  (select peak from SOURCE_FILTER_K4_ROOT_9N_DATA  union select peak from SOURCE_FILTER_K4_ROOT_N9_DATA)  B where a.peak = b.peak" )




nrow(filter_uniq_macs_k27_leaf_DATA)
nrow(filter_uniq_macs_k27_root_DATA)
nrow(filter_uniq_macs_k4_leaf_DATA )
nrow(filter_uniq_macs_k4_root_DATA )

 

write.table(filter_uniq_macs_k27_leaf_DATA, file =paste0(STEP1_PATH, "union_mix_k27_leaf.txt") ,  append = FALSE, quote = F,sep = "\t",row.names = F,col.names = TRUE, fileEncoding = "UTF-8" )
write.table(filter_uniq_macs_k27_root_DATA, file =paste0(STEP1_PATH, "union_mix_k27_root.txt") ,  append = FALSE, quote = F,sep = "\t",row.names = F,col.names = TRUE, fileEncoding = "UTF-8" )
write.table(filter_uniq_macs_k4_leaf_DATA,  file =paste0(STEP1_PATH, "union_mix_k4_leaf.txt") ,  append = FALSE, quote = F,sep = "\t",row.names = F,col.names = TRUE, fileEncoding = "UTF-8" )
write.table(filter_uniq_macs_k4_root_DATA,  file =paste0(STEP1_PATH, "union_mix_k4_root.txt") ,  append = FALSE, quote = F,sep = "\t",row.names = F,col.names = TRUE, fileEncoding = "UTF-8" )
	 
	 
write.table(SOURCE_FILTER_K27_LEAF_9N_DATA, file =paste0(STEP1_PATH, "filter_K27_LEAF_9N_count.txt") ,  append = FALSE, quote = F,sep = "\t",row.names = F,col.names = TRUE, fileEncoding = "UTF-8" )
write.table(SOURCE_FILTER_K27_LEAF_N9_DATA, file =paste0(STEP1_PATH, "filter_K27_LEAF_N9_count.txt") ,  append = FALSE, quote = F,sep = "\t",row.names = F,col.names = TRUE, fileEncoding = "UTF-8" )
write.table(SOURCE_FILTER_K27_ROOT_9N_DATA, file =paste0(STEP1_PATH, "filter_K27_ROOT_9N_count.txt") ,  append = FALSE, quote = F,sep = "\t",row.names = F,col.names = TRUE, fileEncoding = "UTF-8" )
write.table(SOURCE_FILTER_K27_ROOT_N9_DATA, file =paste0(STEP1_PATH, "filter_K27_ROOT_N9_count.txt") ,  append = FALSE, quote = F,sep = "\t",row.names = F,col.names = TRUE, fileEncoding = "UTF-8" )
write.table(SOURCE_FILTER_K4_LEAF_9N_DATA , file =paste0(STEP1_PATH, "filter_K4_LEAF_9N_count.txt") ,  append = FALSE, quote = F,sep = "\t",row.names = F,col.names = TRUE, fileEncoding = "UTF-8" )
write.table(SOURCE_FILTER_K4_LEAF_N9_DATA , file =paste0(STEP1_PATH, "filter_K4_LEAF_N9_count.txt") ,  append = FALSE, quote = F,sep = "\t",row.names = F,col.names = TRUE, fileEncoding = "UTF-8" )
write.table(SOURCE_FILTER_K4_ROOT_9N_DATA , file =paste0(STEP1_PATH, "filter_K4_ROOT_9N_count.txt") ,  append = FALSE, quote = F,sep = "\t",row.names = F,col.names = TRUE, fileEncoding = "UTF-8" )
write.table(SOURCE_FILTER_K4_ROOT_N9_DATA , file =paste0(STEP1_PATH, "filter_K4_ROOT_N9_count.txt") ,  append = FALSE, quote = F,sep = "\t",row.names = F,col.names = TRUE, fileEncoding = "UTF-8" )


######### STEP2 : calc chisq test

fun.chisq<-function(x1,x2)  
{ 
chisq.test(c(x1,x2))$p.value    
}


SOURCE_FILTER_DATA=list(
SOURCE_FILTER_K27_LEAF_9N_DATA,
SOURCE_FILTER_K27_LEAF_N9_DATA,
SOURCE_FILTER_K27_ROOT_9N_DATA,
SOURCE_FILTER_K27_ROOT_N9_DATA,
SOURCE_FILTER_K4_LEAF_9N_DATA ,
SOURCE_FILTER_K4_LEAF_N9_DATA ,
SOURCE_FILTER_K4_ROOT_9N_DATA ,
SOURCE_FILTER_K4_ROOT_N9_DATA )

### mapping SOURCE_FILTER_DATA order 
K4K27_LEAFROOT_N99N=c(
"K27_LEAF_9N",
"K27_LEAF_N9",
"K27_ROOT_9N",
"K27_ROOT_N9",
"K4_LEAF_9N",
"K4_LEAF_N9",
"K4_ROOT_9N",
"K4_ROOT_N9"
)


for (Psfd in 1:length(SOURCE_FILTER_DATA) )
{

    CIRCLE_DATA =as.data.frame(SOURCE_FILTER_DATA[Psfd])
	OUT_NAME = K4K27_LEAFROOT_N99N[Psfd]
	colnames(CIRCLE_DATA) = c("peak","mix_ref","mix_alt","X9N_ref","X9N_alt")
    CIRCLE_MOD_DATA=sqldf(" 
	select 
		peak,
		case when mix_ref  = 0 then 1 else mix_ref  end as mix_ref,
		case when mix_alt  = 0 then 1 else mix_alt  end as mix_alt,
		case when X9N_ref  = 0 then 1 else X9N_ref  end as X9N_ref,
		case when X9N_alt  = 0 then 1 else X9N_alt  end as X9N_alt
	from CIRCLE_DATA ")

	mapply(fun.chisq,CIRCLE_MOD_DATA$mix_ref,CIRCLE_MOD_DATA$mix_alt)->CIRCLE_MOD_DATA$mix_p_value
    p.adjust(CIRCLE_MOD_DATA$mix_p_value,method="fdr")->CIRCLE_MOD_DATA$mix_Q_value 

	mapply(fun.chisq,CIRCLE_MOD_DATA$X9N_ref,CIRCLE_MOD_DATA$X9N_alt)->CIRCLE_MOD_DATA$X9N_p_value
	p.adjust(CIRCLE_MOD_DATA$X9N_p_value,method="fdr")->CIRCLE_MOD_DATA$X9N_Q_value 

	
	P_Q_VALUE_DATA=sqldf("
		select 
			CASE 
				WHEN mix_Q_value > 0.05 THEN 3
				WHEN mix_Q_value <= 0.05 AND mix_ref > mix_alt THEN 2
				WHEN mix_Q_value <= 0.05 AND mix_ref < mix_alt THEN 4
				ELSE 9
			END  mix_TYPE_zwj01
			,CASE 
				WHEN X9N_Q_value > 0.05 THEN 3
				WHEN X9N_Q_value <= 0.05 AND X9N_ref > X9N_alt THEN 2
				WHEN X9N_Q_value <= 0.05 AND X9N_ref < X9N_alt THEN 4
				ELSE 9
			END  x9N_TYPE_zwj01
			,peak
			,mix_ref
			,mix_alt
			,X9N_ref
			,X9N_alt
			,mix_p_value
			,mix_Q_value
			,CASE 
				WHEN mix_Q_value > 0.05 THEN 3
				WHEN mix_Q_value <= 0.05 AND mix_ref > mix_alt THEN 2
				WHEN mix_Q_value <= 0.05 AND mix_ref < mix_alt THEN 4
				ELSE 9
			END  mix_TYPE_zwj
			,X9N_p_value
			,X9N_Q_value
			,CASE 
				WHEN X9N_Q_value > 0.05 THEN 3
				WHEN X9N_Q_value <= 0.05 AND X9N_ref > X9N_alt THEN 2
				WHEN X9N_Q_value <= 0.05 AND X9N_ref < X9N_alt THEN 4
				ELSE 9
			END  x9N_TYPE_zwj
		from CIRCLE_MOD_DATA
	
	")
	 
	 mix_TYPE_zwj01_cnt=  sqldf(paste0("
	 select mix_TYPE_zwj01  , count(*) as cnt 
	 from P_Q_VALUE_DATA 
	 group by mix_TYPE_zwj01
	 "))
	 
	 
	 x9N_TYPE_zwj01_cnt =sqldf(paste0("
	 select x9N_TYPE_zwj01  , count(*) as cnt 
	 from P_Q_VALUE_DATA 
	 group by x9N_TYPE_zwj01
	 "))
	 
	write.table(P_Q_VALUE_DATA ,     file =paste0(STEP2_PATH,"Chisq_",OUT_NAME, ".txt") ,      append = FALSE, quote = F,sep = "\t",row.names = F,col.names = TRUE, fileEncoding = "UTF-8" )
	write.table(mix_TYPE_zwj01_cnt , file =paste0(STEP2_PATH,"static_Chisq_",OUT_NAME, "_mix_TYPE.txt") ,  append = FALSE, quote = F,sep = "\t",row.names = F,col.names = TRUE, fileEncoding = "UTF-8" )
	write.table(x9N_TYPE_zwj01_cnt , file =paste0(STEP2_PATH,"static_Chisq_",OUT_NAME, "_x9N_TYPE.txt") ,   append = FALSE, quote = F,sep = "\t",row.names = F,col.names = TRUE, fileEncoding = "UTF-8" )
 

}
 

