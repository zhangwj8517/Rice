rm(list=ls())
memory.limit(60000)  ## requre  60G memory for run 
WORK_SPACE="D:/wheat_bam/methylation_compare/result_2017_11_01"
setwd(WORK_SPACE)
library(methylKit)   


CONTEXT=c("CG","CHG","CHH")
WTJ=c("WTJ7","WTJ12","WTJ14","WTJ24","WTJ34","WTJ37")
PARENT=c("9311","NIP","N9","9N")

TARGET_PATH="F:/methylation_methylkit_compare_20171102/"


Calculate <- function( StrCONTEXT, StrWTJ,StrPARENT  )
{
INFILE_CONTEXT_WTJ = paste0(StrCONTEXT,"_",StrWTJ,".txt")
INFILE_CONTEXT_PARENT = paste0(StrCONTEXT,"_",StrPARENT,".txt")
LIST_CONTEXT_WTJ  = paste0( StrCONTEXT,"_",StrWTJ)
LIST_CONTEXT_PARENT= paste0( StrCONTEXT,"_",StrPARENT)

OUTFILE_DMC_UNITE  =paste0(TARGET_PATH,"DMC_UNITE_" ,StrCONTEXT,"_",StrWTJ,"_",StrPARENT,"_methylated.txt")   ##  "unite_DMC_CHH_WTJ7_Nip_mydiff.txt"
OUTFILE_DMC_myDiff =paste0(TARGET_PATH,"DMC_myDiff_",StrCONTEXT,"_",StrWTJ,"_",StrPARENT,"_methylated.txt")   ##  "DMC_CHH_WTJ7_Nip_mydiff.txt"
OUTFILE_DMC_hyper  =paste0(TARGET_PATH,"DMC_hyper_" ,StrCONTEXT,"_",StrWTJ,"_",StrPARENT,"_methylated.txt")   ##  "DMC_CHH_WTJ7_Nip_hyper_methylated.txt"
OUTFILE_DMC_hypo   =paste0(TARGET_PATH,"DMC_hypo_"  ,StrCONTEXT,"_",StrWTJ,"_",StrPARENT,"_methylated.txt")   ##  "DMC_CHH_WTJ7_Nip_hypo_methylated.txt"
OUTFILE_DMC_differ =paste0(TARGET_PATH,"DMC_differ_",StrCONTEXT,"_",StrWTJ,"_",StrPARENT,"_methylated.txt")   ##  "DMC_CHH_WTJ7_Nip_differentialy_methylated.txt"

OUTFILE_DMR_UNITE  =paste0(TARGET_PATH,"DMR_UNITE_" ,StrCONTEXT,"_",StrWTJ,"_",StrPARENT,"_methylated.txt")   ##  "unite_DMR_CHH_WTJ7_Nip_mydiff.txt"
OUTFILE_DMR_myDiff =paste0(TARGET_PATH,"DMR_myDiff_",StrCONTEXT,"_",StrWTJ,"_",StrPARENT,"_methylated.txt")   ##  "DMR_CHH_WTJ7_Nip_mydiff.txt"
OUTFILE_DMR_hyper  =paste0(TARGET_PATH,"DMR_hyper_" ,StrCONTEXT,"_",StrWTJ,"_",StrPARENT,"_methylated.txt")   ##  "DMR_CHH_WTJ7_Nip_hyper_methylated.txt"
OUTFILE_DMR_hypo   =paste0(TARGET_PATH,"DMR_hypo_"  ,StrCONTEXT,"_",StrWTJ,"_",StrPARENT,"_methylated.txt")   ##  "DMR_CHH_WTJ7_Nip_hypo_methylated.txt"
OUTFILE_DMR_differ =paste0(TARGET_PATH,"DMR_differ_",StrCONTEXT,"_",StrWTJ,"_",StrPARENT,"_methylated.txt")   ##  "DMR_CHH_WTJ7_Nip_differentialy_methylated.txt"



file.list=list(INFILE_CONTEXT_WTJ,INFILE_CONTEXT_PARENT)                                                     
myobj=methRead(file.list,sample.id=list(LIST_CONTEXT_WTJ,LIST_CONTEXT_PARENT),assembly="msu7.0",treatment=c(1,0))
#####DMC##############
meth=unite(myobj, destrand=FALSE)
write.table(meth, OUTFILE_DMC_UNITE, sep='\t') 
rm(myobj)    ### release memory 
myDiff=calculateDiffMeth(meth)      
rm(meth)     ### release memory                                
write.table(myDiff, OUTFILE_DMC_myDiff, sep='\t')                                     
myDiff25p.hyper <-getMethylDiff(myDiff,differenc=25,qvalue=0.01,type="hyper")  
myDiff25p.hyper                                                                 
write.table(myDiff25p.hyper,OUTFILE_DMC_hyper,sep='\t')                    
myDiff25p.hypo <-getMethylDiff(myDiff,differenc=25,qvalue=0.01,type="hypo")    
myDiff25p.hypo                                                                  
write.table(myDiff25p.hypo,OUTFILE_DMC_hypo,sep='\t')                      
myDiff25p <-getMethylDiff(myDiff,differenc=25,qvalue=0.01)                     
myDiff25p                                                                       
write.table(myDiff25p,OUTFILE_DMC_differ,sep='\t')                  

### release memory 
rm(myDiff)  
rm(myDiff25p.hyper)
rm(myDiff25p.hypo)
rm(myDiff25p)

#####DMR##############  

file.list=list(INFILE_CONTEXT_WTJ,INFILE_CONTEXT_PARENT)                                                     
myobj=methRead(file.list,sample.id=list(LIST_CONTEXT_WTJ,LIST_CONTEXT_PARENT),assembly="msu7.0",treatment=c(1,0))
tiles <- tileMethylCounts(myobj, win.size = 200,step.size = 200)    
rm(myobj)    ### release memory            
meth_tiles=unite(tiles, destrand=FALSE)  
rm(tiles)    ### release memory    
write.table(meth_tiles, OUTFILE_DMR_UNITE, sep='\t')                                              
myDiff_DMR=calculateDiffMeth(meth_tiles)   
rm(meth_tiles) ### release memory    
write.table(myDiff_DMR, OUTFILE_DMR_myDiff, sep='\t')                    
myDiff25p.hyper_DMR <-getMethylDiff(myDiff_DMR,differenc=25,qvalue=0.01,type="hyper")  
myDiff25p.hyper_DMR                                                                 
write.table(myDiff25p.hyper_DMR,OUTFILE_DMR_hyper,sep='\t')   
myDiff25p.hypo_DMR <-getMethylDiff(myDiff_DMR,differenc=25,qvalue=0.01,type="hypo")    
myDiff25p.hypo_DMR                                                                  
write.table(myDiff25p.hypo_DMR,OUTFILE_DMR_hypo,sep='\t')     
myDiff25p_DMR <-getMethylDiff(myDiff_DMR,differenc=25,qvalue=0.01)                     
myDiff25p_DMR                                                                       
write.table(myDiff25p_DMR,OUTFILE_DMR_differ,sep='\t') 

### release memory
rm(myDiff_DMR)   
rm(myDiff25p.hyper_DMR)
rm(myDiff25p.hypo_DMR)
rm(myDiff25p_DMR)


}




for(Pcont in 1:length(CONTEXT))
{
	for( Pwtj in 1:length(WTJ))
	{
		for (Parent in 1: length(PARENT))
		{
		    ##  NIP, 9311 mapping all WTJ
		    ##  N9 mapping "WTJ7","WTJ12","WTJ14"
			##  9N mapping "WTJ24","WTJ34","WTJ37"
			if (PARENT[Parent] == "NIP" ||  PARENT[Parent] =="9311") 
			{
				print (paste0(CONTEXT[Pcont],"_",WTJ[Pwtj],"_",PARENT[Parent]))
				Calculate(CONTEXT[Pcont],WTJ[Pwtj],PARENT[Parent])
			}
			else if ( (PARENT[Parent] == "N9")   &&    (WTJ[Pwtj]  == "WTJ7" || WTJ[Pwtj]  ==  "WTJ12" || WTJ[Pwtj]  == "WTJ14") )
			{
				print (paste0(CONTEXT[Pcont],"_",WTJ[Pwtj],"_",PARENT[Parent]))
				Calculate(CONTEXT[Pcont],WTJ[Pwtj],PARENT[Parent])
			}
			else if ((PARENT[Parent] == "9N")   &&  (WTJ[Pwtj] == "WTJ24" || WTJ[Pwtj]  == "WTJ34" || WTJ[Pwtj]  == "WTJ37") )
			{
				print (paste0(CONTEXT[Pcont],"_",WTJ[Pwtj],"_",PARENT[Parent]))
				Calculate(CONTEXT[Pcont],WTJ[Pwtj],PARENT[Parent])
			}

		}
	}

}




