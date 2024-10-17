library(dplyr)
path=getwd()
setwd(path)
folders=dir(pattern="_R")
files_path=list()
for(j in 1:length(folders)){
	files_path[j]=paste(path,"/",folders[j],sep="")
}
C=list()
for(i in 1:(length(files_path))){
	assign(paste("TPM",i,sep=""),
	       read.table(file = paste0(files_path[[i]],"/",folders[i],"_abund.tab"),header = T,sep = '\t')%>%
		       select(,one_of("Gene.ID","TPM"))%>%
		       group_by(Gene.ID)%>%summarise(TPM=sum(TPM)))
	C[[i]]=get(paste("TPM",i,sep=""))
}
PRJ=Reduce(function(x,y) dplyr::left_join(x,y,by="Gene.ID"),C)%>%as.data.frame()
rownames(PRJ)=PRJ[,1]
PRJ=PRJ[,-1]
names(PRJ)=c(folders)
setwd(path)
write.table(PRJ,file="mergeTPM.tsv")
