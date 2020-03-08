require(writexl)

mlst <- read.delim(snakemake@input[["mlstresults"]],as.is=T,header=F)#read mlst table
cat("red mlst data \n",file=snakemake@log[[1]],sep="")
mlst=mlst[,c(1,3)]#takes 1st and 3rd column (MLSType)
colnames(mlst)=c("ID","MLST_Type")#changes heading
#mlst=read.delim("/home/Share/2018/results/mlst.tab",as.is=T,header=F)
iscfe1 <- read.delim(snakemake@input[["iscfe1_summary"]],as.is=T)#read IS table
cat("red iscfe1 data \n",file=snakemake@log[[1]],sep="",append=T)
colnames(iscfe1)=c("ID","Nr.ISCfe1","Coverage")#changes heading
pcr_list=c(snakemake@input[["pcr"]])#pathrs to all pcr files
cat("red pcr_list data \n",file=snakemake@log[[1]],sep="",append=T)

#pcr_list=c("/home/Share/2018/results//deu-18CS0023/pcr.tab", "/home/Share/2018/results//deu-18CS0371/pcr.tab","/home/Share/2018/results/test3/pcr.tab", "/home/Share/2018/results//deu-> CS0024/pcr.tab")
snp_list=c(snakemake@input[["eleven_snps"]])#pathrs to all snp files
#snp_list=c("/home/Share/2018/results/deu-18CS0023/cfv_snp_markers.tab", "/home/Share/2018/results//deu-18CS0371/cfv_snp_markers.tab","/home/Share/2018/results/test3/cfv_snp_markers.tab", "/home/Share/2018/results//deu-18CS0024/cfv_snp_markers.tab")
cat("red snp_list data \n",file=snakemake@log[[1]],sep="",append=T)

#gets samples IDs from file names
tmp=strsplit(pcr_list,"/")
samleID=rep(NA,length(tmp))
for( i in 1:length(tmp)){
  pos=length(tmp[[i]])-1
  samleID[[i]]=tmp[[i]][[pos]]
}

cat("converted sample IDs \n",file=snakemake@log[[1]],sep="",append=T)


datapcr <- lapply(pcr_list,read.delim,comment.char = "#",header=F,as.is=T)#reads all pcr results
datasnp<- lapply(snp_list,read.delim,comment.char = "#",header=F,as.is=T)#reads all snp results
cat("read lists \n",file=snakemake@log[[1]],sep="",append=T)



for ( i in 1:length(datapcr)){#loops over all pcr and snp results
  pcrtable=datapcr[[i]]
  pcrtable[nchar(pcrtable[,4])==0,4]="neg"#write neg in case of no amplicifcation
  shortpcr=matrix(ncol=nrow(pcrtable),nrow=1)#shorter table just one row per sample
  shortpcr[1,1]=samleID[[i]]#adds sample ID
  shortpcr[1,2:ncol(shortpcr)]=pcrtable[2:nrow(pcrtable),4]#takes column 4 in row 1
  colnames(shortpcr)=c("ID",pcrtable[2:nrow(pcrtable),1])
  datapcr[[i]]=shortpcr#replaces data

  snptable=datasnp[[i]]
  snptable[snptable[,4]==".",4]=snptable[snptable[,4]==".",3]#replaces . with reference nucl
  snpstring=paste0(snptable[,4], collapse = "-")#concatenates bases into string
  shortsnps=matrix(ncol=2,nrow=1)#creates simples table
  colnames(shortsnps)=c("ID","cf_SNP_Markers")#chages heading
  shortsnps[1,1]=samleID[[i]]#add ID
  shortsnps[1,2]=snpstring#add string of SNPS
  datasnp[[i]]=shortsnps
}
cat("cretaed pcr table and snptable \n",file=snakemake@log[[1]],sep="",append=T)

pcrs= do.call("rbind",datapcr)#combines all samples in one table
snps=do.call("rbind",datasnp)#combines all samples in one table

Subspecies=snps#table for final decission
colnames(Subspecies)=c("ID","Subspecies")
Subspecies[,2]="-"#standard is - if decission not knownw
Subspecies=Subspecies[order(Subspecies[,1]),]#reorder row

descission=merge(snps,iscfe1,by.x=1,by.y=1,all.x=T,all.Y=T)#combines data for postion
cat("cretaed decissions \n",file=snakemake@log[[1]],sep="",append=T)

Subspecies[descission["cf_SNP_Markers"]=="T-T-C-T-A-C-A-C-C-A-G" & descission["Nr.ISCfe1"]>0,2]="cfv"#if SNPS are T-T-C-T-A-C-A-C-C-A-G and at leat 1 IS
Subspecies[descission["cf_SNP_Markers"]=="C-C-T-C-G-T-C-T-T-G-A" & descission["Nr.ISCfe1"]==0,2]="cff"#if SNPS are C-C-T-C-G-T-C-T-T-G-Aand at leat 1 IS
  all=merge(Subspecies,snps,by.x=1,by.y=1,all.x=T,all.Y=T)#combines data in table
  all=merge(all,iscfe1,by.x=1,by.y=1,all.x=T,all.Y=T)
  all=merge(all,mlst,by.x=1,by.y=1,all.x=T,all.Y=T)
  all=merge(all,pcrs,by.x=1,by.y=1,all.x=T,all.Y=T)
cat("merged \n",file=snakemake@log[[1]],sep="",append=T)

write_xlsx(all,snakemake@output[["allxls"]])
cat("wrote xls \n",file=snakemake@log[[1]],sep="",append=T)


write.table(all,snakemake@output[["allcsv"]],row.names=F,sep="\t",quote=F)


cat("finnished successfully \n",file=snakemake@log[[1]],sep="",append=T)
