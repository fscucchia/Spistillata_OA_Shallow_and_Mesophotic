
files=c(
'Home/STAR-spis-DE1/deep-shallow/deep7.6-shallow7.6.merged-DE-FPM.lite.txt',
'Home/STAR-spis-DE1/deep-shallow/deep7.8-shallow7.8.merged-DE-FPM.lite.txt',
'Home/STAR-spis-DE1/deep-shallow/deep8.2-shallow8.2.merged-DE-FPM.lite.txt'
) #DESeq output tables

files_description=c('comparison','comparison','comparison')

outidx=c(
'Home/STAR-spis-DE1/deep-shallow/enrichment/7_6',
'Home/STAR-spis-DE1/deep-shallow/enrichment/7_8',
'Home/STAR-spis-DE1/deep-shallow/enrichment/8_2'
)

f_selected_genes=''

# common column
column_id='ID'

# coumn per test
regexpr_logfc=c('.*\\.logFC','.*\\.log2FoldChange')
regexpr_fdr=c('FDR','padj')
regexpr_cpm=c('^CPM\\.','^FPM\\.')
#regexpr_condition=c('^[A-Za-z0-9_]+\\.[A-Za-z0-9_]+')
regexpr_condition=c('[A-Za-z0-9_\\.]+(?=\\.)') # perl positive lookahead

##########################################################################

isDE = function(x,logfc,direction){
	de=FALSE
	logfc1=x[1]
	fdr1=x[2]
	if((!is.na(logfc1))&(!is.na(fdr1))){
		if(as.numeric(fdr1)<0.05){
			if(direction == 'regulated'){
				de=TRUE
			}else if(direction=='upregulated'){
				if((as.numeric(logfc1) >= logfc)&(as.numeric(logfc1) >= 0)){de=TRUE}
			}else if(direction=='downregulated'){
				if((as.numeric(logfc1) <= logfc)&(as.numeric(logfc1) <= 0)){de=TRUE}
			}
		}
	}
	return(de)
}
	
select_background = function(x){
	select1=FALSE
	if((!is.na(x[1]))&(!is.na(x[2]))){select1=TRUE}
	return(select1)
}

selectGenes=F
if(nchar(f_selected_genes)>0){
	selectGenes=T
	t_selectedGenes=read.csv(f_selected_genes,sep="\t",stringsAsFactors=F,header=F)
	selectedGenes = t_selectedGenes[,1]
}

for(f in 1:length(files)){
	f1=files[f]
	
	cutoff_logfc=c(1,-1,0,0,0)
	cutoff_direction=c('upregulated','downregulated','upregulated','downregulated','regulated')
	
	t1=read.csv(f1,sep="\t",stringsAsFactors=F,header=T)
	t1 = t1[(!is.na(t1[,column_id]))&(t1[,column_id]!=""),]
	
	if(selectGenes){
		select1=sapply(t1[,column_id],function(x) x %in% selectedGenes)
		t1=t1[select1,]
	}
	
	if(nrow(t1)>0){
		column_logfc = names(t1)[sapply(names(t1),function(x) grepl(paste(regexpr_logfc,collapse='|'),x,perl=T))]
		column_fdr   = names(t1)[sapply(names(t1),function(x) grepl(paste(regexpr_fdr,collapse='|'),x,perl=T))]
		condition    = sapply(column_logfc,function(x) regmatches(x,regexpr(paste(regexpr_condition,collapse='|'),x,perl=T)))
		columns_CPM = names(t1)[sapply(names(t1),function(x) grepl(paste(regexpr_cpm,collapse='|'),x,perl=T))]
		
		j=0
		for(i in 1:length(column_logfc)){
			cat(files_description[f]," ",column_logfc[i],"\n")
			
			cpms=t1[,columns_CPM]
			mean_cpm=apply(cpms,1,mean)
			
			for(k in 1:length(cutoff_logfc)){
				select = apply(t1[,c(column_logfc[i],column_fdr[i])],1,FUN=isDE,logfc=cutoff_logfc[k],direction=cutoff_direction[k]) # select DE genes as the foreground list
				
				# selection of a background list which includes all expressed genes :
				background1 = apply(t1[,c(column_logfc[i],column_fdr[i])],1,FUN=select_background) 
				
				df=data.frame(ID=t1[,column_id],
							  condit=paste0(files_description[f],' ',condition[i],' ',cutoff_direction[k],ifelse(cutoff_logfc[k]!=0,paste0('.Log2FC',cutoff_logfc[k]),'')),
							  length=mean_cpm,bkgrnd=background1,select=select)
				j=j+1
				names(df)[2:length(names(df))]=sapply(names(df)[2:length(names(df))],function(x) paste0(x,'.',j))
				if(j>1){df1=cbind(df1,df[,2:ncol(df)])
				}else{df1=df}
			}
		}
		
		write.table(df1,file=paste0(outidx[f],'.',files_description[f],'.txt'),sep="\t",quote=F,row.names=F)
		
		# statistics:
		genes_DEs        = df1[names(df1)[sapply(names(df1),function(x) grepl('select\\.\\d+',x,perl=T))]]
		genes_background = df1[names(df1)[sapply(names(df1),function(x) grepl('bkgrnd\\.\\d+',x,perl=T))]]
		comparison       = df1[names(df1)[sapply(names(df1),function(x) grepl('condit\\.\\d+',x,perl=T))]]
		
		statistics1 = data.frame( comparison = apply(comparison,2,unique)  , genes_DEs = apply(genes_DEs,2,sum)  , genes_background = apply(genes_background,2,sum) )
		
		write.table(statistics1,file=paste0(outidx[f],'.',files_description[f],'.statistics.txt'),sep="\t",quote=F,row.names=F)
	}else{
		cat("no genes selected\n")
		break
	}
	
}

