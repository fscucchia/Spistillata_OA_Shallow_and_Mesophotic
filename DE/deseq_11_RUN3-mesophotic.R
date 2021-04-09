################################################
# data:

RData1='/Home/STAR-spis-DE1/spis1_params.RData'
load(RData1)

file_filter_in= ''
file_filter_out=''
fout2='/Home/STAR-spis-DE1/mesophotic'
workingDir1='/Home/STAR-spis-DE1/mesophotic'
include_factor_values=c('Mesophotic') #'Shallow')
include_factor='depth'
f_design='/Home/design/design2.csv'
f_design_sep=","
design_columns = c('sample','pH','depth','colony')
include_gene_regexpr="^gene\\d+|^SpisGene\\d+"
min1=5
filtering_option=1

# loaded parameters of setRawCountsTable_1.R output .RData file:
# fout1,title,file_gene_len,xx,x,m2,factors

vars = ls(all.names = TRUE)

### DE params ######################################

# following setting relevant for a specific baseline set automatically based on alphabetic order. It is possible to change baseline using relevel() or defined the levels of the factor.
# compare_type ={'name','contrast','contrast_vector'}
# when compare_type = 'name', then for example: "compare='banana_vs_orange'", or any other value in resultsNames(dds).
# when compare_type = 'contrast', and compare=vector, which is a numeric vector of the same size as resultsNames(dds), for example c(0,0,0,-1,0,0,0).
# when compare_type = 'contrast', and compare is a list of: a list and a vector, then for example: "compare=list(list('apple_vs_orange','banana_vs_orange'),c(1,-1))", sould find "apple_vs_banana". 
# Moreover, such contranst can be more complex, for example: contrast=list and listValues=vector, for example: contrast=list("apple_vs_orange", c("potato_vs_orange","tomato_vs_orange")), listValues=c(1,-1/2) .
# In such case we compare the effect of "apple" to the average of some other effects, "porato" and "tomato".
# when compare type = LRT (Likelihood Ratio Test), compare=NA.

compare_type=list()
compare=list()
compare_title=list()

compare_type[[1]]="name";		compare[[1]]="pH_7.6_vs_8.2";				compare_title[[1]]="pH_7.6_vs_8.2"
compare_type[[2]]="name";		compare[[2]]="pH_7.8_vs_8.2" ;				compare_title[[2]]="pH_7.8_vs_8.2"
compare_type[[3]]="contrast";	compare[[3]]=list(list("pH_7.6_vs_8.2","pH_7.8_vs_8.2"),c(1,-1));	compare_title[[3]]="pH_7.6_vs_7.8"

compare_title_levels=c("pH_7.6_vs_8.2","pH_7.8_vs_8.2","pH_7.6_vs_7.8")
compare_title_levels_groups=c(1,1,1,2,3,3)

relevel_column=c('pH')
relevel_factor=list('8.2') # control

# ordering and filtering #################################

design_t=read.csv(f_design,sep=f_design_sep,stringsAsFactors=F,header=T)
design_t = design_t[,design_columns]
for(i in 1:ncol(design_t)){design_t[,i]=sapply(design_t[,i],function(x) gsub("\\s","",x,perl=TRUE))}
names(xx) == titles # must all be TRUE
idx1=sapply(titles,function(x) match(x,design_t$sample))
design_t1 = design_t[idx1,]
design_t$sample == titles # must all be TRUE

# removal of unwanted samples:

keep_samples = design_t1[,include_factor] %in% include_factor_values
counts0 = xx[,keep_samples]
design_t2 = design_t1[keep_samples,]

names(xx)[!keep_samples] # check which samples were excluded
names(xx)[keep_samples]  # check which samples were included
names(counts0) == design_t2$sample # all must be TRUE
titles=design_t2$sample

# factors naming:

sum(design_t2$sample == titles)
sum(design_t2$sample != titles)

#############################################

counts1_=counts0
if(nchar(include_gene_regexpr)>0){
	toInclude=sapply(rownames(counts1_),function(x) grepl(include_gene_regexpr,x,perl=T))
	counts1_ = counts1_[toInclude,]
	sink(file=paste0(fout2,'.log.txt'),append = FALSE)
	cat("including",sum(toInclude),'genes\n')
	sink()
}
if(nchar(file_filter_out)>0){
	filterOut=read.csv(file_filter_out,sep="\t",stringsAsFactors=F,header=T)
	nrow(filterOut)
	toRemove=(rownames(counts1_) %in% rownames(filterOut))
	sink(file=paste0(fout2,'.log.txt'),append = FALSE)
	cat("filtering",sum(toRemove),'genes\n')
	sink()
	counts1_ = counts1_[!toRemove,]
}
if(nchar(file_filter_in)>0){
	filterIn=read.csv(file_filter_in,sep="\t",stringsAsFactors=F,header=T)
	nrow(filterIn)
	toInclude=(rownames(counts1_) %in% rownames(filterIn))
	sink(file=paste0(fout2,'.log.txt'),append = TRUE)
	cat("selecting",sum(toInclude),'genes\n')
	sink()
	counts1_ = counts1_[toInclude,]
}
sink(file=paste0(fout2,'.log.txt'),append = TRUE)
if((nchar(file_filter_out)>0) || (nchar(file_filter_in)>0)){
	cat("gene count after further filtering unwanted taxa:\n")
	print(nrow(counts1_))
	cat('read counts before taxa filtering:\n')
	print(apply(counts0,2,sum))
	cat('read counts after taxa filtering:\n')
	print(apply(counts1_,2,sum))
}
sink()

# filtering low count genes:

if(filtering_option==1){
	# option 1:
	keep = rowSums(counts1_>min1) >= 2
}else if(filtering_option==2){
	# option 2:
	library(edgeR)
	keep = rowSums(cpm(counts1_)>min1) >= 2
}else if(filtering_option==3){
	# option 3:
	rs = rowSums(counts1_)
	theta = 0.2
	keep = (rs > quantile(rs, probs=theta))
}
counts1 = counts1_[keep,]

sink(file=paste0(fout2,'.log.txt'),append = TRUE)
cat("genes with higher than cutoff read counts:")
print(table(keep))
cat("original gene count = ",nrow(counts0),"\n")
cat("final gene count    = ",nrow(counts1),"\n")
cat("final gene count = \n")
print(apply(counts1,2,sum))
sink()

write.table(counts1,file=paste0(fout2,'.filtered-counts.txt'),row.names=T,col.names=T,sep='\t',quote=FALSE)
save.image(file = paste0(fout2,".0.RData"))

#####################################################################################
# DE

library('DESeq2')
library('ggplot2')
library('vegan')
library('plotly')
library('reshape2')

sum(design_t2$sample == titles)
sum(design_t2$sample != titles)
#colData=as.data.frame(apply(design_t2,c(1,2),as.factor))
colData=design_t2
rownames(colData)=titles

for(i in 1:length(relevel_column)){
	colData[,relevel_column[i]] = factor(colData[,relevel_column[i]])
	colData[,relevel_column[i]] = relevel(colData[,relevel_column[i]],ref=relevel_factor[[i]])
}

write.table(colData,file=paste0(fout2,'.design.txt'),row.names=T,col.names=T,sep='\t',quote=FALSE)

modelMatrixTest=as.data.frame(model.matrix( ~ pH,data=colData))
dds_ = DESeqDataSetFromMatrix(countData= counts1,colData= colData,design = ~ pH)
write.table(modelMatrixTest,file=paste0(fout2,'.model-matrix.txt'),row.names=T,col.names=T,sep='\t',quote=FALSE)

dds_
dds_ = dds_[rowSums(counts(dds_)) > 1, ]
dds_

dds = DESeq(dds_)
sink(file=paste0(fout2,'.comparisons.txt'),append = FALSE)
resultsNames(dds)
sink()

# DE results:
results1=list()
results0=list()
for(i in 1:length(compare)){
	
	title1=compare_title[[i]]
	cat(title1,"\n")
	
	if(compare_type[[i]] == 'name'){
        cat("name\n")
		res = results(dds,name=compare[[i]]) # from resultsNames()
	}else if(compare_type[[i]] == 'contrast'){
		if((is.vector(compare[[i]]))&(!is.list(compare[[i]]))){
			# when compare_type = 'contrast', and compare=vector, which is a numeric vector of the same size as resultsNames(dds), for example c(0,0,0,-1,0,0,0).
			res = results(dds,contrast=compare[[i]])
		}else if((length(compare[[i]])==2)&(is.list(compare[[i]][1]))){
            cat("contrast 2 elements\n")
			# when compare_type = 'contrast', and compare is a list of a list and a vector, then for example: "compare=list(list('apple_vs_orange','banana_vs_orange'),c(1,-1))", sould find "apple_vs_banana". 
			# Moreover, such contranst can be more complex, for example: contrast=list and listValues=vector, for example: contrast=list("apple_vs_orange", c("potato_vs_orange","tomato_vs_orange")), listValues=c(1,-1/2) .
			res = results(dds,contrast=compare[[i]][[1]],listValues=compare[[i]][[2]])
		}else{
			cat("unsupported contrast structure\n")
			break
		}
    }else if(compare_type[[i]] == 'LRT'){
		cat("LRT\n")
		dds_lrt = DESeq(dds_, test="LRT", reduced = ~ 1)
		res = results(dds_lrt)
	}
    
	pdf(paste0(fout2,'.',title1,'.DE.MA.pdf'))
	plotMA(res,main=title1,ylim=c(-5,5))
	dev.off()
	
	res$ID = row.names(res)
	res1=res[!is.na(res$padj),]
	write.table(res1,file=paste(fout2,title1,'DE','txt',sep='.'),sep="\t",row.names=F,quote=F)
	results1[[i]]=res
		
	if(i==1){append1=F}else{append1=T}
	sink(paste0(fout2,'.pairs.DE-summary.txt'), append=append1, split=F)
	cat("contrast: ",title1,"\n")
	cat("# genes with valid DE results = ",nrow(res1),"\n")
	print(summary(res,alpha=0.05))
	cat("------------------\n")
	sink()
}

fpm1=fpm(dds)
fpm2=data.frame(fpm1)
names(fpm2)=sapply(names(fpm2),function(counts1) paste('FPM',counts1,sep='.'))
fpm2$ID=row.names(fpm2)
write.table(fpm2,file=paste(fout2,'FPM','txt',sep='.'),sep="\t",row.names=F,quote=F)

if(nchar(file_gene_len)>0){
	t_length0=read.csv(file_gene_len,sep="\t",stringsAsFactors=F,header=T) ##
	t_length=t_length0[,c('gene','length')] ##
	fpkm1=fpm(dds)
	fpkm2=data.frame(fpkm1)
	names(fpkm2)=sapply(names(fpkm2),function(counts1) paste('FPKM',counts1,sep='.'))
	fpkm2$ID=row.names(fpkm2)
	nrow(fpkm2)
	#fpkm3=merge(fpkm2,t_length,by.x='ID',by.y='Symbol',all.x=T,all.y=F) ##
	fpkm3=merge(fpkm2,t_length,by.x='ID',by.y='gene',all.x=T,all.y=F) ##
	nrow(fpkm3)
	fpkm4=sweep(fpkm3[,colnames(fpkm3)!='ID' & colnames(fpkm3)!='length'],1,fpkm3$length,FUN=function(counts1,yy) counts1/(yy/1000.0))
	fpkm4$length=fpkm3$length
	fpkm4$ID=fpkm3$ID
	
	write.table(fpkm4,file=paste(fout2,'FPKM','txt',sep='.'),sep="\t",row.names=F,quote=F)
}

# NMDS:

fpm2_=fpm2[,1:(ncol(fpm2)-1)]
fpm2_=data.frame(apply(fpm2_,c(1,2),function(x) log10(x+1)))

pdf(paste0(fout2,'.NMDS','.pdf'))
dist=c('bray','euclidean')
setwd(workingDir1)
for(d in 1:length(dist)){
	dist1 = dist[d]
	nmds1=metaMDS(t(fpm2_),k=3,dist=dist1,maxit=20000,try=1000)
	df_mds1=as.data.frame(nmds1$points)
	df_mds1$names=rownames(nmds1)
	df_mds1$sample = colData$sample
	df_mds1$pH = colData$pH
	df_mds1$depth = colData$depth
	
	gg1 = ggplot(df_mds1, aes(x=MDS1, y=MDS2, color=factor(pH),shape=factor(colony)))
	gg1 = gg1 + geom_point(size=4,alpha=0.7,stroke=1.5)
	gg1 = gg1 + xlab(paste0("MDS1"))
	gg1 = gg1 + ylab(paste0("MDS2"))
	#gg1 = gg1 + scale_shape_identity()
	gg1 = gg1 + geom_text(aes(label=sample),hjust=2.0, vjust=2.0,size=2,colour='black')
	gg1 = gg1 + ggtitle(paste('dist=',dist1,'; log10 FPM NMDS (Vegan metaMDS) ; stress=',nmds1$stress))
	print(gg1)
	
	#plotly1 = plot_ly(df_mds1, x = ~MDS1, y = ~MDS2, z = ~MDS3, color = ~factor(Treatment),symbol = ~factor(Diurnal),symbols=c('circle','square','x','diamond','o'),text = ~ Day) %>% add_markers() %>% add_text() %>% 
	#	layout(scene = list( aspectratio = list(x = 1,y = 1,z = 1) ))
	#htmlwidgets::saveWidget(as_widget(plotly1),paste0("plot_ly.",dist1,".NMDS_3D.html"),selfcontained =F)
}
dev.off()

# PCA:
#rld = rlog(dds, blind=FALSE)
#vsd = vst(dds, blind=FALSE)
#data = plotPCA(vsd,intgroup=c("factor1"), returnData=TRUE)
#percentVar = round(100 * attr(data, "percentVar"))
#pdf(paste0(fout2,'.PCA','.','pdf'))
#ggplot(data, aes(PC1, PC2, color=factor1)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance"))
#dev.off()

# dispersion:

pdf(paste0(fout2,'.dispersion','.','pdf'))
plotDispEsts(dds)
dev.off()

# PCA:
rld = rlog(dds, blind=FALSE)
#vsd = vst(dds, blind=FALSE)
data = plotPCA(rld,intgroup=c("pH"), returnData=TRUE)
percentVar = round(100 * attr(data, "percentVar"))
pdf(paste0(fout2,'.PCA-rlog','.','pdf'))
ggplot(data, aes(PC1, PC2, color=factor(pH,levels=c('8.2','7.8','7.6')),shape=colony)) + 
	geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
	ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
	geom_text(aes(label=name),hjust=0.5, vjust=-0.0,size=2,colour='black')
dev.off()


# PCA

pca1 = prcomp(fpm2_, scale = F)
pca2=as.data.frame(pca1$rotation)
pca2=pca2[,1:3]
pca2$sample =colData$sample
pca2$pH=colData$pH
pca2$depth=colData$depth
comp = as.data.frame(summary(pca1)$importance)

gg1 = ggplot(pca2, aes(x=PC1, y=PC2, color=factor(pH),shape=factor(colony)))
gg1 = gg1 + geom_point(size=4,alpha=0.7,stroke=1.5)
gg1 = gg1 + xlab(paste0("PC1 var explained=",round(100*comp['Proportion of Variance','PC1'],digits=1),"%"))
gg1 = gg1 + ylab(paste0("PC2 var explained=",round(100*comp['Proportion of Variance','PC2'],digits=1),"%"))
gg1 = gg1 + geom_text(aes(label=sample),hjust=0.0, vjust=0.0,size=2,colour='black')
gg1 = gg1 + ggtitle(paste('log10FPM PCA'))
pdf(paste0(fout2,'.PCA-logFPM','.pdf'))
gg1
dev.off()


###############################################################################
# arrangements

for(i in 1:length(results1)){
	t0=data.frame(results1[[i]])
	t0$FC=sapply(t0$log2FoldChange,function(counts1) ifelse(counts1>=0,2^abs(counts1),-1*(2^abs(counts1))))
	names(t0)=sapply(names(t0),function(counts1) ifelse(counts1=='ID','ID',paste(compare_title[[i]],counts1,sep='.')) )
	#cat(names(t0),"\n")
	if(i==1){
		t1=t0
	}else{
		t1=merge(t1,t0,by='ID',all.x=T,all.y=T)
	}
	cat(nrow(t1),"\n")
}
nrow(t0)
nrow(t1)
m1_=merge(t1,fpm2,by='ID',all.x=T,all.y=T)
nrow(m1_)

if(nchar(file_gene_len)>0){
	m2_=m1_
	m2_=merge(m1_,fpkm4,by='ID',all.x=T,all.y=F)
}else{m2_=m1_}

write.table(m2_,file=paste(fout2,'merged-DE-FPM','txt',sep='.'),sep="\t",row.names=F,quote=F)
DE_table1 = m2_[,!grepl("\\.pvalue$|\\.baseMean$|\\.stat$|\\.lfcSE$|\\.FC$",names(m2_),perl=T)]
write.table(DE_table1,file=paste(fout2,'merged-DE-FPM.lite','txt',sep='.'),sep="\t",row.names=F,quote=F)

#############
# CPM frequencies:

molten1=melt(fpm2,id='ID')
molten1$log10FPM = sapply(molten1$value,function(x) log10(x+0.0001))
molten1$sample = sapply(molten1$variable,function(x) sub('FPM\\.X?','',x,perl=T))
#molten1$sample = sapply(molten1$variable,function(x) sub('FPM\\.','',x,perl=T))
molten2=merge(molten1,design_t1,by.x='sample',by.y='sample',all.x=T,all.y=F)
nrow(molten1)
nrow(molten2)
pdf(paste0(fout2,'.histogram','.pdf'))
gg2=ggplot(molten2,aes(x=log10FPM,colour=as.factor(sample)))
gg2 = gg2 + geom_density(fill=NA) + theme_bw() +xlab('log10 normalized counts')
gg2
dev.off()

###############
# DE summary:

cols_padj=names(m2_)[grepl("\\.padj$",names(m2_),perl=T)]
cols_logfc=names(m2_)[grepl("\\.log2FoldChange$",names(m2_),perl=T)]
for(i in 1:length(cols_padj)){
    DE_class=apply(m2_[,c(cols_padj[i],cols_logfc[i])],1,function(x) ifelse(is.na(x[1]),"NA",ifelse(x[1]>0.05,"ns",ifelse(x[2]>0,'up','down'))))
    DE_class0=as.data.frame(table(DE_class))
    DE_class0$comparison=sub("\\.padj$","",cols_padj[i])
    if(i==1){DE_class1=DE_class0}else{DE_class1=rbind(DE_class1,DE_class0)}
}
if(length(compare_title_levels)>0){
	gg3=ggplot(data=DE_class1)+geom_col(aes(x=factor(comparison,levels=compare_title_levels),y=Freq,fill = factor(DE_class,levels=c('up','down','ns','NA')))) +
		theme(axis.text.x = element_text(angle = 45, hjust = 1,color=factor(compare_title_levels_groups))) + guides(fill=guide_legend(title="type")) + xlab("") +
		scale_fill_manual(values=c('brown3','chartreuse4','cornflowerblue','gray'))
}else{
	gg3=ggplot(data=DE_class1)+geom_col(aes(x=factor(comparison),y=Freq,fill = factor(DE_class))) +
		theme(axis.text.x = element_text(angle = 45, hjust = 1)) + guides(fill=guide_legend(title="type")) + xlab("") +
		scale_fill_manual(values=c('brown3','chartreuse4','cornflowerblue','gray'))
}
pdf(paste0(fout2,'.DE-stat','.pdf'),height=7,width=12)
gg3
dev.off()

write.table(DE_class1,sep="\t",quote=F,file=paste0(fout2,'.DE-stat','.txt'),row.names=F)

################

save.image(file = paste0(fout2,".RData"))
save(file = paste0(fout2,".lite.RData"),list=c('fpm2','results1','DE_table1','m2_','colData'))
