f1='/data/home/STAR-spis-SNPs/gatk1-joint/pH2-joint-all.txt'
out1='/data/home/STAR-spis-SNPs/gatk1-joint/pH2-joint-all.2'
out2='/data/home/STAR-spis-SNPs/gatk1-joint/pH2-joint-all.3'

remove1='.5.filtered1.DP6.snps.pass.2.txt'

step2_samplesLociCoverage="/data/home/STAR-spis-SNPs/gatk1-joint/mosDepth/merge.txt"

f_design='/data/home/design/design2.csv'
design_columns = c('sample','pH','depth','colony')
f_design_sep=','

selectedSamples=list()
selectedSamples[[1]]=c(
'A1','A4','A7','B1','B4','B7','C3','C4','C7',
'VI_2_1','VI_2_4','VI_2_7','V_2_4','V_2_7','V2_1',
'X_2_1','X_2_4','X_2_8'
)
selectedSamples[[2]]=c(
'A1','A4','A7','B1','B4','B7','C3','C4','C7'
)
selectedSamples[[3]]=c(
'VI_2_1','VI_2_4','VI_2_7','V_2_4','V_2_7','V2_1',
'X_2_1','X_2_4','X_2_8'
)

selectedSamples_snpset=list()
selectedSamples_snpset[[1]]=c(
'A1','A4','A7','B1','B4','B7','C3','C4','C7',
'VI_2_1','VI_2_4','VI_2_7','V_2_4','V_2_7','V2_1',
'X_2_1','X_2_4','X_2_8'
)
selectedSamples_snpset[[2]]=c(
'A1','A4','A7','B1','B4','B7','C3','C4','C7'
)
selectedSamples_snpset[[3]]=c(
'VI_2_1','VI_2_4','VI_2_7','V_2_4','V_2_7','V2_1',
'X_2_1','X_2_4','X_2_8'
)

selected_title0=c('all','shallow','mesophotic')

f_expression = c(NA)
expression_regFPKM = '^FPM\\.|^FPKM\\.'
expression_ID = 'ID'

######### step 1 ##################################################
# find all loci of SNPs found in at least one of the samples. These loci will be saved in a table which will be used for finding RNA-Seq converage for each loci and sample (including samples with no SNP call at the loci)

library(data.table)
library(gridExtra)
library(grid)

t1=fread(f1,sep="\t",header=T,stringsAsFactors=F)
names(t1) = sapply(names(t1),function(x) gsub(remove1,"",x))
nrow(t1)

cols = names(t1)[grepl("CHROM_POS",names(t1))]
chrPos=apply(t1[,cols,with=FALSE],1,function(x) x[!is.na(x)][1])
t1$chrPos = chrPos
sum(duplicated(chrPos))

sumNA=function(x){ sum(!is.na(x))}
t1[,sumNA := sumNA(.SD),by=chrPos,.SDcols=cols] # select speciifc cols (.SDcols), group by specific key (chrPos, which is unique for all rows), add a column sumNA
table(t1$sumNA)

for(i in 1:length(cols)){
	cat(cols[i],sumNA(t1[,cols[i],with=FALSE]),"\n")
}

# TRANSITION (1 for a bi-allelic transition (SNP), 0 for bi-allelic transversion (SNP), -1 for INDELs and multi-allelics) 
#table(t1$TRANSITION.SH60D[!is.na(t1$TRANSITION.SH60D)])

g1 = t1[,lapply(.SD,function(x)ifelse(is.na(x),0,1)),by=chrPos,.SDcols=cols]

t1[,chr := sub("_\\d+$","",chrPos,perl=T),by=chrPos]
t1[,pos := as.integer(sub("^.*_","",chrPos,perl=T)),by=chrPos]
nrow(t1)

t1[,pos1 := as.integer(pos+1),by=chrPos]
bed1=t1[,c('chr','pos','pos1','chrPos'),with=FALSE]
bedFile1 = paste0(out1,".bed")
fwrite(bed1,quote=F,sep="\t",row.names=F,col.names=F,file=bedFile1)

save.image(file=paste0(out1,"step1.RData"))

######### step 2 ##################################################
# find the converage of each SNP loci tested here with the script bams-mos-depth_1.sh

######## step 3 ##############################################
# merge the SNP data with the coverage data obtained in the previous step

library(data.table)
library(ggplot2)

t2 = read.csv(step2_samplesLociCoverage,sep=" ",stringsAsFactors=F,header=T)
t2_ids = t2[,grepl("^ID\\..*",names(t2),perl=T)]
t2_cov = t2[,grepl("^coverage\\..*",names(t2),perl=T)]
ids = apply(t2_ids,1,function(x) ifelse(length(unique(x))==1,x[1],NA))
sum(is.na(ids)) # must be zero
length(ids)
t2_cov$minCov = apply(t2_cov,1,min)
t2_cov$ID = ids

n1=names(t1)[grepl('^REF\\.',names(t1),perl=T)]
refs1 = as.data.frame(t1[,n1,with=F])
ref1 = apply(refs1,1,function(x){x=as.character(x); u=unique(x[!is.na(x)]); return(ifelse(length(u)==1,u[1],NA))})
t1$REF = ref1
sum(is.na(ref1)) # must be zero

n2=names(t1)[grepl('^CHROM_POS\\.',names(t1),perl=T)]
chrpos1 = as.data.frame(t1[,n2,with=F])
chrp1 = apply(chrpos1,1,function(x){x=as.character(x); u=unique(x[!is.na(x)]); return(ifelse(length(u)==1,u[1],NA))})
sum(is.na(chrp1)) # must be zero

cols_t1=names(t1)[grepl("^chrPos$|^chr$|^pos$|^REF$|^ALT\\.|^DP\\.",names(t1),perl=T)]
t1b = t1[,cols_t1,with=F]

nrow(t1)
nrow(t1b) # must be equal to t1
nrow(t2_cov) # must be equal to t1
m1=merge(as.data.frame(t1b),as.data.frame(t2_cov),by.x='chrPos',by.y='ID',all.x=T,all.y=T)
nrow(m1) # must be equal to t1

write.table(m1,paste0(out1,"snps1.txt"),sep="\t",quote=F,row.names=F)

alt0 =m1[,names(m1)[grepl("^ALT\\.",names(m1),perl=T)]]
rownames(alt0) = m1$chrPos
alt1=apply(alt0,c(1,2),function(x) ifelse(is.na(x),0,ifelse(grepl("\\w+,\\w+",x,perl=T),2,1)))

save.image(file=paste0(out1,"step3.RData"))

### SNPs analysis ##############

library(data.table)
library(gridExtra)
library(grid)
library(gdsfmt)
library(SNPRelate)
library(reshape2)
library(ggplot2)
library(vegan)
library(gridExtra)
library(grid)
library(plyr)

nrow(alt1)
nrow(m1)
info1=m1[,c('chr','pos','REF')]
rownames(info1)=m1$chrPos

ids_depth=list()
d=30
ids_depth[[d]]=m1[m1$minCov>d,'chrPos']
length(ids_depth[[d]])

alt2=alt1[rownames(alt1) %in% ids_depth[[d]],]
info2=info1[rownames(info1) %in% ids_depth[[d]],]
z=sapply(colnames(alt2),function(x) sub('^ALT\\.','',x,perl=T))
colnames(alt2)=z

design_t=read.csv(f_design,sep=f_design_sep,stringsAsFactors=F,header=T)
design_t = design_t[,design_columns]
rownames(design_t)=design_t$sample

for(i in 1:length(selectedSamples)){
	selected_title = selected_title0[i]
	cat(selected_title,"\n")
	selected_samples = selectedSamples[[i]]
	selected_samples_snpset = selectedSamples_snpset[[i]] # the Snps set test is separated for cases where the same individuals are used more than once
	
	alt3=alt2[,selected_samples ]
	alt3_snpSet=alt2[,selected_samples_snpset ] # the Snps set test is separated for cases where the same individuals are used more than once
	info3=info2
	
	# get snp set
	gds1_snpSet=paste0(out2,".",selected_title,".depth",d,".temp.gds")
	snpgdsCreateGeno(gds1_snpSet, genmat = as.matrix(alt3_snpSet),
		sample.id = colnames(alt3_snpSet), snp.id = rownames(alt3_snpSet),
		snp.chromosome = info3$chr,
		snp.position = info3$pos,
	   snpfirstdim=TRUE)
	genofile_snpSet = snpgdsOpen(gds1_snpSet)
	genofile_snpSet
	set.seed(1000)
	snpset = snpgdsLDpruning(genofile_snpSet, ld.threshold=0.5,autosome.only=FALSE) # Try different LD thresholds for sensitivity analysis
	snpset.id = unlist(unname(snpset))
	snpgdsClose(genofile_snpSet)
	
	gds1=paste0(out2,".",selected_title,".depth",d,".gds")
	snpgdsCreateGeno(gds1, genmat = as.matrix(alt3),
		sample.id = colnames(alt3), snp.id = rownames(alt3),
		snp.chromosome = info3$chr,
		snp.position = info3$pos,
	   snpfirstdim=TRUE)
	genofile = snpgdsOpen(gds1)
	genofile
	
	write.table(alt3,file=paste0(gds1,'.txt'),sep="\t",row.names=T,quote=F)
	
	#  IBS heatmap guidelines:
	u=length(selectedSamples[[i]])/2
	rect1 = matrix(data=rep(1:u,times=u),nrow=u,ncol=u,byrow=TRUE)
	rect1[lower.tri(rect1, diag = FALSE)]=NA
	rect2=melt(rect1)
	rect2 = rect2[!is.na(rect2$value),]
	rect2$xmin=rect2$Var1*2 - 1.5
	rect2$xmax=rect2$Var1*2 + 0.5
	rect2$ymin=rect2$Var2*2 - 1.5
	rect2$ymax=rect2$Var2*2 + 0.5
	
	# IBS heatmap:
	IdentityByState = snpgdsIBS(genofile, num.thread=1,autosome.only=FALSE)
	ibs0 = IdentityByState$ibs
	ibs0[upper.tri(ibs0, diag = FALSE)]=NA
	ibs1=as.data.frame(ibs0)
	names(ibs1)=IdentityByState$sample.id
	ibs1$ID = IdentityByState$sample.id
	ibs2=melt(ibs1,id.vars='ID')
	ibs3=merge(ibs2,design_t[,c('sample','pH','colony')],by.x="ID",by.y="sample",all.x=T,all.y=F)
	ibs4=merge(ibs3,design_t[,c('sample','pH','colony')],by.x="variable",by.y="sample",all.x=T,all.y=F)
	ibs4$sample.1=apply(ibs4,1,function(x) paste0(x[4],".",x[5]))
	ibs4$sample.2=apply(ibs4,1,function(x) paste0(x[6],".",x[7]))
	ibs4 = ibs4[!is.na(ibs4$value),]
	nrow(ibs4)
	nrow(IdentityByState$ibs)*ncol(IdentityByState$ibs)
	cols0 = c('salmon','salmon4','palegreen4','purple','lightskyblue','pink','mediumblue','navajowhite4','orangered3')
	labelCol = mapvalues(design_t[ibs1$ID,'colony'],from=unique(design_t[ibs1$ID,'colony']),
						to= cols0[1:length(unique(design_t[ibs1$ID,'colony']))])
	ibs4_ = ggplot() +
			geom_raster(data = ibs4, mapping = aes(x =factor(variable,levels=selected_samples), y = factor(ID,levels=selected_samples), fill = value)) + theme_bw() + ggtitle("SNPs Indentity By State (IBS)") +
			geom_rect(data = rect2,mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),color='black',alpha=0) +
			geom_text(data = ibs4,aes(x =variable, y = factor(ID,levels=selected_samples),label = round(value, 2)),color='red') +
			theme(axis.text.x = element_text(angle = 90, hjust = 1,colour=labelCol),axis.text.y = element_text(colour=labelCol),panel.grid=element_blank(),legend.position="left",axis.title.x=element_blank(),axis.title.y=element_blank()) + coord_fixed()
	pdf(paste0(gds1,'.IdentityByState.pdf'))
	print(ibs4_)
	dev.off()
	
	if(!is.na(f_expression[i])){
		expression0 = read.csv(f_expression[i],sep="\t",header=T,stringsAsFactors=F)
		expression1 = expression0[,grepl(expression_regFPKM,names(expression0),perl=T)]
		rownames(expression1) =  expression0[,expression_ID]
		names(expression1) = sapply(names(expression1),function(x) sub(expression_regFPKM,'',x,perl=T))
		match1 = match(colnames(alt3),colnames(expression1))
		print(sum(is.na(match1)))
		print(match1)
		expression2 = expression1[,colnames(expression1)[match1]]
		
		testCor=apply(expression2,c(1,2),function(x) if(x!=0){return(log10(x))}else{return(-6)})
		#testCor=expression2
		cor0 = cor(testCor,method='spearman')
		cor0[upper.tri(cor0, diag = FALSE)]=NA
		cor0=as.data.frame(cor0)
		cor0$ID=rownames(cor0)
		cor1=melt(cor0,id.vars='ID')
		cor1=cor1[!is.na(cor1$value),]
		cor2=merge(cor1,design_t[,c('sample','depth','colony')],by.x="ID",by.y="sample",all.x=T,all.y=F)
		cor2heatmap = ggplot() +
				geom_raster(data = cor2, mapping = aes(x = variable, y = factor(ID,levels=selected_samples), fill = value)) +  theme_bw()  + geom_rect() +
				geom_rect(data = rect2,mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),color='black',alpha=0) + ggtitle("Spearman correlation of log10 FPKM") +
				theme(axis.text.x = element_text(angle = 90, hjust = 1,colour=labelCol),axis.text.y = element_text(colour=labelCol), panel.grid=element_blank(),legend.position="right",axis.title.x=element_blank(),axis.title.y=element_blank()) + coord_fixed()
		pdf(paste0(gds1,'.expressionCor.pdf'),height=4,width=8)
		print(grid.arrange(ibs4_,cor2heatmap,nrow=1, ncol=8, layout_matrix = rbind(c(1,1,1,1,2,2,2,2))))
		dev.off()
				
		}
	
	snpgdsClose(genofile)
}


save.image(file=paste0(out1,"step4.RData"))
