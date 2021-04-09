
################################################
# data:

table_type='star' # table_type='hetseq' / table_type='rsem' / table_type='generic'  / table_type='star' / table_type='featureCounts'

f_files_and_factors='$HOME1/design/design2.csv'
f_files_and_factors_sep=","
cols_factors=c('sample','pH','depth')
cols_factors_regexpr_select=c('.*','.*','.*')

dirout1='$HOME1/DE1'

title="spis1"

file_gene_len='$HOME1/stylophoraP-liftover3_concat.gtf.gene-length.txt'
file_gene_len_cols=c('gene','length')

include_gene_regexpr="^gene\\d+|^SpisGene\\d+"

################################################
# populate variables:

fout1=paste0(dirout1,'/',title)

files_and_factors = read.csv(f_files_and_factors,sep=f_files_and_factors_sep,stringsAsFactors=F,header=T) # columns names should look like: sample,file,factor1,factor2,...
select0=list()
for(i in 1:length(cols_factors)){
    select0[[i]]=sapply(files_and_factors[,cols_factors[i]],function(x) grepl(cols_factors_regexpr_select[i],x))
    if(i==1){ select1 = select0[[i]]
    }else{    select1 = select1 & select0[[i]]}
}
files_and_factors1 = files_and_factors[select1,]

files = files_and_factors1[,'file']
factors=list()
for(i in 1:length(cols_factors)){
    factors[[i]] = files_and_factors1[,cols_factors[i]]
}
titles=files_and_factors1[,'sample']

################################################
# read and merge data:

if(length(files)>1){
	t1_=list()
	duplicated1=list()
	for(i in 1:length(files))
	{
		sample1=paste(files[i],sep="")
		cat(sample1,"\n")
		if(table_type=='htseq'){
			t1 = read.csv(file=sample1,sep='\t',header=F)
			t1_[[i]]=t1
			names(t1)=c('gene_id','expected_count')
			t1 = t1[,c('gene_id','expected_count')]
			t1=t1[1:(nrow(t1)-5),]
		}
		else if(table_type=='star'){
			t1 = read.csv(file=sample1,sep='\t',header=F,skip=4,stringsAsFactors=F)
			t1_[[i]]=t1[,1:2]
			names(t1)=c('gene_id','expected_count')
			t1 = t1[,c('gene_id','expected_count')]
		}
		else if(table_type=='featureCounts'){
			t0 = read.csv(file=sample1,sep='\t',header=T,skip=1,stringsAsFactors=F)
			t1=t0
			t1_[[i]]=t1[,c(1,7)]
			cat("names\n")
			print(names(t1_[[i]]))
			names(t1_[[i]])=c('gene_id','expected_count')
			t1 = t1_[[i]]
			
			# get gene Length
			lengths0 = t0[,c('Geneid','Length')]
			names(lengths0) = c('Symbol',as.character(titles[i]))
			if(i==1){
				lengths1_ = lengths0
			}else{
				lengths1_ = merge(lengths1_,lengths0,by='Symbol')
			}
		}
		else{
			t1 = read.csv(file=sample1,sep='\t',header=T)
			t1_[[i]]=t1
			t1 = t1[,c('gene_id','expected_count')]
		}
		names(t1)=c('Symbol',as.character(titles[i]))
		
		##
		duplicated_symbols = t1$Symbol[duplicated(t1$Symbol)]
		t1=t1[!(t1$Symbol %in% duplicated_symbols),]
		duplicated1[[i]]=duplicated_symbols
		##
		
		if(i==1){m=t1
		}else{m=merge(m,t1,by='Symbol',all.x=TRUE,all.y=TRUE)}
	}
	m2=cbind(m[,1],data.frame(apply(m[,2:ncol(m)],c(1,2),as.integer)))
	names(m2)=names(m)

	write.table(m2,file=paste(fout1,'mergedCounts.txt',sep='.'),row.names=T,col.names=T,sep='\t',quote=FALSE)
	
	# check sum of counts:
	
	if(table_type=='htseq'){
		for(i in 1:length(t1_)){
			h1=t1_[[i]][(nrow(t1_[[i]])-4):nrow(t1_[[i]]),]
			names(h1)=c('subset',as.character(titles[i]))
			if(i==1){htseqReport=h1}
			else{htseqReport=merge(htseqReport,h1,by='subset',all.x=TRUE,all.y=TRUE)}
		}
		write.table(htseqReport,file=paste(fout1,'htseqReport.txt',sep='.'),row.names=T,col.names=T,sep='\t',quote=FALSE)
	}
}else{
	m2=read.csv(file=files[1],sep='\t',header=T)
}

if(table_type=='featureCounts'){
	lengths1_$length = apply(lengths1_[,2:ncol(lengths1_)],1,mean)
	t_length = lengths1_[,c('Symbol','length')]
	write.table(t_length,file=paste0(fout1,'meanLength-featureCounts.txt'),sep="\t",quote=F,row.names=F)
	cat(nrow(m2),"\n")
	m2=merge(m2,t_length,by='Symbol',all.x=F,all.y=F)
	cat(nrow(m2),"\n")
}else if(nchar(file_gene_len)>0){
	t_length=read.csv(file=file_gene_len,sep='\t',header=T)
	t_length=t_length[,file_gene_len_cols]
	names(t_length)=c('Symbol','length')
	cat(nrow(m2),"\n")
	m2=merge(m2,t_length,by='Symbol',all.x=F,all.y=F)
	cat(nrow(m2),"\n")
}

select1 = sapply(m2$Symbol,function(x) grepl(include_gene_regexpr,x,perl=T))
sum(select1)
sum(!select1)
m2=m2[select1,]

# counts table for deseq:
x=m2[,2:ncol(m2)]
row.names(x)=m2[,1]
xx=x[,!grepl('^length$',names(x),perl=T)]
nrow(xx)
xx=xx[apply(xx,1,function(xx) sum(is.na(xx))==0),]
nrow(xx)

save.image(file=paste0(fout1,'_params.RData'))
