
indir='/Home/pH_experiment_adults/output/mapping_diamond'
files_pattern='(.*)(?=\\.m8$)'
design_file='/Home/pH_experiment_adults/scripts/design/design2.csv'

design_factor='depth'
design_factor_levels=c("Shallow","Mesophotic")
design_sample='sample'

cols1=c('query' ,'subject',	'identity','alignment_length','mismatches','gap_openings','query_start','query_end','subject_start','subject_end','E_value','bit_score')
clades_pattern=c('^([Ss]ymb|[Cc]ladocopium)(.*)','^(?!([Ss]ymb|[Cc]ladocopium))(.*)')
clades_name=c('Symbiodinium','Metazoa')
outidx='/Home/pH_experiment_adults/output/mapping_diamond/diamond'

#################################

library(data.table)
library(plyr)
library(ggplot2)
library(reshape2)

files0 = list.files(indir)
names0 = sapply(files0,function(x) regmatches(x,regexpr(files_pattern,x,perl=T))[1])
names1=names0[!is.na(names0)]
files1=files0[!is.na(names0)]

t_design=read.csv(design_file,sep=",",stringsAsFactors=F)
t_design1 = t_design[,c(design_sample,design_factor)]
names(t_design1) = c('sample','factor1')
t_design1$factor1 = factor(t_design1$factor1,levels=design_factor_levels)

speciesFreq=list()
for(i in 1:length(files1)){
    cat(names1[i],"\n")
    f1=paste0(indir,"/",files1[i])
    blast1=fread(f1)
    names(blast1)=cols1
    blast1=blast1[,c('query' ,'subject','identity','alignment_length','E_value','bit_score')]
    blast1$i =1:nrow(blast1)
    top1=function(x){return(1:length(x))}
    blast1[, priority := top1(i),by=query]
    blast1 = blast1[priority ==1,]
    blast1[,species := regmatches(subject,regexpr('^[A-Za-z0-9]+',subject,perl=TRUE))]
    
    cat("freq:\n")
    speciesFreq[[i]]=as.data.frame(table(blast1$species))
    speciesFreq[[i]]=speciesFreq[[i]][order(-speciesFreq[[i]]$Freq),]
    names(speciesFreq[[i]])=c('species',paste0(names1[i]))
    }

save.image(file=paste0(outidx,'.1.RData'))

for(i in 1:length(speciesFreq)){
    if(i==1){
        speciesFreq1=speciesFreq[[i]]
    }else{
        speciesFreq1=merge(speciesFreq1,speciesFreq[[i]],by='species',all.x=TRUE,all.y=TRUE)
    }
}
speciesFreq2=speciesFreq1[order(-apply(speciesFreq1[,2:ncol(speciesFreq1)],1,sum)),]
samples1 = names(speciesFreq2)[2:ncol(speciesFreq2)]
factors1 = mapvalues(samples1,from=c(t_design1[,1]),to=c(as.character(t_design1[,2])))
titles1 = sapply(1:length(factors1),function(x) paste0(factors1[x],'.',samples1[x]))

whichClade=function(x){
	y='other'
	for(w in 1:length(clades_name)){
		z=regmatches(x,regexpr(clades_pattern[w],x,perl=T))
		if(length(z)>0){y=clades_name[w]; break}
	}
	return(y)
}
samples1=names(speciesFreq2)[2:ncol(speciesFreq2)]
speciesFreq2$clade = sapply(speciesFreq2$species,whichClade)
speciesFreq2 = speciesFreq2[,c('species','clade',samples1)]
write.table(speciesFreq2,sep="\t",file=paste0(outidx,".species-freq.txt"),row.names=F,quote=F)

speciesFreq_byClade=list()
speciesFreq_byClade_per=list()
for(k in 1:length(clades_name)){
    speciesFreq_byClade[[k]]=speciesFreq2[speciesFreq2$clade==clades_name[k],]
    speciesFreq_byClade_per[[k]] =sweep(speciesFreq_byClade[[k]][,samples1],2,apply(speciesFreq_byClade[[k]][,samples1],2,sum),"/")
    speciesFreq_byClade_per[[k]]$species=speciesFreq_byClade[[k]]$species
    speciesFreq_byClade_per[[k]] = speciesFreq_byClade_per[[k]][,c('species',samples1)]
    write.table(speciesFreq_byClade_per[[k]],sep="\t",file=paste0(outidx,".species-freq.",clades_name[k],".txt"),row.names=F,quote=F)
    
    mlt1=melt(speciesFreq_byClade_per[[k]],id.vars='species')
    mlt2=merge(mlt1,t_design1,by.x='variable',by.y='sample',all.x=T,all.y=T)
    pdf(paste0(outidx,".species-freq.",clades_name[k],".pdf"),width=10)
    hm1 = ggplot(data = mlt2, mapping = aes(x = variable, y = species, fill = value)) +
        geom_raster() + ggtitle('proportion of mapped reads') +
        theme(axis.text.x=element_text(angle = 45)) + facet_grid(cols = vars(factor1),scale='free')
    print(hm1)
    dev.off()
    
    mlt1b=melt(speciesFreq_byClade[[k]],id.vars=c('species','clade'))
    mlt2b=merge(mlt1b,t_design1,by.x='variable',by.y='sample',all.x=T,all.y=T)
    #mlt2b$value=sapply(mlt2b$value,function(x) log10(x+1))
    pdf(paste0(outidx,".species count.",clades_name[k],".pdf"),width=10)
    hm2 = ggplot(data = mlt2b, mapping = aes(x = variable, y = species, fill = value)) +
        geom_raster() + ggtitle('count of mapped reads') +
        theme(axis.text.x=element_text(angle = 45)) + facet_grid(cols = vars(factor1),scale='free') +
        scale_fill_gradientn(colours = c('white','gray','black'),limits=c(q[1],q[2]),oob = scales::squish)
    print(hm2)
    dev.off()
}
