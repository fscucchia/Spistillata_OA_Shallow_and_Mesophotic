
params_file=''

#######################################################################################

# for command line Rscript run
if(nchar(params_file)==0){
	# get loop start - end (to allow parallel runs ...)
	args = commandArgs(trailingOnly=T)
	if(length(args)!=8){ 
		stopifnot(F)
	}else{
		args0 = commandArgs(F)
		parseArgs = function(x) strsplit(sub("^--", "", x), "=")
		args1 = as.data.frame(do.call("rbind", parseArgs(args0)))
		args2 = as.list(as.character(args1$V2))
		names(args2) = args1$V1
		check1=function(list1,k){if(k %in% names(list1)){return(list1[[k]])}else{cat(paste0('\'',k,'\' doesn\'t exit\n')); q("no",1,F)}}
		
		start1      = as.integer(check1(args2,'start'))
		end1        = as.integer(check1(args2,'end'))
		outDir1     = as.character(check1(args2,'outDir'))
		f1          = as.character(check1(args2,'DE_file'))
		f_genesToGo = as.character(check1(args2,'geneToGo_file'))
		title0      = as.character(check1(args2,'title'))
		f1_id1      = as.character(check1(args2,'column_id'))
		splitID     = as.character(check1(args2,'splitID'))
		if(splitID=='T'){splitID=TRUE}else{splitID=FALSE}
	}
	select_range=c(start1:end1)
}else{
	source(params_file)
	select_range=c(start1:end1)
}

stopifnot(file.exists(f1))
stopifnot(file.exists(f_genesToGo))
stopifnot(dir.exists(outDir1))

overide1=FALSE

#######################################################################################
library(goseq)

ttt_clustOrDE = read.csv(f1,sep="\t",stringsAsFactors=F,header=T)
ttt_clustOrDE = ttt_clustOrDE[ (!is.na(ttt_clustOrDE[,f1_id1 ]))&(ttt_clustOrDE[,f1_id1 ] != ""),]
comparisons0=names(ttt_clustOrDE)[grepl("^select\\..*",names(ttt_clustOrDE),perl=T)]
comparisons1=sapply(comparisons0,function(xx) gsub("^select\\.","",xx,perl=T))
comparisons0
comparisons1

genesToGo=read.csv(f_genesToGo,sep="\t",stringsAsFactors=F,header=T)
gene2cat = genesToGo[,c('ID','term_id')]
nrow(gene2cat)
term_id2name = genesToGo[!duplicated(genesToGo$term_id),c('term_id','term_name')]

k=1
test=F
for(i in 1:length(select_range)){
	cols_select=comparisons0[sapply(comparisons0,function(xx) grepl(paste0('^select\\.',select_range[i],'(_\\d+)*$'),xx,perl=T))]
	for(j in 1:length(cols_select)){
		col_select=cols_select[j]
		comparison1b=gsub('^select\\.|_\\d+$','',col_select,perl=T) # this line deals with cases where the selection group, e.g. select.1, has several subgroups, e.g. select.1_1,select.1_2,select.1_3 ... , and all these subgroups have the same background/length/condit columns
		col_length=paste0('length.',comparison1b)
		col_bkgrnd=paste0('bkgrnd.',comparison1b)
		col_condit=paste0('condit.',comparison1b)
		cat('col_select=',col_select,"\n")
		cat(col_length,"\n")
		cat(col_bkgrnd,"\n")
		cat(col_condit,"\n")
		cond=ttt_clustOrDE[ttt_clustOrDE[,col_bkgrnd]==T,][1,col_condit]
		cat('cond=',cond,"\n")
		title1=paste0(title0,'_',i,'.',cond,'.',col_select)
		
		outDir=paste0(outDir1,'/',cond)
		if((dir.exists(outDir))&(!overide1)){next}
		dir.create(outDir, showWarnings = FALSE) # in the first iteration a direcrory 'outDir' will be creared and then be used in its next iterations
		setwd(outDir)
		
		t_DE     = ttt_clustOrDE[ttt_clustOrDE[,col_bkgrnd]==T,c(f1_id1,col_select,col_length)]
		t_DE$score=1
		f_de=paste0(outDir,'/',title1,'.DE.txt')
		write.table(t_DE,file=f_de,sep="\t",row.names=F,quote=F)
		
		save.image(paste0(outDir,'/',title1,'.1.Rdata'))
		
		if(nrow(t_DE) > 25){
			
			# entichment: ####################
			
			genes_bias = as.integer(t_DE[,col_length]) #+ 4000 # TODO - this is only temporary - delete the random var here
			genes=sapply(t_DE[,col_select],function(xx) ifelse(xx,1,0))
			names(genes)=t_DE[,f1_id1]
			
			pdf(paste0(outDir,'/',title1,'.goseq.pdf'))
			pwf = nullp(DEgenes = genes, genome =  NULL,        id = NULL,     bias.data = genes_bias, plot.fit=T)	
			dev.off()
			
			GO.wall=   goseq(pwf, genome =  NULL, id = t_DE[,f1_id1]  ,gene2cat=gene2cat)
			GO.nobias= goseq(pwf, genome =  NULL, id = t_DE[,f1_id1], ,gene2cat=gene2cat, method="Hypergeometric")
			#GO.samp=   goseq(pwf, genome =  GO.db.build, id = GO.db.db, method="Sampling",repcnt=1000)
			
			goseqDf=function(gseq,maxAdjPval){
				gseq$adjPval = p.adjust(gseq$over_represented_pvalue, method="BH")
				g = gseq[gseq$adjPval<maxAdjPval,]
				return(g)
			}
			
			GO.wall.1=   goseqDf(GO.wall,0.15)
			GO.nobias.1= goseqDf(GO.nobias,0.15)
			
			# add gene ids to table: ######################
			
			DEgenes = names(genes)[genes == 1]
			
			DEgenesInCat.wall = c()
			for(z in 1:nrow(GO.wall.1)){
				cat1 = GO.wall.1[z,'category']
				DEgenesInCat.wall[z] = paste(DEgenes[DEgenes %in% gene2cat[gene2cat$term_id == cat1,'ID']],collapse='; ')
			}
			GO.wall.2 = GO.wall.1
			GO.wall.2$DEgenes = DEgenesInCat.wall
			
			DEgenesInCat.nobias = c()
			for(z in 1:nrow(GO.nobias.1)){
				cat1 = GO.nobias.1[z,'category']
				DEgenesInCat.nobias[z] = paste(DEgenes[DEgenes %in% gene2cat[gene2cat$term_id == cat1,'ID']],collapse='; ')
			}
			GO.nobias.2 = GO.nobias.1
			GO.nobias.2$DEgenes = DEgenesInCat.nobias
			
			#gosPerGenes    = getgo(names(genes),NULL, NULL)
			
			# add term names of save tables: #######################
			
			if('term' %in% names(GO.wall)){
				write.table(GO.wall.2,file=paste0(outDir,'/',title1,'.goseq_wall.txt'),sep="\t",row.names=F,quote=F)
				write.table(GO.nobias.2,file=paste0(outDir,'/',title1,'.goseq_nobias.txt'),sep="\t",row.names=F,quote=F)
			}else{
				GO.wall.3   = merge(GO.wall.2,  term_id2name,by.x='category',by.y='term_id',all.x=T,all.y=F)
				GO.nobias.3 = merge(GO.nobias.2,term_id2name,by.x='category',by.y='term_id',all.x=T,all.y=F)
				write.table(GO.wall.3,file=paste0(outDir,'/',title1,'.goseq_wall.txt'),sep="\t",row.names=F,quote=F)
				write.table(GO.nobias.3,file=paste0(outDir,'/',title1,'.goseq_nobias.txt'),sep="\t",row.names=F,quote=F)
			}
		}else{
			sink(paste0(outDir,'/',title1,'.log.txt'))
			cat("not enough data\n")
			sink()
		}
		
		save.image(paste0(outDir,'/',title1,'.2.Rdata'))
	}
}

file.create(paste0(outDir,'/OK'))