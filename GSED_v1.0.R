### GSED: Gene set enrichment in disease categories
# Copyright (C) 2016  Dong, Xiao
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Affero General Public License for more details.

# Xiao Dong, 2016.03.16, biosinodx@gmail.com, xiao.dong@einstein.yu.edu
# 2016.06.09 - v1.0 release version of u0.0.1
# 2016.06.09 - u0.0.1 add command line options
# 2016.03.16 - u0.0.0 initial version

# Usage
# Rscript gsed.R workingdir projectname inputgenelist.txt backgroundgenelist.txt repeattime diseaseclass
# workingdir: e.g. ./
# projectname: anyname, e.g. test
# inputgenelist.txt: input gene list, should be tab limited file with header, the first colunmn is the gene names of interest.
# backgroundgenelist.txt: input background gene list, should be tab limited file with header, the first colunmn is the gene names of interest, 2rd colunmn gene start, 3rd colunmn gene end; (I prepared one for protein coding genes, and uploaded to the depository "background_genelist.txt")
# repeattime: recommend 2000
# diseaseclass: disease - gene class; (Simon C Johnson prepared one, and I uploaded to the depository "disease_genename_simon_agingcell_2015.RData")

verion='v1.0'
print(paste("GSED: Gene set enrichment in disease categories, verion", verion))
Args <- commandArgs(TRUE)

# project name
#PN='ps1'
PN=Args[2]
# input gene list, should be tab limited file with header, the first colunmn is the gene names of interest.
#genelist_candidate='./data/ps1_simon.txt'
genelist_candidate=Args[3]
# input background gene list, should be tab limited file with header, the first colunmn is the gene names of interest, 2rd colunmn gene start, 3rd colunmn gene end
#genelist_background='./data/background_genelist.txt'
genelist_background=Args[4]
# input number of random gene set should be generated
#numrandom=2000
numrandom=as.numeric(Args[5])
# load gene class assignment
# gdclass='data/disease_genename_simon_agingcell_2015.RData'
gdclass=Args[6]
# setworking dir
#wd='~/projects/2015-diseaseenrichment/workdir'
wd=Args[1]

### analysis
setwd(wd)
dir.create(PN)
dir.create(paste(PN,'/randomset',sep=''))

genelist_background=read.table(genelist_background, header=T)
genelist_candidate=read.table(genelist_candidate, header=T)

# generate percentiles of gene length in background gene list, may take a few minutes
x=genelist_background
x$l=x[,3]- x[,2]
a=levels(x[,1])
b=vector()
for(i in 1:length(a)){
	b[i]=mean(x$l[x$Associated_Gene_Name==a[i]])
}
blist=data.frame(a,b)
colnames(blist)=c('gene','length')

q=quantile(blist$length, 1:100/100)
x=vector()
for(i in 1:nrow(blist)){
	x[i]=which(blist$length[i]<=q)[1]
}
blist$quantile = x
genelist=genelist_candidate
g = merge(genelist[,1], blist, by=1)

# generate random gene sets which saves to current working dir, this may take an hour for 1000 random gene sets
tmp=as.factor(g$quantile)
tmp=summary(tmp, maxsum=nrow(g))

for(k in 1:numrandom){
randomdat=vector()
for(i in 1:length(tmp)){
	x=blist[names(tmp)[i]==blist$quantile,]
	x=x[sample(nrow(x),tmp[i]),]
	randomdat=rbind(randomdat,x)
}
write.table(randomdat, paste('./',PN,'/randomset/randomgene_set', k, '.txt' , sep=''), col.names=T, row.names=F, sep='\t', quote=F)
}


# 
counting=function(c_genelist,c_categorylist){
	c_out=vector()
	for(i in 1:length(c_categorylist)){
		c_out[i] <- nrow(merge(c_genelist, c_categorylist[[i]], by=1))
	}
	names(c_out)=names(c_categorylist)
	return(c_out)
}

load(gdclass)
datout=vector()
datout=rbind(datout, counting(g, disease_genename))
rownames(datout)='observed'
randomout=vector()
for(i in 1:numrandom){
	x=read.table(paste('./',PN,'/randomset/randomgene_set', i, '.txt' , sep=''), header=T)
	randomout=rbind(randomout, counting(as.vector(x[,1]), disease_genename))
}
rownames(randomout)=paste('randomset',1:numrandom,sep='')
pvalue=vector()
for(i in 1:length(disease_genename)){
	pvalue[i]=(length(which(datout[1,i] <= randomout[,i]))+1)/(numrandom+1)
}

fileout=rbind(datout, pvalue, randomout)
write.table(fileout,paste('./',PN,'/results.txt',sep=''), col.names=T, row.names=T, sep='\t', quote=F)

for(i in 1:length(disease_genename)){
	pdf(file=paste('./',PN,'/DensityPlot_',names(disease_genename[i]),'.pdf', sep=''))
	plot(density(randomout[,i]), main=paste(names(disease_genename[i]),', pvalue=', pvalue[i],sep=''),xlab='# genes (distribution - null hypothesis; red dashed line - observed)', xlim=c(0,max(c(randomout[,i], datout[,i]))))
	abline(v=datout[,i], col='red', cex=100, lty=2)
	dev.off()
}

