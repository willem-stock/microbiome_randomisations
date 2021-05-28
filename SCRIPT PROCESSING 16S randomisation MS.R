#Willem Stock- randomisation 

library(vegan);library(ggplot2);library(picante); library(bipartite)

#functions####

# c-score function
Cscr<-function(x) {return (C.score(decostand(x,method='pa'), normalise = T))}

# rank variability
rnkvar<-function(x){
  x[x==0]<-NA
  return (apply(apply(x,1,function (y) {return(rank(y,na.last="keep"))}),1,function (z) {return (var(z,na.rm = T))}))
}

# calculate the different indexes for a community dataset
indexes<-function(w,type){
  
  x<-w[,colSums(w)>0]
  faprotax1_on_1_sel<-dplyr::select(faprotax1_on_1,dplyr::one_of(colnames(x)))
  funct_table_dia<-as.matrix(x) %*% t(faprotax1_on_1_sel)
  
  #alpha
  ##otu
  alpha1OTU<-mean(diversity(x))#shannon
  ##funct
  alpha2func<-mean(diversity(funct_table_dia))
  ##phylo
  alpha3phy<-mean(mpd(x,tree_dist),na.rm =T)#mean pairwise distance, not ab weighted
  #beta
  ##otu
  beta1OTU<-mean(vegdist(x[rowSums(x)>0,], method="bray", binary=T))
  beta2OTU<-Cscr(x)#c-scores
  beta3OTU<-mean(rnkvar(x),na.rm =T)#rank_var
  ##func
  beta4func<-mean(vegdist(funct_table_dia[rowSums(funct_table_dia)>0,], method="bray", binary=T))
  ##phylo
  beta5phy<-mean(comdist(x,tree_dist),na.rm =T)
  #beta_with_ab
  ##otu
  brel1OTU<-mean(vegdist(x, method="bray", binary=F))
  ##func
  brel2func<-mean(vegdist(funct_table_dia[rowSums(funct_table_dia)>0,], method="bray", binary=F))
  ##phylo
  brel3phy<-mean(comdist(x,tree_dist,abundance.weighted =T))
  
  
  return(c('alpha1OTU'=alpha1OTU, 'alpha2func'=alpha2func, 'alpha3phy'=alpha3phy,
           'beta1OTU'=beta1OTU,'Cscore'=beta2OTU,'rank_var'=beta3OTU,
           'beta4func'=beta4func,'beta5phy'=beta5phy,
           'brel1OTU'=brel1OTU,'brel2func'=brel2func,'brel3phy'=brel3phy,'type'=type))
}


#data ####
load("randomisation_example_data.RData")
#isolate_table: dataframe of rarified bacterial communities in the diatom isolate cultures (samples as rows)
#sediment_table: dataframe of rarified bacterial communities in the source samples (samples as rows)
#faprotax1_on_1: dataframe of the functional annotations derived from FAPROTAX, created by using an identity matrix as input
#tree_16S:16S based phylogenetic tree
##note: to stratify within communities coming from the same source, the first part (before the underscore) of the sample name in isolate identifier is identical to the second part of the sample name  in the sediment identifier 

#functional annotations table
faprotax1_on_1_sel<-dplyr::select(faprotax1_on_1,dplyr::one_of(colnames(isolate_table)))
funct_table_isolate<-isolate_table %*% t(faprotax1_on_1_sel)

#cophenetic phylogenetic tree
tree_dist<-cophenetic(tree_16S)


#1 random ####
#permaswap algoritm with strata (within location)

grps<-as.factor(sapply(strsplit(row.names(isolate_table),'_'),unlist)[1,])
grps_sed<-as.factor(sapply(strsplit(row.names(sediment_table),'_'),unlist)[2,])
nr_times<-5 #number of randomisations

OTUswap<-permatswap(isolate_table[which(levels(grps)[1]==grps),],method="quasiswap",mtype="count",times=nr_times)
for (i in levels(grps)){print(i)
  OTUswap_b<-permatswap(isolate_table[which(i==grps),],method="quasiswap",mtype="count",times=nr_times)
  for (j in 1:nr_times){
    OTUswap$perm[[j]]<-rbind(OTUswap$perm[[j]],OTUswap_b$perm[[j]])
  }}

#2 OTU lottery from source community####

selectsed<-sediment_table[,which(colnames(sediment_table)%in%names(which(colSums(isolate_table)>0)))]
OTUlottery<-data.frame(isolate_table)
OTUlottery$iteration<-0
nr_times<-5 #number of randomisations
for (a in 1:nr_times){
reshuf<-isolate_table[1:2,]
for (i in 1:nrow(isolate_table)){
  row<-isolate_table[i,]
  row<-row[row>0]
  sedrow<-selectsed[which(grps_sed==grps[i]),]
  names(row)<-sample(names(sedrow), length(row), replace = FALSE, prob = sedrow)
  row<-t(data.frame(row)); 
  reshuf<-plyr::rbind.fill.matrix(reshuf,row)
}
reshuf<-reshuf[-c(1:2),]
row.names(reshuf)<-row.names(isolate_table)
reshuf<-data.frame(reshuf)
reshuf$iteration<-a
OTUlottery<-rbind(OTUlottery,reshuf)
}

OTUlottery[is.na(OTUlottery)]<-0

#calculate indexes####
#observed communities
#alpha
##otu
alpha1OTU<-mean(diversity(isolate_table))#shannon
##funct
alpha2func<-mean(diversity(funct_table_isolate))
##phylo
alpha3phy<-mean(mpd(isolate_table,tree_dist),na.rm =T)#mean pairwise distance, not ab weighted
isolate_table2<-isolate_table
isolate_table2<-isolate_table2[,colSums(isolate_table2)>0]
#beta
##otu
beta1OTU<-mean(vegdist(isolate_table, method="bray", binary=T))
beta2OTU<-Cscr(isolate_table)#c-scores
beta3OTU<-mean(rnkvar(isolate_table),na.rm =T)#rank_var
##func
beta4func<-mean(vegdist(funct_table_isolate[rowSums(funct_table_isolate)>0,], method="bray", binary=T))
##phylo
beta5phy<-mean(comdist(isolate_table,tree_dist))

#beta_with_ab
##otu
brel1OTU<-mean(vegdist(isolate_table, method="bray", binary=F))
##func
brel2func<-mean(vegdist(funct_table_isolate[rowSums(funct_table_isolate)>0,], method="bray", binary=F))
##phylo
brel3phy<-mean(comdist(isolate_table,tree_dist,abundance.weighted =T))


index_results<-data.frame('alpha1OTU'=alpha1OTU, 'alpha2func'=alpha2func, 'alpha3phy'=alpha3phy,
                          'beta1OTU'=beta1OTU,'Cscore'=beta2OTU,'rank_var'=beta3OTU,
                          'beta4func'=beta4func,'beta5phy'=beta5phy,
                          'brel1OTU'=brel1OTU,'brel2func'=brel2func,'brel3phy'=brel3phy, 'type'='original')
index_results$type<-as.character(index_results$type)

#1 random
for (counter in 1:attr(OTUswap, "times")){ 
  index_results<-rbind(index_results,indexes(as.matrix(as.data.frame(OTUswap$perm[counter])),type='random'))
}
index_results01<-index_results

#2 OTU lottery
for (a in 1:max(OTUlottery$iteration)) {
  index_results<-rbind(index_results,indexes(subset(OTUlottery[,-ncol(OTUlottery)],OTUlottery$iteration==a),type='OTUlottery'))
}
index_results02<-index_results


#plot indexes####
index_results_joined<-rbind(data.frame(index_results01),subset(data.frame(index_results02),type!='original')) #join multiple simulations

index_results_joined[1:11]<-sapply(index_results_joined[1:11],function (x) (as.numeric(as.character(x))))
index_res_long<-reshape2::melt(index_results_joined,idvar='type');index_res_long$groups<-NA
index_res_long$groups[grep("alpha", index_res_long$variable)]<-'alpha';index_res_long$groups[grep("beta", index_res_long$variable)]<-'beta'
index_res_long$groups[grep("brel", index_res_long$variable)]<-'beta_ab'
index_res_long2<-subset(index_res_long,groups=='alpha'|groups=='beta'|groups=='beta_ab')

ggplot() + geom_histogram(data=subset(index_res_long2,type=='OTUlottery'),aes(value), bins=50, fill='green')+
  geom_histogram(data=subset(index_res_long2,type=='random'),aes(value), bins=50, fill='blue')+
  geom_vline(data=subset(index_res_long2,type=='original'),aes(xintercept=value), color='darkred', size=1.2)+
  facet_wrap(variable~groups, scales='free')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))

