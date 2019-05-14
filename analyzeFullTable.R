fulltab<-read.delim("Documents/LA1589_SA1_2_29_32_fulltable.txt")
uniqgenes<-as.character(unique(fulltab$affected_gene))
uniqgenesrows<-sapply(uniqgenes, function(x) grep(x,fulltab$affected_gene))
#test<-data.frame(uniqgenes,uniqgenesrows)
uniqgenesexprdataorigin<-sapply(uniqgenes, function(x) unique(as.character(fulltab[grep(x,fulltab$affected_gene),c("Expr_Data_Origin")])))
uniqgenesexprdataoriginnumbs<-sapply(uniqgenes, function(x) length(unique(as.character(fulltab[grep(x,fulltab$affected_gene),c("Expr_Data_Origin")]))))

#just subset those that have correlation in all 5 accessions 
geneswith5ormoreexpr<-uniqgenes[uniqgenesexprdataoriginnumbs>=5]
genesinallTab<-fulltab[grep(paste(geneswith5ormoreexpr,collapse = "|"),fulltab$affected_gene),]

#subset those those that have at least 3 corr that is 0.99>x>=0.975 or 1 corr that is 0.99 (keep in 2 seperate tabs)
  # take only those that have 1 that is 0.99; then if it doesn't see if it has at least 3 that are 0.99>x>=0.975
#make intermediate tab with only 0.975,find uniq affected genes inside, for each find number of exprdataorigins and only take those that are x>=3

Genesinall099<-genesinallTab[grep(paste(unique(genesinallTab[grep(0.99,genesinallTab$corrthresh),c("affected_gene")]),collapse = "|"),genesinallTab$affected_gene),]

inter0975<-genesinallTab[grep(0.975,genesinallTab$corrthresh),]
uniq0975<-unique(inter0975$affected_gene)
exprdataorignumbs0975<-sapply(uniq0975, function(g) length(unique(inter0975[grep(g,inter0975$affected_gene),"Expr_Data_Origin"])))
with3ormoreat0975<-uniq0975[exprdataorignumbs0975>=3]
Genesinall0975<-genesinallTab[grep(paste(with3ormoreat0975,collapse = "|"),genesinallTab$affected_gene),]
Genesinall0975<-Genesinall0975[!grepl(paste(Genesinall099$affected_gene,collapse = "|"),Genesinall0975$affected_gene),]


#SUM TOTAL: 1) affected gene has 0.95 corr in all 5 Expr_datasets, 2) has at least 1 at x>=0.99 or if NOT has at least 3 at 0.99>x>=0.975 



#WASTE OF TIME 
  ## use [[]] to acess when putting together as data frame for sapply/lapply outputs 
  #find how many gois/what gois correlate at 0.95,975,99 for each uniqgene ; make its own column for each ;
      #have to look at corrthresh value and then count commas in gois   