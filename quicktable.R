#testchange
genes<-c("Solyc07g008670","Solyc09g005750","Solyc06g083660","Solyc09g063080","Solyc07g032710","Solyc03g006840")
corrinter<-corr975tab[grep(paste(genes,collapse = "|"),corr975tab$GOI),c(1,3)]
SVinter<-SV_intersect3[grep(paste(corrinter[,2],collapse = "|"),SV_intersect3[,1]),]

final<-c()
for (r in 1:nrow(SVinter)) {
  gois <- paste(corrinter[grep(SVinter[r,1],corrinter[,2]),1], collapse = ",")
  final <- rbind(cbind(gois,SVinter[r,]),final)
}

#takes genesvector and subset correlation table of desired threshold by GOIs
#seperately rbind SVintersectCDS wtih SVintersectRE with the desired RE_column (for both)
   # OVATE or SUN subset in RE column, and add column for intersect CDS
# subset SVintersectcombined with only those that have affected genes that are matched in the correlated genes of the subsetted 
                                                                                                                   # corr table
# go through SVinter and column bind gois from genegroup (the only gois remaining in subsetted corrtable) whose correlated gene
      # match that of the affected gene in SVinter (it is possible to have multiple gois who are correlated to this gene so need
      # to paste/collapse gois into that cell)
#Also read in Rel_Expr_file or Meristems/Anthesis and find what affected gene protein is a cbind 
# make sure columns are in logical order (GOI(s),correlated/affected gene, Protein, SVstart, SVend, accessions, intersect CDS, RelExpr)  
# write table 
#see if you can just give it a folder or something and it will solve it (don't need to specify every file)

condensedtab <- function(genesv,corrtable,SVs_intersect_CDS_file,SVs_intersect_RE_file,Anthesis) {
  corrtable_file<-"Correlation_Table"
  SVs_intersect_CDS_file<-"1e6maxSV/SVs_intersect_CDS.txt"
  SVs_intersect_RE_file<-"1e6maxSV/5000fiveprimewindow/SVs_intersect_RE.txt"
  humtab_file<-"../../humtab"
  Rel_Expr_file<-"../../../Desktop/SiftStrucVars3.0/Anthesis.txt"
  humtab<-read.delim(humtab_file)
  genesv<-as.character(humtab$genesv)
  corrtable<-read.delim(corrtable_file)
  corrtable<-corrtable[grep(paste(genesv,collapse = "|"),corrtable$GOI),c("GOI","correlated_gene")]
  SV_intersect_CDS<-read.delim(SVs_intersect_CDS_file)
  SV_intersect_CDS<-SV_intersect_CDS[grep("OVATE",SV_intersect_CDS[,13]),c(1:5,13)]
  SV_intersect_CDS<-cbind(SV_intersect_CDS, intersect_CDS = rep("yes",nrow(SV_intersect_CDS)))
  if (SVs_intersect_RE_file == 0) {
    SV_inter<-SV_intersect_CDS
  } else {
    SV_intersect_RE<-read.delim(SVs_intersect_RE_file)
    SV_intersect_RE<-SV_intersect_RE[grep("OVATE",SV_intersect_RE[,7]),]
    intersect_CDS<-gsub("Yes","No in 5'UTR",SV_intersect_RE$intersect_5_prime)
    SV_intersect_RE<-cbind(SV_intersect_RE[,c(1:5,7)], intersect_CDS)
    SV_inter<-rbind(SV_intersect_CDS,SV_intersect_RE)
  } 
  grep1<-grep(paste(corrtable$correlated_gene[1:2538],collapse = "|"),SV_inter$affected_gene)
  grep2<- grep(paste(corrtable$correlated_gene[2539:4500],collapse = "|"),SV_inter$affected_gene) 
  grep3<- grep(paste(corrtable$correlated_gene[4501:6500],collapse = "|"),SV_inter$affected_gene) 
  grep4<- grep(paste(corrtable$correlated_gene[6501:8500],collapse = "|"),SV_inter$affected_gene) 
  SV_inter<-SV_inter[unique(c(grep1,grep2,grep3,grep4)),]
  Anthesis<-read.delim(Rel_Expr_file)
  final<-c()
  for (r in 1:nrow(SV_inter)) {
    gois <- corrtable[grep(SV_inter[r,4],corrtable$correlated_gene),1]
    goiHumNames<-c()
    for (goi in gois) {
      goiHumName<-paste(goi,"(",humtab[grep(goi,humtab$genesv),2],")",sep = "")
      goiHumNames<-append(goiHumNames,goiHumName)
    }
    protein<- as.character(Anthesis[grep(SV_inter[r,4],Anthesis$gene.name),c("description")])
    final <- rbind(cbind(gois = paste(goiHumNames,collapse = ","),SV_inter[r,],protein),final)
  }
  final<-cbind(final[,1:7],"sv_size_(bp)" = final$V3-final$V2,final[,8:9])
  return(final)
}

write.table(final,"../../../Downloads/corr99_50000wind_1e6maxSV_OVATEonly.txt",row.names = F,sep = "\t")

final<-final[!(grepl("DEL",final$svtype) & grepl("\\+",final$intersect_CDS)),]
write.table(final,"../corr0975withplusdels",row.names = F,sep = "\t")













#####

x<-group3599
test<-c()
for (r in 1:nrow(x)) {
  test<-rbind(paste(x[r,]$V2,x[r,]$V3,collapse = ","),test)
}

#####
group<-group3
groupnew<-c()
for (r in 1:nrow(group)) {
  if (group[r,]$intersect_CDS == "no") {
    if (any(grep(group[r,]$affected_gene,RE5000$affected_gene))) {
      groupnew<-rbind(group[r,],groupnew)
    } 
  } else { groupnew <- rbind(group[r,],groupnew)}  
}

####
group<-group35
groupnew<-c()
for (r in 1:nrow(group)) {
  if (any(grep(group[r,]$affected_gene,corr99$correlated_gene))){
    groupnew<-rbind(group[r,],groupnew)
  }
}

##### Making automated condensed table script 
# use it when put all output folders (corrs) in same directory (theoretically can be used in SiftStrucVars3.0 itself)
# do not have genes repeated (if greater than 0.99 don't have it in 0.95 or 0.975 etc.)
# just setwd() to appropriate directory 
corrs<-sort(list.dirs()[grep("corr$",list.dirs())],decreasing = T) 
library(stringr)
corrsnumbs<-str_extract(corrs,"(?<=\\./).+(?=c)")
detach(package:stringr)
final<-c()
for (i in 1:length(corrs)) {
  setwd(corrs[i])  
  corrtable_file<-"Correlation_Table"
  SVs_intersect_CDS_file<-"1e6maxSV/SVs_intersect_CDS.txt"
  SVs_intersect_RE_file<-"1e6maxSV/5000fiveprimewindow/SVs_intersect_RE.txt"
  humtab_file<-"../../humtab"
  Rel_Expr_file<-"../../../Desktop/SiftStrucVars3.0/Anthesis.txt"
  humtab<-read.delim(humtab_file)
  genesv<-as.character(humtab$genesv)
  corrtable<-read.delim(corrtable_file)
  if (i > 1) { 
     corrtable<-corrtable[corrtable$correlation <= corrsnumbs[(i-1)],]
  }
  corrtable<-corrtable[grep(paste(genesv,collapse = "|"),corrtable$GOI),]
  SV_intersect_CDS<-read.delim(SVs_intersect_CDS_file)
  SV_intersect_CDS<-SV_intersect_CDS[grep("OVATE",SV_intersect_CDS[,13]),c(1:5,13)]
  SV_intersect_CDS<-cbind(SV_intersect_CDS, intersect_CDS = rep("yes",nrow(SV_intersect_CDS)))
  # fix like try and except in python 
  if (SVs_intersect_RE_file == 0) {
    SV_inter<-SV_intersect_CDS
  } else {
    SV_intersect_RE<-read.delim(SVs_intersect_RE_file)
    SV_intersect_RE<-SV_intersect_RE[grep("OVATE",SV_intersect_RE[,7]),]
    intersect_CDS<-gsub("Yes","No in 5'UTR",SV_intersect_RE$intersect_5_prime)
    SV_intersect_RE<-cbind(SV_intersect_RE[,c(1:5,7)], intersect_CDS)
    SV_inter<-rbind(SV_intersect_CDS,SV_intersect_RE)
  } 
  grep1<-grep(paste(corrtable$correlated_gene[1:2538],collapse = "|"),SV_inter$affected_gene)
  grep2<- grep(paste(corrtable$correlated_gene[2539:4500],collapse = "|"),SV_inter$affected_gene) 
  grep3<- grep(paste(corrtable$correlated_gene[4501:6500],collapse = "|"),SV_inter$affected_gene) 
  grep4<- grep(paste(corrtable$correlated_gene[6501:8500],collapse = "|"),SV_inter$affected_gene) 
  SV_inter<-SV_inter[unique(c(grep1,grep2,grep3,grep4)),]
  Anthesis<-read.delim(Rel_Expr_file)
  onecorr<-c()
  for (r in 1:nrow(SV_inter)) {
    gois <- corrtable[grep(SV_inter[r,4],corrtable$correlated_gene),c("GOI","correlation")]
    goiHumNames<-c()
    for (it in 1:nrow(gois)) {
      goiHumName<-paste(gois[it,c("GOI")],"(",humtab[grep(gois[it,c("GOI")],humtab$genesv),2],")","[",gois[it,c("correlation")],"]",sep = "")
      goiHumNames<-append(goiHumNames,goiHumName)
    }
    protein<- as.character(Anthesis[grep(SV_inter[r,4],Anthesis$gene.name),c("description")])
    onecorr <- rbind(cbind(gois = paste(goiHumNames,collapse = ","),SV_inter[r,],protein),onecorr)
  }
  onecorr <-cbind("corrthresh" = rep(corrsnumbs[i]),onecorr[,1:7],"sv_size_(bp)" = onecorr$V3-onecorr$V2,onecorr[,8:9])
  final<-rbind(onecorr,final)
  setwd("../")
}