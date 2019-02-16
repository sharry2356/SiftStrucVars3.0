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
  corrtable<-read.delim(corrtable)
  corrtable<-corrtable[grep(paste(genesv,collapse = "|"),corrtable$GOI),c("GOI","correlated_gene")]
  SV_intersect_CDS<-read.delim(SVs_intersect_CDS_file)
  SV_intersect_CDS<-SV_intersect_CDS[grep("OVATE",SV_intersect_CDS[,12]),c(1:4,12)]
  SV_intersect_CDS<-cbind(SV_intersect_CDS, intersect_CDS = rep("yes",nrow(SV_intersect_CDS)))
  if (SVs_intersect_RE == 0) {
    SV_inter<-SV_intersect_CDS
  }
  else {
    SV_intersect_RE<-read.delim(SVs_intersect_RE_file)
    SV_intersect_RE<-SV_intersect_RE[grep("OVATE",SV_intersect_RE[,6]),c(1:4,6)]
    SV_intersect_RE<-cbind(SV_intersect_RE, intersect_CDS = rep("no",nrow(SV_intersect_RE)))
    SV_inter<-rbind(SV_intersect_CDS,SV_intersect_RE)
  } 
  grep1<-grep(paste(corrtable$correlated_gene[1:2538],collapse = "|"),SV_inter$affected_gene)
  grep2<- c(grep1,grep(paste(corrtable$correlated_gene[2539:3972],collapse = "|"),SV_inter$affected_gene)) 
  SV_inter<-SV_inter[unique(grep2),]
  Anthesis<-read.delim(Rel_Expr_file)
  final<-c()
  for (r in 1:nrow(SV_inter)) {
    gois <- paste(corrtable[grep(SV_inter[r,3],corrtable$correlated_gene),1], collapse = ",")
    protein<- as.character(Anthesis[grep(SV_inter[r,3],Anthesis$gene.name),c("description")])
    final <- rbind(cbind(gois,SV_inter[r,],protein),final)
  }
  return(final)
}

write.table(final,"../../../Downloads/corr99_50000wind_1e6maxSV_OVATEonly.txt",row.names = F,sep = "\t")















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

#####
group<-group3599
groupnew<-c()
for (r in 1:nrow(group)) {
  protein <- relexpr[grep(group[r,]$affected_gene,relexpr$gene),]$Protein
  groupnew <- rbind(cbind(group[r,],protein),groupnew)
}
