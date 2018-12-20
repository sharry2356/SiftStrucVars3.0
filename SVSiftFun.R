###############Make the following things parameters for the functions
############### Correlation threshold; Max # of genes (semicolons) in an SV
############### For defaults, your current values are fine

#remember to +1 or -1 to max SV_length based on 0 or 1 base as well 
# remove "NAMES=" from accessions?
Visa_SVtobed <- function(filename,max_SVlength=50000){
  library(stringr)
  
  svLines <- readLines(filename)
  svLines<-svLines[!grepl("^#",svLines)]
  
  ends <- as.numeric(str_extract(svLines, "(?<=END=)[-0123456789]+(?=;)"))
  chromNames <- str_extract(svLines, "^.+ch\\d{2}(?=\t)")
  starts <- as.numeric(str_extract(svLines, "(?<=ch\\d{2}\t)\\d+(?=\t)"))-1 
  svType <- str_extract(svLines, "(?<=SVTYPE=)\\w+(?=;)")
  accessions <- str_extract(svLines, "SNAME=.*;")
  
  AllSVs <- data.frame(chromNames, starts, ends, accessions, svType, stringsAsFactors = F)
  bedformat<-unique(AllSVs[grep("DEL",AllSVs$svType),1:4])
  bedformat<-bedformat[ bedformat$ends - bedformat$starts < max_SVlength,]
  bedformat$chromNames<-str_remove(bedformat$chromNames,"SL\\d+\\.\\d+")
  bedformat$starts<-as.integer(bedformat$starts)
  bedformat$ends<-as.integer(bedformat$ends)
  startends<-c()
  for (r in 1:nrow(bedformat)) {
    startend<-paste(bedformat$starts[r],bedformat$ends[r],sep = ",") 
    startends<-append(startends,startend) 
  }
  uniqstartends<-unique(startends)
  for (sei in uniqstartends) {
    indicesstartend<-grep(sei,startends)
    accessions<-bedformat$accessions[c(indicesstartend)]
    accessions<-paste(accessions,collapse = ",")
    for (index in indicesstartend) {
      bedformat$accessions[index]<-accessions
    }
  }
  bedformat<-unique(bedformat)
  write.table(bedformat, "bedformat_SVs.bed", col.names = F, quote=F, row.names=F, sep="\t")
  detach(package:stringr)
}


####this isn't universal yet...make the fixing of the formatting (more unifrom) of marksgenes its own function
###(the first part before the first for loop)
AllCorrelations_IMPROVED<- function(goiList_FILE,expression_table_FILE,expr_cols,corrthresh = 0.995){
  expression_table <- read.csv(expression_table_FILE)
  for (col in expr_cols) {
    expression_table[,col]<-as.numeric(expression_table[,col])
  }
  expression_table[,1]<-as.character(expression_table[,1])
  goilogic <- readLines(goiList_FILE) == ""
  goiList <- readLines(goiList_FILE)
  goiList <- goiList[!goilogic]
  library(stringr)
  Cs<-c()
  for (GOI in goiList) {
    GOIrow<-grep(GOI,expression_table[,1])
    if (length(GOIrow)==1) {
      for(r in 1:nrow(expression_table)){
        correlation<-cor(as.numeric(expression_table[GOIrow,expr_cols]),as.numeric(expression_table[r,expr_cols]))
        if (correlation > corrthresh) {
          Cs<-rbind(Cs,data.frame(GOI,correlation, correlated_gene = expression_table[r,1]))
        } else{next}
      }
    } else{next}
  }  
  Cs$correlated_gene<-str_extract(Cs$correlated_gene,"Solyc\\d*g\\d*")
  write.table(Cs,"Correlation_Table",sep = "\t")
  CsandGOIs <- unique(c(as.character(Cs[,3]), goiList))
  writeLines(CsandGOIs, "GOIs_and_Correlatedgenes_list.txt")
  detach(package:stringr)
}

###19636902 is average size of 100 genes for each SV.   min is 657876         Make max_SVlength???  9300396 (from SV_deletion_analysis)
####myabe need to subtract 1 from starts if its not in bed (not 0 based)
Lemon_SVtobed<- function(filename,max_SVlength=50000) {
  svTabtest <- read.delim(filename, 
                          header=T, sep="\t",
                          colClasses = c("character", "character", "character", 
                                         "character", "character", "character",
                                         "character","character", "logical", 
                                         "character", "logical","logical",
                                         "logical","logical","logical"))
  svTabtest <- svTabtest[grepl("DEL", svTabtest$SVtype),]
  library(stringr)
  chromNames<-c()
  starts<-c()
  ends<-c()
  accessions<-c()
  jump<- str_split(svTabtest$Positions,"-")
  for(r in 1:length(jump)){
    inter<- str_split(jump[[r]],":")
    chromNames<- append(inter[[1]][1],chromNames)
    starts<- append(inter[[1]][2],starts)
    ends<- append(inter[[2]][2],ends)
    accessions<- append(svTabtest$accessions[r],accessions)
  }  
  bedformat<-data.frame(chromNames,starts = as.numeric(starts)-1,ends = as.numeric(ends),accessions)
  bedformat<-bedformat[ bedformat$ends - bedformat$starts < max_SVlength,]
  bedformat$chromNames<-str_remove(bedformat$chromNames,"SL\\d+\\.\\d+")
  bedformat$starts<-as.integer(bedformat$starts)
  bedformat$ends<-as.integer(bedformat$ends)
  write.table(bedformat, "bedformat_SVs.bed", col.names = F, quote=F, row.names=F, sep="\t")
  detach(package:stringr)
}

#which SV
SVs_to_Bedformat<- function(senderofSV, filename, max_SVlength=50000){
  if (senderofSV=="Lemon") {
    return(Lemon_SVtobed(filename,max_SVlength))
  } else if (senderofSV=="Visa") {
    return(Visa_SVtobed(filename,max_SVlength))
  } else {print("Names so far are either 'Visa' or 'Lemon' ")}
}

#change annotated genome to 0 based (NOT NEEEDED, bedtoos recognizes gff as 1 based)
annotatedgenome_togff <- function(annotatedgenome_FILE) {
  library(stringr)
  annotatedgenomeinter<-readLines(annotatedgenome_FILE)
  annotatedgenomeinter<-annotatedgenomeinter[!grepl("^#", annotatedgenomeinter)]
  writeLines(annotatedgenomeinter, "annotatedgenomeinter.txt",)
  annotatedgenomeinter<-read.delim("annotatedgenomeinter.txt",header = F)
  annotatedgenomeinter$V4<-as.integer(annotatedgenomeinter$V4-0)
  annotatedgenomeinter$V5<-as.integer(annotatedgenomeinter$V5)
  annotatedgenomeinter$V1<-str_remove(annotatedgenomeinter$V1,"SL\\d+\\.\\d+")
  write.table(annotatedgenomeinter,"annotatedgenome.gff", col.names = F, quote=F, row.names=F, sep="\t")
  detach(package:stringr)
}

#strip down to integers for later mutationmath func
# input: outputs from bedtools    output: table to be analyzed by function mutationmath 
Strip_to_Integers<- function(Full_CDS_sets_FILE, overlapping_SVs_and_GOIsandCs_CDS_FILE){
  library(stringr)
  Full_CDS_sets <- read.delim(Full_CDS_sets_FILE, header = F)
  overlapping_SVs_and_GOIsandCs_CDS <- read.delim(overlapping_SVs_and_GOIsandCs_CDS_FILE, header = F)
  overlapping_SVs_and_GOIsandCs_CDS$V2<-overlapping_SVs_and_GOIsandCs_CDS$V2 + 1 
  Full_CDS_sets$V4<-Full_CDS_sets$V4 + 0
  overlapping_SVs_and_GOIsandCs_CDS<- cbind(Rel_Id = rep(1:nrow(overlapping_SVs_and_GOIsandCs_CDS)), overlapping_SVs_and_GOIsandCs_CDS)
  Full_CDS_sets[,9] <- str_extract(Full_CDS_sets[,9],"Solyc\\d*g\\d*")
  genes_with_CDS_overlapped_by_SVs<- unique(Full_CDS_sets[,9])
  detach(package:stringr)
  Stripped_for_Analysis<- c()
  for (r in 1:length(genes_with_CDS_overlapped_by_SVs)) {
    namelogic_CDS<- grepl(genes_with_CDS_overlapped_by_SVs[r], Full_CDS_sets[,9])
    genesubset_CDS<- Full_CDS_sets[namelogic_CDS,1:9]
    genesubset_CDS_integers_combined<-c()
    for (z in 1:nrow(genesubset_CDS)){
      start_to_end<- genesubset_CDS[z,4]:genesubset_CDS[z,5]
      genesubset_CDS_integers_combined <- append(start_to_end, genesubset_CDS_integers_combined)
    }
    namelogic_SVs<- grepl(genes_with_CDS_overlapped_by_SVs[r], overlapping_SVs_and_GOIsandCs_CDS[,14])
    genesubset_SVs<- unique(overlapping_SVs_and_GOIsandCs_CDS[namelogic_SVs,3:5])
    genesubset_SVs<-cbind(genesubset_SVs[,1:2],"integers_inCDS_notin_SVs" = vector(, nrow(genesubset_SVs)),
                          "affected_gene"= rep(genes_with_CDS_overlapped_by_SVs[r], nrow(genesubset_SVs)), 
                          "combinedCDSlength" = rep(length(genesubset_CDS_integers_combined),nrow(genesubset_SVs)),
                          "strand" = rep(genesubset_CDS[1,7],nrow(genesubset_SVs)), 
                          "Low CDSbp" = rep(min(genesubset_CDS_integers_combined),nrow(genesubset_SVs)),
                          "High CDSbp" = rep(max(genesubset_CDS_integers_combined),nrow(genesubset_SVs)),
                          "accessions" = genesubset_SVs[,3])
    for(y in 1:nrow(genesubset_SVs)){
      start_to_end<- genesubset_SVs[y,1]:genesubset_SVs[y,2]
      intersect_CDS<-intersect(start_to_end,genesubset_CDS_integers_combined)
      notinSV_but_inCDS <- list(setdiff(genesubset_CDS_integers_combined,intersect_CDS))
      genesubset_SVs[y,3]<-list(notinSV_but_inCDS)
    }
    Stripped_for_Analysis<- rbind(genesubset_SVs, Stripped_for_Analysis)
  }
  return(Stripped_for_Analysis)
}
#Analyze output from Strip_to_Integers to determine the mutation of each SV on corresponding gene 
# input: Stripped_for_Analysis output: SV and gene mutation table (as described in above line)  
MutationMath<- function(Stripped_for_Analysis){
  SVs_intersect_CDS<-c()
  for(r in 1:nrow(Stripped_for_Analysis)){
    del_gene_length<-length(Stripped_for_Analysis[r,3][[1]])  
    if (del_gene_length == 0) {
      whole_gene_deleted<-"yes"
      in_frame<- "no"
      framshift<- "no"
      percentofgene_deleted<- 100
      start_deleted<- "yes"
      end_deleted<- "yes"
      likely_knockout<- "yes"
    }
    else{
      whole_gene_deleted<-"no"
      percentofgene_deleted<- (Stripped_for_Analysis[r,5]-del_gene_length)/(Stripped_for_Analysis[r,5])*100
      if ("+" %in% Stripped_for_Analysis[r,6]) {
        if(Stripped_for_Analysis[r,7] %in% Stripped_for_Analysis[r,3][[1]]){
          start_deleted<-"no"
        }
        else{
          start_deleted<-"yes"
        }
        if(Stripped_for_Analysis[r,8] %in% Stripped_for_Analysis[r,3][[1]]){
          end_deleted<-"no"
        }
        else{
          end_deleted<-"yes"
        }
      }
      else{
        if(Stripped_for_Analysis[r,8] %in% Stripped_for_Analysis[r,3][[1]]){
          start_deleted<-"no"
        }
        else{
          start_deleted<-"yes"
        }
        if(Stripped_for_Analysis[r,7] %in% Stripped_for_Analysis[r,3][[1]]){
          end_deleted<-"no"
        }
        else{
          end_deleted<-"yes"
        }
      }
      if(start_deleted=="yes"){
        likely_knockout<-"likely (check sequence)"
        framshift<-"likely (check sequence)"
        in_frame<-"maybe (check sequence)"
      }
      else if (end_deleted=="yes") {
        likely_knockout<-"not likely (check sequence)"
        framshift<-"not relevent"
        in_frame<-"not relevent"
      }
      else{
        mult_of_3<- del_gene_length/3
        if (mult_of_3%%1==0) {
          in_frame<-"yes"
          framshift<-"no"
          likely_knockout<- "no"
        }
        else{
          in_frame<-"no"
          framshift<-"yes"
          likely_knockout<-"yes"
        }
      }
    }
    del_analyzed<- cbind(Stripped_for_Analysis[r,c(1,2,4,9)],whole_gene_deleted,  percentofgene_deleted,in_frame, framshift, start_deleted, end_deleted, likely_knockout)
    SVs_intersect_CDS<-rbind(del_analyzed,SVs_intersect_CDS)
  }
  write.table(SVs_intersect_CDS, "SVs_intersect_CDS.txt", sep = "\t")
  writeLines(as.character(unique(SVs_intersect_CDS[,3])), "SVs_intersect_CDS_genelist.txt",sep = "\n")
}

# combine SVs that intersect 5' and promoter region (both have corresponding genes affected too)
SVs_intersect_RE_FUNC<- function(SVs_2000_from_promoter_and_genes_FILE, SVs_only_intersect_fives_FILE) {
  SVs_2000_from_promoter_and_genes<-read.delim(SVs_2000_from_promoter_and_genes_FILE, header = F)
  SVs_only_intersect_fives<-read.delim(SVs_only_intersect_fives_FILE, header = F)
  SVs_2000_from_promoter_and_genes<- cbind(SVs_2000_from_promoter_and_genes, "intersect_5_prime" = rep("NO (within 
                                            specified bp upstream of 5 prime)", nrow(SVs_2000_from_promoter_and_genes)))
  SVs_only_intersect_fives<- cbind(SVs_only_intersect_fives, "intersect_5_prime" = rep("Yes", 
                                                                                       nrow(SVs_only_intersect_fives)))
  inter<-rbind(SVs_2000_from_promoter_and_genes,SVs_only_intersect_fives)
  library(stringr)
  SVs_intersect_RE<-data.frame(inter[,2:3], "affected_gene" = str_extract(inter[,13], "Solyc\\d*g\\d*"), 
                               "accessions" = inter[,4], "intersect_5_prime" = inter[,14])
  SVs_intersect_RE$V2<- SVs_intersect_RE$V2 + 1
  write.table(unique(SVs_intersect_RE),"SVs_intersect_RE.txt", sep = "\t")
  writeLines(as.character(unique(SVs_intersect_RE[,3])), "SVs_intersect_RE_genelist.txt",sep = "\n")
}
###need to add 1 to SVs to make 1 based again
###Generate list of genes to look up and view expression of and place downloads into appropriate directories  
Generate_genelistsFUNC<- function(SVs_intersect_CDS_FILE, SVs_intersect_RE_FILE) {  
  SVs_intersect_CDS<-read.delim(SVs_intersect_CDS_FILE)
  SVs_intersect_RE<-read.delim(SVs_intersect_RE_FILE)
  writeLines(as.character(unique(SVs_intersect_CDS[,3])), "SVs_intersect_CDS_genelist.txt",sep = "\n")
  writeLines(as.character(unique(SVs_intersect_RE[,3])), "SVs_intersect_RE_genelist.txt",sep = "\n")
}

###Read expression files downloaded from solgenomics (make sure each type in seperate folder)
###NEW
ReadFiles<- function(Meristems_file,Anthesis_file){
  Meristems_Expression<-read.delim(Meristems_file)
  Anthesis_Expression<-read.delim(Anthesis_file)
  Expression<-cbind(Meristems_Expression[,1:92], Anthesis_Expression[,2:179])
  return(Expression)
}


####OLD
ReadFilesONLINEDOWNLOAD_OLD_OLD_OLD<- function(Meristem_path_to_expression, Anthesis_path_to_expression) {
  original_dir<-getwd()
  setwd(Meristem_path_to_expression)
  geneFiles<-list.files(path = Meristem_path_to_expression, pattern = "custom_list")
  gene_expression<-c()
  for (geneexpressionfile in geneFiles) {
    readin<-read.delim(geneexpressionfile)
    gene_expression<-rbind(readin, gene_expression)
  }
  Meristems_expression<-gene_expression
  
  setwd(Anthesis_path_to_expression)
  geneFiles<-list.files(path = Anthesis_path_to_expression, pattern = "custom_list")
  gene_expression<-c()
  for (geneexpressionfile in geneFiles) {
    readin<-read.delim(geneexpressionfile)
    gene_expression<-rbind(readin, gene_expression)
  }
  Anthesis_expression<-gene_expression
  Meristems_expression<-Meristems_expression[order(as.character(Meristems_expression[,1])),]
  Anthesis_expression<-Anthesis_expression[order(as.character(Anthesis_expression[,1])),]
  Expression<-cbind(Meristems_expression[,1:92], Anthesis_expression[,2:179])
  setwd(original_dir)
  return(Expression)
}

####Analyze expression data to make it more comparable 
Expression_Analysis_Func<- function(Expression, zero.rm =T) {
  if (zero.rm == T) {
    Expression<-as.matrix(Expression)
    Expression<-gsub("^0$|\\s\\s0.00$",NA,Expression)
    Expression<-data.frame(Expression)
  }
  for (column in 2:(ncol(Expression)-1)) {
    Expression[,column]<-as.numeric(as.character(Expression[,column]))
  }
  Expression$gene.name<-as.character(Expression$gene.name)
  Expression$description<-as.character(Expression$description)
  TMcols<-grep("\\.TM",names(Expression))
  SIMcols<-grep("\\.SIM",names(Expression))
  FMcols<-grep("\\.FM",names(Expression))
  otherMeristemcols<-grep("\\.EVM|\\.MVM|\\.LVM|\\.SYM",names(Expression))
  Anthesis_psc_cols<-c(173, 189, 237)
  fiveDPA_psc_cols<-c(174,190,238)
  tenDPA_psc_cols<-c(175,191,239)
  Anthesis_non_cols<-setdiff(grep("\\.Anthesis",names(Expression)),c(173, 189, 237))
  five_dpa_non_cols<-setdiff(grep("\\.5DPA",names(Expression)),c(174,190,238))
  ten_dpa_non_cols<-setdiff(grep("\\.10DPA",names(Expression)),c(175,191,239))
  psc_non_A_5_10_cols<- setdiff(grep("Pericarp|Septum|Columella",names(Expression)), c(Anthesis_psc_cols,
                                                                                       fiveDPA_psc_cols, tenDPA_psc_cols))
  Expression_Analyzed<-c()
  for (r in 1:nrow(Expression)){
    TM<-mean(as.numeric(Expression[r,TMcols]),na.rm = zero.rm)
    SIM<-mean(as.numeric(Expression[r,SIMcols]),na.rm = zero.rm)
    FM<-mean(as.numeric(Expression[r,FMcols]),na.rm = zero.rm)
    TM_SIM_FM<-mean(c(TM,SIM,FM),na.rm = zero.rm)
    other<-mean(as.numeric(Expression[r,otherMeristemcols]),na.rm = zero.rm)
    TMSIMFM_over_other<- TM_SIM_FM/other
    Anthesis_psc<-mean(as.numeric(Expression[r,Anthesis_psc_cols]),na.rm = zero.rm)
    fiveDPA_psc<-mean(as.numeric(Expression[r,fiveDPA_psc_cols]),na.rm = zero.rm)
    tenDPA_psc<-mean(as.numeric(Expression[r,tenDPA_psc_cols]),na.rm = zero.rm)
    fiveDPA_psc_over_Anthesis_psc_Ratio<-fiveDPA_psc/Anthesis_psc
    fiveDPA_psc_over_tenDPA_psc_Ratio<-fiveDPA_psc/tenDPA_psc
    A_5_10_psc<-mean(c(Anthesis_psc,fiveDPA_psc,tenDPA_psc),na.rm = zero.rm)
    A_5_10_non_psc<-mean(as.numeric(Expression[r,c(Anthesis_non_cols,five_dpa_non_cols,ten_dpa_non_cols)]),na.rm = zero.rm)
    A510psc_over_non_psc<- A_5_10_psc/A_5_10_non_psc
    psc_non_A510<- mean(as.numeric(Expression[r,psc_non_A_5_10_cols]),na.rm = zero.rm)
    A510psc_over_non_A510<- A_5_10_psc/psc_non_A510
    row<- data.frame("gene" = Expression[r,1],"Protein" = Expression[r,270],"TM_average" = TM, 
                     "SIM_average" = SIM, "FM_average" = FM, "TM_SIM_FM_average" = TM_SIM_FM, 
                     "other_Meristem_average" = other, "TM_SIM_FM/other_Ratio" = TMSIMFM_over_other, 
                     Anthesis_psc, fiveDPA_psc, tenDPA_psc,fiveDPA_psc_over_Anthesis_psc_Ratio,
                     fiveDPA_psc_over_tenDPA_psc_Ratio,"Anthesis_5DPA_10DPA_psc_average" = A_5_10_psc, 
                     "Anthesis_5DPA_10DPA_NON_psc_average" = A_5_10_non_psc, 
                     "Anthesis_5DPA_10DPA_psc/Anthesis_5DPA_10DPA_NON_psc_Ratio" = A510psc_over_non_psc, 
                     "NON_Anthesis_5DPA_10DPA_psc_average" = psc_non_A510, 
                     "Anthesis_5DPA_10DPA_psc/NON_Anthesis_5DPA_10DPA_psc_Ratio" = A510psc_over_non_A510)
    Expression_Analyzed<-rbind(row,Expression_Analyzed)
  }
  return(Expression_Analyzed)
}


### Loop around Relevent Expression 


Rel_Expr_Bash_Loop<- function(Expression_Analyzed) {
  bashargs<-lapply(commandArgs(trailingOnly = T), function(x) eval(parse(text = x)))
  names(bashargs)<-c("TM_SIM_FM","TM_SIM_FM_over_other_meristem_ratio",
                     "Anthesis_fiveDPA_tenDPA_psc","A510psc_over_A510_NONpsc","A510psc_over_NON_A510_psc",
                     "fiveDPA_psc_over_A_psc_and_tenDPA_psc")
  Param_comb<-expand.grid(bashargs)
  for (r in 1:nrow(Param_comb)) {
    Expression_Analyzed<-Relevent_Expression(Expression_Analyzed,Param_comb[r,]$TM_SIM_FM,Param_comb[r,]$TM_SIM_FM_over_other_meristem_ratio,
                        Param_comb[r,]$Anthesis_fiveDPA_tenDPA_psc,Param_comb[r,]$A510psc_over_A510_NONpsc,
                        Param_comb[r,]$A510psc_over_NON_A510_psc,Param_comb[r,]$fiveDPA_psc_over_A_psc_and_tenDPA_psc)       
  }
  ReleventExpressionTab<-Expression_Analyzed
  return(ReleventExpressionTab)
}


###
Relevent_Expression<- function(Expression_Analyzed,TM_SIM_FM,TM_SIM_FM_over_other_meristem_ratio,
                               Anthesis_fiveDPA_tenDPA_psc,A510psc_over_A510_NONpsc,A510psc_over_NON_A510_psc,
                               fiveDPA_psc_over_A_psc_and_tenDPA_psc) {
  ReleventExpressionTab<-data.frame(Expression_Analyzed, "Relevent_Expression" = vector(,length = nrow(Expression_Analyzed)))
  NAtoF<- function(logic) {
    if (is.na(logic)) {
      return(FALSE)
    } else {return(logic)}
  }
  for (r in 1:nrow(Expression_Analyzed)) {
    if (NAtoF(any(c(Expression_Analyzed$TM_average[r],Expression_Analyzed$SIM_average[r],
                    Expression_Analyzed$FM_average[r]) > TM_SIM_FM))) {
      if (NAtoF(Expression_Analyzed$TM_SIM_FM.other_Ratio[r] > TM_SIM_FM_over_other_meristem_ratio)) {
        ReleventExpressionTab$Relevent_Expression[r]<-"OVATE"  
      }    
    }  
    if (NAtoF(any(c(Expression_Analyzed$Anthesis_psc[r],Expression_Analyzed$fiveDPA_psc[r],
                    Expression_Analyzed$tenDPA_psc[r]) > Anthesis_fiveDPA_tenDPA_psc))) {
      if (NAtoF(Expression_Analyzed$Anthesis_5DPA_10DPA_psc.Anthesis_5DPA_10DPA_NON_psc_Ratio[r] > 
                A510psc_over_A510_NONpsc &
                Expression_Analyzed$Anthesis_5DPA_10DPA_psc.NON_Anthesis_5DPA_10DPA_psc_Ratio[r] > 
                A510psc_over_NON_A510_psc)){
        if (NAtoF(all(c(Expression_Analyzed$fiveDPA_psc_over_Anthesis_psc_Ratio[r], 
                Expression_Analyzed$fiveDPA_psc_over_tenDPA_psc_Ratio[r]) > fiveDPA_psc_over_A_psc_and_tenDPA_psc))) {
          ReleventExpressionTab$Relevent_Expression[r]<-paste(ReleventExpressionTab$Relevent_Expression[r],
                                                              "SUN", sep = ",") 
        }
      }
    }
  }
  ReleventExpressionTab$Relevent_Expression<-gsub("FALSE,SUN","SUN",ReleventExpressionTab$Relevent_Expression)
  identity<-paste("RelExpr",TM_SIM_FM, TM_SIM_FM_over_other_meristem_ratio, Anthesis_fiveDPA_tenDPA_psc, A510psc_over_A510_NONpsc, 
                  A510psc_over_NON_A510_psc, fiveDPA_psc_over_A_psc_and_tenDPA_psc, sep = "_")
  names(ReleventExpressionTab)[max(ncol(ReleventExpressionTab))]<-identity
  return(ReleventExpressionTab)
}    

#grep gene name once; do not loop over parameter column names b/c this will cause it to need to grep genes constantly (more genes than parameters)
AddReleventExpressiontoSVs<- function(Relevent_ExpressionFILE, SVs_intersect_CDSorRE_FILE) {
  Relevent_Expression<-read.delim(Relevent_ExpressionFILE,)
  SVs_intersect<- read.delim(SVs_intersect_CDSorRE_FILE)
  RelExpr_names<-names(Relevent_Expression)[grepl("RelExpr",names(Relevent_Expression))]
  SVs_intersect<-data.frame(SVs_intersect,(sapply(RelExpr_names, function(x) vector(,length = nrow(SVs_intersect)))))
  for (r in 1:nrow(SVs_intersect)) {
    Relrow<-grep(as.character(SVs_intersect$affected_gene[r]), Relevent_Expression$gene)
    SVs_intersect[r,RelExpr_names]<-unlist(Relevent_Expression[Relrow,as.character(RelExpr_names)])
  }
  return(SVs_intersect)
}
  
  
  
  
  
