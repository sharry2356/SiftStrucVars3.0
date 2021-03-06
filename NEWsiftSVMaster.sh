#!/bin/bash

# make sure outupts are correct 
# create inputs for files on this bash script 
# delete unnecesary files 
#negative correlations

#Inputs---FILE NAMES 
#Except variables "SV_type" and "Expression_cols" which refer to what the formating of the SV_Table and which columns have expression data of Expression_data, respectively; so far for SV_type, "Visa" is the name for a vcf file and "Lemon" is the name for a file like the one Zach Lippman sent (the postdoc who did it was named Zach Lemmon); MAKE SURE Expression_cols IS A STRING OF A R VECTOR AS SEEN WITH THE "RELEVENTEXPRESSIONVARS A FEW LINES DOWN; skipCorr="False" or "F" when not trying to do correlation     
annotatedgenome="ITAG2.4_gene_models.gff3"
SV_table="master_vcf_5_overlap_len_20_to_1m.csv"
SV_type="NewVisa" 
Expression_data="markGeneList_RPKM_OvateCoexpr.csv" 
Expression_cols="c(2,4,5,3,7,6)"
goiList="TRMs_wM8.txt"
skipCorr="T"
#Parameter values 
# I believe 5' windows can only be in numeric notation (does not support scientific notation with "e") NEED TO CHECK THIS FIRST 
correlations=(0.99)
max_SVlengths=(1e6)
fivePrimeWindows=(5000)
#ReleventExpressionsVars--- MAKE SURE IN FORMAT OF R AND IN QUOTES--- order listed as arguments into Rscript ARE IMPORTANT (if adding more variables maintain order)  
TM_SIM_FM="c(5)" 
TM_SIM_FM_over_other_meristem_ratio="c(0.1)"
Anthesis_fiveDPA_tenDPA_psc="c(10)"
A510psc_over_A510_NONpsc="c(0.75)"
A510psc_over_NON_A510_psc="c(0.75)"
fiveDPA_psc_over_A_psc_and_tenDPA_psc="c(0.5)"

###############################################################################################################################################################################
#Everything Below is code (Do not need to adjust anything) !!!!!!



#Add file names to R scripts
if [[ "$skipCorr" == "False" ]] || [[ "$skipCorr" == "F" ]] 
then
  echo "YOU ARE SKIPPING CORRELATION!!!!"
  echo "You will get 2 errors regarding 'Correlation Table'; This is normal, disregard them; Errors after these  are also normal if for some parameters no SV overlap was found"
  correlations=(1.0)
else
  sed -i -E "s/\".*\",.*,/\"$goiList\",expression_table_FILE=\"$Expression_data\",$Expression_cols,/" Rscripts/siftSVscorr.R 
fi
sed -i -E "s/\".*\",/\"$SV_type\",filename=\"$SV_table\",/" Rscripts/siftSVsmaxSV.R 
sed -i -E "s/annotatedgenome_FILE=\".*\"/annotatedgenome_FILE=\"$annotatedgenome\"/" Rscripts/siftSVsag_togf.R

#Start running pipeline 
Rscript Rscripts/siftSVsag_togf.R

for max_SVlength in ${max_SVlengths[*]}
do
  sed -i -E "s/max_SVlength\\=[0-9].*\\)/max_SVlength\\=$max_SVlength\\)/" Rscripts/siftSVsmaxSV.R
  Rscript Rscripts/siftSVsmaxSV.R   
  mkdir -p ${max_SVlength}maxSV
  mv -t ${max_SVlength}maxSV bedformat_SVs.bed
done 	

for correlation in ${correlations[*]}
do
  if [[ "$skipCorr" == "False" ]] || [[ "$skipCorr" == "F" ]] 
  then
    cp $goiList GOIs_and_Correlatedgenes_list.txt     
  else  
    sed -i -E "s/[0-9]\\.[0-9]+/$correlation/" Rscripts/siftSVscorr.R 
    Rscript Rscripts/siftSVscorr.R
  fi  
  mkdir -p ${correlation}corr
  mv -t ${correlation}corr GOIs_and_Correlatedgenes_list.txt Correlation_Table
  cd ${correlation}corr  
  grep -F -f GOIs_and_Correlatedgenes_list.txt ../annotatedgenome.gff > GOIs_and_Cs_annotatedgenome.bed
  # subset mRNA and CDS lines from GOIs_and_Cs_annotatedgenome.bed
  grep ID=mRNA GOIs_and_Cs_annotatedgenome.bed > GOIs_and_Cs_mRNA.bed 
  grep ID=CDS GOIs_and_Cs_annotatedgenome.bed > GOIs_and_Cs_CDS.bed
  grep ID=five_prime_UTR GOIs_and_Cs_annotatedgenome.bed > GOIs_and_Cs_five_prime_UTR.bed
  wc -l Correlation_Table *.* > wc_l_corr.txt
  for max_SVlength in ${max_SVlengths[*]}
  do	
    # **** Find the CDS that intersect SVs and then concatonate with full set of CDS of those genes/mRNAs
    bedtools window -a GOIs_and_Cs_CDS.bed -b ../${max_SVlength}maxSV/bedformat_SVs.bed -w 0 -u  > GOIsandCs_CDS_that_intersect_SVs.bed
    grep -E -e "Solyc[0-9]+g[0-9]+" GOIsandCs_CDS_that_intersect_SVs.bed -o | uniq > genes_with_CDS_overlapped_by_SVs.txt
    grep -F -f genes_with_CDS_overlapped_by_SVs.txt GOIs_and_Cs_CDS.bed > Full_CDS_sets.bed
    # which SVs intersect which CDS (and size of overlap)
    bedtools intersect -a  ../${max_SVlength}maxSV/bedformat_SVs.bed -b GOIs_and_Cs_CDS.bed -wo  > overlapping_SVs_and_GOIsandCs_CDS.bed
    # Process for finding SVs that are 2000 bp upstream promoter or intersect 5'UTR (and don't intersect CDS of any GOIs/Cs 
    bedtools intersect -a ../${max_SVlength}maxSV/bedformat_SVs.bed -b GOIs_and_Cs_CDS.bed -v > SVs_no_intersect_CDS.bed
    bedtools intersect -a SVs_no_intersect_CDS.bed -b GOIs_and_Cs_five_prime_UTR.bed -v > SVs_no_interesct_CDS_or_5_prime.bed 
    bedtools intersect -a SVs_no_intersect_CDS.bed -b GOIs_and_Cs_five_prime_UTR.bed -wa -wb > SVs_only_intersect_fives.bed
    # Put all into SV directory inside correlation directory 
    Rscript ../Rscripts/siftSVsmumath.R
    head -n 1 ../Anthesis.txt > Anthesis_CDS.txt 
    head -n 1 ../Meristems.txt > Meristems_CDS.txt 
    grep -F -f SVs_intersect_CDS_genelist.txt ../Anthesis.txt >> Anthesis_CDS.txt
    grep -F -f SVs_intersect_CDS_genelist.txt ../Meristems.txt >> Meristems_CDS.txt
    Rscript ../Rscripts/siftSVsCDSRelExpr.R $TM_SIM_FM $TM_SIM_FM_over_other_meristem_ratio $Anthesis_fiveDPA_tenDPA_psc $A510psc_over_A510_NONpsc $A510psc_over_NON_A510_psc $fiveDPA_psc_over_A_psc_and_tenDPA_psc
    ###HERE### Use this to include the SV wc -l from the outside folder and don't forget to move into mkdir 
    wc -l ../${max_SVlength}maxSV/bedformat_SVs.bed GOIsandCs_CDS_that_intersect_SVs.bed genes_with_CDS_overlapped_by_SVs.txt Full_CDS_sets.bed overlapping_SVs_and_GOIsandCs_CDS.bed SVs_no_intersect_CDS.bed SVs_no_interesct_CDS_or_5_prime.bed SVs_only_intersect_fives.bed SVs_intersect_CDS.txt SVs_intersect_CDS_genelist.txt Relevent_Expression_SVs_intersect_CDS.txt Anthesis_CDS.txt Meristems_CDS.txt > wc_l_SV.txt 
    mkdir -p ${max_SVlength}maxSV
    mv -t ${max_SVlength}maxSV SVs_no_interesct_CDS_or_5_prime.bed SVs_only_intersect_fives.bed Relevent_Expression_SVs_intersect_CDS.txt SVs_intersect_CDS.txt wc_l_SV.txt 
    cd ${max_SVlength}maxSV
    for windowSize in ${fivePrimeWindows[*]}
    do
      #changed A and B in bedtools call below (adjust pipeline accordingly)  
      bedtools window -a ../GOIs_and_Cs_five_prime_UTR.bed -b SVs_no_interesct_CDS_or_5_prime.bed -sw -r 0 -l $windowSize > SVs_upstream_of_genes.bed
      Rscript ../../Rscripts/siftSVsintregulatory.R 
      head -n 1 ../../Anthesis.txt > Anthesis_RE.txt 
      head -n 1 ../../Meristems.txt > Meristems_RE.txt 
      grep -F -f SVs_intersect_RE_genelist.txt ../../Anthesis.txt >> Anthesis_RE.txt
      grep -F -f SVs_intersect_RE_genelist.txt ../../Meristems.txt >> Meristems_RE.txt
      Rscript ../../Rscripts/siftSVsRERelExpr.R $TM_SIM_FM $TM_SIM_FM_over_other_meristem_ratio $Anthesis_fiveDPA_tenDPA_psc $A510psc_over_A510_NONpsc $A510psc_over_NON_A510_psc $fiveDPA_psc_over_A_psc_and_tenDPA_psc
      ###HERE### don't forget to move into mkdir 
      wc -l SVs_upstream_of_genes.bed SVs_intersect_RE.txt SVs_intersect_RE_genelist.txt Relevent_Expression_SVs_intersect_RE.txt Anthesis_RE.txt Meristems_RE.txt > wc_l_fiveprimewind.txt
      mkdir -p ${windowSize}fiveprimewindow
      mv -t ${windowSize}fiveprimewindow SVs_intersect_RE.txt Relevent_Expression_SVs_intersect_RE.txt wc_l_fiveprimewind.txt    
    done 
    rm SVs_no_interesct_CDS_or_5_prime.bed SVs_only_intersect_fives.bed SVs_upstream_of_genes.bed SVs_intersect_RE_genelist.txt Anthesis_RE.txt Meristems_RE.txt 
    cd .. 
  done 
  rm GOIsandCs_CDS_that_intersect_SVs.bed genes_with_CDS_overlapped_by_SVs.txt Full_CDS_sets.bed overlapping_SVs_and_GOIsandCs_CDS.bed SVs_no_intersect_CDS.bed SVs_intersect_CDS_genelist.txt Anthesis_CDS.txt Meristems_CDS.txt
  rm GOIs_and_Correlatedgenes_list.txt GOIs_and_Cs_annotatedgenome.bed GOIs_and_Cs_mRNA.bed GOIs_and_Cs_CDS.bed GOIs_and_Cs_five_prime_UTR.bed
  cd .. 
done

# rm maxSV and annotated genome intermediate files in 1st for loop 
rm -r *maxSV annotatedgenomeinter.txt annotatedgenome.gff
