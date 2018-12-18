#!/bin/bash

#instead of doing intermediates at the end do it along the way and delete intermediates files along the way too (>> to an appended intermedial wc -l list) !!!

 

correlations=(0.95 0.975 0.99)
max_SVlengths=(5e6 1e6 500000 1e5)
fivePrimeWindows=(10000 5000 1000)
#ReleventExpressionsVars--- MAKE SURE IN FORMAT OF R AND IN QUOTES--- order listed as arguments into Rscript ARE IMPORTANT (MAINTAIN SAME ORDER)  
TM_SIM_FM="c(20,15,10)" 
TM_SIM_FM_over_other_meristem_ratio="c(1.25)"
Anthesis_fiveDPA_tenDPA_psc="c(15,10)"
A510psc_over_A510_NONpsc="c(1.25)"
A510psc_over_NON_A510_psc="c(1.25)"
fiveDPA_psc_over_A_psc_and_tenDPA_psc="c(0.5)"

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
  sed -i -E "s/[0-9]\\.[0-9]+/$correlation/" Rscripts/siftSVscorr.R 
  Rscript Rscripts/siftSVscorr.R 
  mkdir -p ${correlation}corr
  mv -t ${correlation}corr GOIs_and_Correlatedgenes_list.txt Correlation_Table
  cd ${correlation}corr  
  grep -F -f GOIs_and_Correlatedgenes_list.txt ../annotatedgenome.gff > GOIs_and_Cs_annotatedgenome.bed
  # subset mRNA and CDS lines from GOIs_and_Cs_annotatedgenome.bed
  grep ID=mRNA GOIs_and_Cs_annotatedgenome.bed > GOIs_and_Cs_mRNA.bed 
  grep ID=CDS GOIs_and_Cs_annotatedgenome.bed > GOIs_and_Cs_CDS.bed
  grep ID=five_prime_UTR GOIs_and_Cs_annotatedgenome.bed > GOIs_and_Cs_five_prime_UTR.bed
  
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
    mkdir -p ${max_SVlength}maxSV
    mv -t ${max_SVlength}maxSV GOIsandCs_CDS_that_intersect_SVs.bed genes_with_CDS_overlapped_by_SVs.txt Full_CDS_sets.bed overlapping_SVs_and_GOIsandCs_CDS.bed       SVs_no_intersect_CDS.bed SVs_no_interesct_CDS_or_5_prime.bed SVs_intersect_CDS.txt SVs_intersect_CDS_genelist.txt SVs_only_intersect_fives.bed Relevent_Expression_SVs_intersect_CDS.txt
    cd ${max_SVlength}maxSV
    for windowSize in ${fivePrimeWindows[*]}
    do
      bedtools window -a SVs_no_interesct_CDS_or_5_prime.bed -b ../GOIs_and_Cs_five_prime_UTR.bed -sw -r 0 -l $windowSize > SVs_upstream_of_genes.bed
      Rscript ../../Rscripts/siftSVsintregulatory.R 
      head -n 1 ../../Anthesis.txt > Anthesis_RE.txt 
      head -n 1 ../../Meristems.txt > Meristems_RE.txt 
      grep -F -f SVs_intersect_RE_genelist.txt ../../Anthesis.txt >> Anthesis_RE.txt
      grep -F -f SVs_intersect_RE_genelist.txt ../../Meristems.txt >> Meristems_RE.txt
      Rscript ../../Rscripts/siftSVsRERelExpr.R $TM_SIM_FM $TM_SIM_FM_over_other_meristem_ratio $Anthesis_fiveDPA_tenDPA_psc $A510psc_over_A510_NONpsc $A510psc_over_NON_A510_psc $fiveDPA_psc_over_A_psc_and_tenDPA_psc
      mkdir -p ${windowSize}fiveprimewindow
      mv -t ${windowSize}fiveprimewindow SVs_upstream_of_genes.bed SVs_intersect_RE.txt SVs_intersect_RE_genelist.txt Relevent_Expression_SVs_intersect_RE.txt     
    done 
    cd .. 
  done 
  cd .. 
done
