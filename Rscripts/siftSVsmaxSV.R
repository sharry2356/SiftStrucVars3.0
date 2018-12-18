#generate bedformat with SVlength maximum defined 
source("SVSiftFun.R")

SVs_to_Bedformat("Visa","master_SVs_fromSofiaVisa_201807.vcf",max_SVlength=1e5)
