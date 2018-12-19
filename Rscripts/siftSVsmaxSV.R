#generate bedformat with SVlength maximum defined 
source("SVSiftFun.R")

SVs_to_Bedformat(senderofSV="Visa",filename="master_SVs_fromSofiaVisa_201807.vcf",max_SVlength=5e6)
