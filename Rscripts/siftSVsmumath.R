#calculate affect of SVs on CDS of correlated genes 
source("../SVSiftFun.R")

MutationMath(Strip_to_Integers("Full_CDS_sets.bed","overlapping_SVs_and_GOIsandCs_CDS.bed"))