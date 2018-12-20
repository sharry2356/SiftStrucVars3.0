#correlation1
source("SVSiftFun.R")

AllCorrelations_IMPROVED(goiList_FILE="shapeGoIList.txt",expression_table_FILE="markGeneList_RPKM_OvateCoexpr.csv",c(2,3,4,5,6,7),0.99)
