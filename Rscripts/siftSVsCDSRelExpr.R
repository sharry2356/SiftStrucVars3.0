#Relevent Expression for SV intersect CDS genes
source("../SVSiftFun.R")

Relevent_Expression_SVs_intersect_CDS<-Rel_Expr_Bash_Loop(Expression_Analysis_Func(ReadFiles(Meristems_file = "Meristems_CDS.txt", 
                                                        Anthesis_file = "Anthesis_CDS.txt"),zero.rm = T))
write.table(Relevent_Expression_SVs_intersect_CDS, "Relevent_Expression_SVs_intersect_CDS.txt",
            sep = "\t")
SVs_intersect_CDS<-AddReleventExpressiontoSVs("Relevent_Expression_SVs_intersect_CDS.txt","SVs_intersect_CDS.txt")
write.table(SVs_intersect_CDS,"SVs_intersect_CDS.txt", sep = "\t")
