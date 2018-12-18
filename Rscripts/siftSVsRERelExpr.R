source("../../SVSiftFun.R")

Relevent_Expression_SVs_intersect_RE<-Rel_Expr_Bash_Loop(Expression_Analysis_Func(ReadFiles(Meristems_file = "Meristems_RE.txt", 
                                                                                             Anthesis_file = "Anthesis_RE.txt"),zero.rm = T))
write.table(Relevent_Expression_SVs_intersect_RE, "Relevent_Expression_SVs_intersect_RE.txt",
            sep = "\t")
SVs_intersect_RE<-AddReleventExpressiontoSVs("Relevent_Expression_SVs_intersect_RE.txt","SVs_intersect_RE.txt")
write.table(SVs_intersect_RE,"SVs_intersect_RE.txt", sep = "\t")
