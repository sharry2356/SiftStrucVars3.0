# SiftStrucVars3.0
multi-parameter set pipeline

#For biological context and thorough explanation of pipeline's intricacies, please review associated poster (created early in project's development) and ppt provided in this repository (*.pdf and *.pptx respectively)

Essentially, SiftSVFunc.R is the source code and the ./Rscripts folder's files are calling functions from that source code. The ./Rscripts 
are linked together in a a shell script (NEWsiftSVMaster.sh) that loops over certain Rscripts to enable multi-paramter runs and outputs. 
The shell script also serves to utilize the bedtools command line tool which is an essential component of the pipeline. 

These two remaining R scripts in the home wd (quickTable.R and analyzeFullTable.R) are not part of main pipeline but can be run afterward
to help interpret the main pipeline's output.  quickTable.R condenses and summarizes the main pipeline's output into 1 table from the many 
files generated in a multi-parameter run of the pipeline. Furthermore, analyzeFullTable.R essentially pieces together "quicktables" from 
different runs from different file inputs and isolates consistent patterns accross all the runs.  


