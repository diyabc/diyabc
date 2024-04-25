context("depindep dataset (-o)")

setwd(file.path(datasetsdir,"depindep"))
depindepdir = list.files(".")
file.remove(depindepdir[! depindepdir %in% c("headerRF.txt","DSIM_indep_dep_4pop_SNPind_10indsPerPop.snp","depindep_fast_1to10_param_S2.txt","yellow_depindep_fast_1to10_param_S2_ref_table.txt","maf.txt")])

system2(generalexe,paste(cmdpre,"-n \"t:1;c:1;s:1974;f:f\""),stdout = FALSE)
system2(generalexe,paste(cmdpre,"-o depindep_fast_1to10_param_S2.txt -i yellow_depindep_fast_1to10_param_S2.txt -g 10 -m -t 1"),stdout = FALSE)
data <- read.table("yellow_depindep_fast_1to10_param_S2.txt",header = FALSE)
dataref <- read.table("yellow_depindep_fast_1to10_param_S2_ref_table.txt",header = FALSE)
expect_equal(data,dataref)

file.remove(depindepdir[! depindepdir %in% c("headerRF.txt","DSIM_indep_dep_4pop_SNPind_10indsPerPop.snp","depindep_fast_1to10_param_S2.txt","yellow_depindep_fast_1to10_param_S2_ref_table.txt","maf.txt")])
