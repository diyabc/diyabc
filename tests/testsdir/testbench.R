context("multithreaded performance")

setwd(file.path(datasetsdir,"bench"))
indseqdir = list.files(".")
file.remove(indseqdir[! indseqdir %in% c("header.txt","INDSNP_sim_dataset_4POP_001.snp","RNG_state_0000.bin")])

system2(generalexe,paste(cmdpre,"-R \"ALL\" -r 1 -g 1 -m -t 8"),stdout = FALSE)
result1 = system.time(system2(generalexe,paste(cmdpre,"-R \"ALL\" -r 101 -g 20 -m -t 8"),stdout = FALSE))["elapsed"]
result2 = system.time(system2(generalexe,paste(cmdpre,"-R \"ALL\" -r 201 -g 20 -m -t 1"),stdout = FALSE))["elapsed"]
expect_gt(result2,result1*1.5)
