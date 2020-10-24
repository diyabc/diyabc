context("human dataset (snp)")

setwd(file.path(datasetsdir,"human"))
humandir = list.files(".")
file.remove(humandir[! humandir %in% c("header.txt","human_snp_all22chr_maf5.snp","RNG_state_0000.bin")])

system2(generalexe,paste(cmdpre,"-R \"FST1;ML1\" -r 100 -g 50 -m -t 8"),stdout = FALSE)
data <- readRefTable(filename = "reftableRF.bin", header="headerRF.txt")
expect_equal(colnames(data$params),c("N1","N2","N3","N4","t1","t2","d3","Nbn3","d4","Nbn4","N34","t3","d34","Nbn34","t4","Na","ra","t11","t22","t33","t44"))
expect_equal(ncol(data$stats),12)
expect_equal(nrow(data$stats),100)
