context("sequences + microsat")

setwd(file.path(datasetsdir,"sequences-mut"))
sequencesdir = list.files(".")
file.remove(sequencesdir[! sequencesdir %in% c("header.txt","4pops_11loci_1mito.txt","RNG_state_0000.bin")])

system2(generalexe,paste(cmdpre,"-R \"ALL\" -r 100 -g 50 -m -t 8"),stdout = FALSE)
data <- readRefTable(filename = "reftableRF.bin", header="headerRF.txt")
expect_identical(colnames(data$params),c("N1","N2","N3","N4","ti4","DB4","NF4","t3","DB3","NF3","t2","DB2","NF2","tg","tig2","t4","r","tg2","µmic_1","pmic_1","snimic_1","µseq_2","k1seq_2"))
expect_equal(ncol(data$stats),150)
expect_equal(nrow(data$stats),100)
