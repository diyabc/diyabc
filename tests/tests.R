library(testthat)

args = commandArgs(trailingOnly=TRUE)
env <- test_env()
generalexe=args[1]

generalexe <- normalizePath(generalexe)
if (Sys.info()[['sysname']] == "Windows") {
    cmdpre <- "-p .\\ "
} else {
    cmdpre <- "-p ./ "
}

projectdir <- getwd()
source(file.path(projectdir,"readReftable.R"))
datasetsdir <- file.path(projectdir,"datasets")

# test_file("diyabc.R",reporter="progress", env = env)
test_dir("testsdir",reporter="progress", env = env)
