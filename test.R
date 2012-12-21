peptide.strings <- c(
    "AAAGAEAGKATTEEQ",
    "AAEMEEALKGLPIRY",
    "AAFTSSSKAATAKAP",
    "AAGMEAQFLYLYALI",
    "AALAAAAGVPPADKY",
    "AAMLFCAVVIIGVLH",
    "AANVMAASLRKAGKS",
    "AAPLSWSKDIYNYME",
    "AARGWAAHRARANES",
    "ACVKDLVSKYLADNE",
    "EAHACQINSDQKFVD",
    "EAIKGGRHLIFCHSKKK",
    "EAKITMLTNGQCQNI")

source("computingKernel.R")
load("BLOSUM62.2.Rdata")
K3 <- computingKernel(peptide.strings, BLOSUM62.2, 0.11387, "GK3W")
K3.hat <- computingKernel(peptide.strings, BLOSUM62.2, 0.11387)
K3t <- computingKernel2(peptide.strings, BLOSUM62.2, 0.11387, "GK3W")
K3.hatt <- computingKernel2(peptide.strings, BLOSUM62.2, 0.11387)
