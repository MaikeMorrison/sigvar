ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

A <- c(88,85,28,68, 141, 146,44, 116,88,98,33,94,85, 111,21,95,88,85,28,68, 141, 146,
       44, 116,88,98,33,94, 85, 111,21,95,88,85,28,68, 141, 146,44, 116,88,98,33,94, 
       85, 111,21,95,29,57,78,46,38, 108, 141,86, 40,68,97,68,33,80,87, 115,29,57,
       78,46,38, 108, 141,86,40,68,97,68,33,80,87, 115,29,57,78,46,
       38, 108, 141,86,40,68,97,68,33,80,87, 115)
names(A) = c("ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC", "GCG", "GCT","TCA","TCC" ,
             "TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT", "GCA","GCC","GCG","GCT" ,
             "TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC",
             "GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT", "CTA","CTC","CTG","CTT","GTA",
             "GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG", "ATT","CTA","CTC","CTG" ,
             "CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA", "ATC","ATG","ATT","CTA",
             "CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT")

A_error = A
A_error[2] = 37

test_that("get_SBS96_spectrum works", {
  expect_no_warning(get_SBS96_spectrum(transcript = "ENST00000269305.9", ref_genome=ref_genome))
  expect_true(all(get_SBS96_spectrum(transcript = "ENST00000269305.9", ref_genome=ref_genome)==A) )
  expect_error(get_SBS96_spectrum(A_error))
  expect_false( all(get_SBS96_spectrum(transcript = "ENST00000269305.9", ref_genome=ref_genome)==A_error) )
})
