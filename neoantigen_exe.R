#2018.12.18
#neoantigen fasta generation.
args=commandArgs(TRUE)
tissue_path = args[1]
patient_name = args[2]

neoantigen <- read.table(tissue_path, sep = "\t", header =TRUE, fill = TRUE)
#affinity_mut < 500
neoantigen_affinity_select <- neoantigen[neoantigen$affinity_mut < 500, ]
neoantigen_affinity_select <- neoantigen_affinity_select[neoantigen_affinity_select$neoORF_status=="mut|no", ]
#exclude peptides not 9-mers.
size_peptide <- nchar(as.character(neoantigen_affinity_select$peptide_wt), type = "chars", allowNA = FALSE, keepNA = NA)
neoantigen_affinity_select <- neoantigen_affinity_select[size_peptide==9, ]

MT <- neoantigen_affinity_select$peptide_mut
WT <- neoantigen_affinity_select$peptide_wt
sample_id <- neoantigen_affinity_select$sample
transcript_id <- neoantigen_affinity_select$transcript
protein_change <- neoantigen_affinity_select$protein_change
total_matrix <- matrix(nrow = 1, ncol = 1)
for (i in 1:nrow(neoantigen_affinity_select)){
  row1 <- paste(paste(as.character(">"), as.character(sample_id[i]), sep = ""), "WT", as.character(transcript_id[i]), as.character(protein_change[i]), as.character(i), sep = "|")
  row2 <- as.character(WT[i])
  row3 <- paste(paste(as.character(">"), as.character(sample_id[i]), sep = ""), "MUT", as.character(transcript_id[i]), as.character(protein_change[i]), as.character(i), sep = "|")
  row4 <- as.character(MT[i])
  submatrix <- rbind(row1, row2, row3, row4)
  total_matrix <- rbind(total_matrix, submatrix)
}
total_matrix_2 <- as.data.frame(total_matrix[2:nrow(total_matrix), 1:ncol(total_matrix)], drop = FALSE)

#write.table(total_matrix_2, file = "/Users/LinZiao/Desktop/Getz_lab/neoantigen/292T-Tumor-SM-BZ6K2.fasta", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(total_matrix_2, file = paste(str(patient_name), ".fasta", sep = ""), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
