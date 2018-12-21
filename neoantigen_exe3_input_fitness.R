#linux script to make blastdb
#in cga2
#makeblastdb -in iedb.fasta -input_type fasta -dbtype prot -title iedb -parse_seqids -out iedb.fasta

#convert netMHC to input txt format.


args=commandArgs(TRUE)
neoantigen_path = args[1]
hla_path = args[2]
patient_name = args[3]

neoantigen <- read.table(neoantigen_path, sep = "\t", header =TRUE, fill = TRUE)
neoantigen$mutation_id <- paste(neoantigen$gene, neoantigen$protein_change, neoantigen$sample, sep = "_")

#affinity_mut < 500
neoantigen_affinity_select <- neoantigen[neoantigen$affinity_mut < 500, ]
neoantigen_affinity_select <- neoantigen_affinity_select[neoantigen_affinity_select$neoORF_status=="mut|no", ]
#exclude peptides not 9-mers.
size_peptide <- nchar(as.character(neoantigen_affinity_select$peptide_wt), type = "chars", allowNA = FALSE, keepNA = NA)
neoantigen_affinity_select <- neoantigen_affinity_select[size_peptide==9, ]


neoantigen_affinity_select$id <- 1:nrow(neoantigen_affinity_select)
subset_neoantigen <- neoantigen_affinity_select[, c("sample", "mutation_id", "peptide_wt", "peptide_mut", "affinity_wt", "affinity_mut", "hla")]
subset_neoantigen$id <- paste(subset_neoantigen$sample, neoantigen_affinity_select$protein_change, neoantigen_affinity_select$id, sep = "_")
subset_neoantigen$combined_id <- paste(subset_neoantigen$sample, neoantigen_affinity_select$transcript, neoantigen_affinity_select$protein_change, neoantigen_affinity_select$id, sep = "|")

final_subset_id <- c()
for (i in 1:length(wt_final_id)){
  current_final <- wt_final_id[i]
  current_subset <- unlist(strsplit(current_final, "[|]"))
  current_subset_combine <- paste(current_subset[1], current_subset[3], current_subset[4], current_subset[5], sep = "|")
  final_subset_id <- c(final_subset_id, current_subset_combine)
}

#subset_neoantigen
final_subset_neoantigen <- subset_neoantigen[subset_neoantigen$combined_id %in% final_subset_id, ]
#read in HLA table.
hla_table <- read.table(hla_path, sep = "\t", header =FALSE, fill = TRUE)

combine_hla <- as.character(hla_table$V1[1])
for (i in 2:nrow(hla_table)){
  current_hla <- as.character(hla_table$V1)[i]
  combine_hla <- paste(combine_hla, current_hla, sep = ",")
  
}

hla_column <- rep(combine_hla, nrow(final_subset_neoantigen))
final_subset_neoantigen$HLA <- hla_column
final_subset_neoantigen$ID <- 1:nrow(final_subset_neoantigen)

#final input.
final_input <- final_subset_neoantigen[, c("ID", "mutation_id", "combined_id", "sample", "peptide_wt", "peptide_mut", "hla", "affinity_wt", "affinity_mut", "HLA")]
final_input$hla2 <- gsub(":", "", final_input$hla)
final_input$HLA2 <- gsub(":", "", final_input$HLA)
final_input$HLA2 <- gsub("HLA-", "", final_input$HLA2)
chop_score <- rep(1, nrow(final_input))
final_input$chop_score <- chop_score
final_input$combined_id <- NULL
final_input$hla <- NULL
final_input$HLA <- NULL

final_input_2 <- final_input[, c("ID", "mutation_id", "sample", "peptide_wt", "peptide_mut", "hla2", "affinity_wt", "affinity_mut", "HLA2", "chop_score")]
colnames(final_input_2) <- c("ID", "MUTATION_ID", "Sample", "WT.Peptide", "MT.Peptide", "MT.Allele", "WT.Score", "MT.Score", "HLA", "chop_score")
final_input_2$WT.Score <- as.integer(final_input_2$WT.Score)
final_input_2$MT.Score <- as.integer(final_input_2$MT.Score)


#write.table(final_input_2, file = "/Users/LinZiao/Desktop/Getz_lab/neoantigen/SupplementaryDataFile7/input_new/neoantigens_292T.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


write.table(final_input_2, file = paste(str("neoantigens"), paste(patient_name, ".txt", sep = ""), sep = "_"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
