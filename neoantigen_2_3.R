#combine_2_and_3.
#analyze outputs from LAST (version 9).

args=commandArgs(TRUE)
last_path = args[1]
neoantigen_path = args[2]
patient_name = args[3]
hla_path = args[4]

last_output <- read.table(last_path, sep = "\t", header =FALSE, fill = TRUE, skip=32)


#1. check whether WT maps to the proteome exactly.
#first select WT rows.
subject_id <- as.character(last_output$V2)
wt_list <- c()
for (i in 1:nrow(last_output)){
  if (grepl("WT", subject_id[i])){
    wt_list <- c(wt_list, i)
  }
}

last_output_wt <- last_output[wt_list, ]
#next check whether WT maps to the proteome exactly.

process_1 <- last_output_wt[last_output_wt$V3 == 100, ]
process_1_id <- as.character(process_1$V2)

mut_exact_id <- c()
for (i in 1:length(process_1_id)){
  current_wt <- process_1_id[i]
  split_wt <- unlist(strsplit(current_wt, "[|]"))
  current_mut <- paste(split_wt[1], "MUT", split_wt[3], split_wt[4], split_wt[5], sep = "|")
  mut_exact_id <- c(mut_exact_id, current_mut)
}

process_1_total <- c(process_1_id, mut_exact_id)
process_1_table <- last_output[as.character(last_output$V2) %in% process_1_total, ]

"""#2. check whether MUT and WT map to the same locus. 
same_locus_id <- c()
for (i in 1:length(wt_exact_id)){
current_wt <- wt_exact_id[i]
current_mut <- mut_exact_id[i]
if (unique_mapping_1[as.character(unique_mapping_1$V2) == current_wt, "V1"] == unique_mapping_1[as.character(unique_mapping_1$V2) == current_mut, "V1"]){
same_locus_id <- c(same_locus_id, i)
}
}

wt_same_id <- wt_exact_id[same_locus_id]
mut_same_id <- mut_exact_id[same_locus_id]

total_mapping_2 <- c(wt_same_id, mut_same_id)
unique_mapping_2 <- unique_mapping_1[as.character(unique_mapping_1$V2) %in% total_mapping_2, ]"""

#2. check whether MUT is only 1 mismatch.

process_1_mut <- process_1_table[as.character(process_1_table$V2) %in% mut_exact_id, ]
process_2_mut <- process_1_mut[process_1_mut$V5 == 1, ]
mut_final_id <- as.character(process_2_mut$V2)

wt_final_id <- c()
for (i in 1:length(mut_final_id)){
  current_mut <- mut_final_id[i]
  split_mut <- unlist(strsplit(current_mut, "[|]"))
  current_wt <- paste(split_mut[1], "WT", split_mut[3], split_mut[4], split_mut[5], sep = "|")
  wt_final_id <- c(wt_final_id, current_wt)
}

total_final_id <- c(mut_final_id, wt_final_id)


#not yet to check whether WT and MUT map to the same locus. How to define the same locus for multiple mapping??
#for now, the selected peptides are those with at least one mapping satisfiying WT = 100%, MT = 1 mismatch.


#read in sequence information.
neoantigen <- read.table(neoantigen_path, sep = "\t", header =TRUE, fill = TRUE)

neoantigen_affinity_select <- neoantigen[neoantigen$affinity_mut < 500, ]
neoantigen_affinity_select <- neoantigen_affinity_select[neoantigen_affinity_select$neoORF_status=="mut|no", ]
#exclude peptides not 9-mers.
size_peptide <- nchar(as.character(neoantigen_affinity_select$peptide_wt), type = "chars", allowNA = FALSE, keepNA = NA)
neoantigen_affinity_select <- neoantigen_affinity_select[size_peptide==9, ]



MT <- neoantigen_affinity_select$peptide_mut
WT <- neoantigen_affinity_select$peptide_wt
sample_id <- neoantigen_affinity_select$sample
transcript_id <- neoantigen_affinity_select$transcript
gene <- neoantigen_affinity_select$transcript['gene']
protein_change <- neoantigen_affinity_select$protein_change
total_matrix <- matrix(nrow = 1, ncol = 1)
count <- 0
for (i in 1:nrow(neoantigen_affinity_select)){
  row1 <- paste(as.character(sample_id[i]), "WT", as.character(transcript_id[i]), as.character(protein_change[i]), as.character(i), sep = "|")
  if (row1 %in% total_final_id){
    count <- count + 1
    row1 <- paste(paste(as.character(">"), as.character(sample_id[i]), sep = ""), "WT", count, paste(as.character(gene[i]), as.character(protein_change[i]), sample_id[i], sep = "_"), sep = "|")
    row2 <- as.character(WT[i])
    row3 <- paste(paste(as.character(">"), as.character(sample_id[i]), sep = ""), "MUT", count, paste(as.character(gene[i]), as.character(protein_change[i]), sample_id[i], sep = "_"), sep = "|")
    row4 <- as.character(MT[i])
    submatrix <- rbind(row1, row2, row3, row4)
    total_matrix <- rbind(total_matrix, submatrix)
  }
}
total_matrix_2 <- as.data.frame(total_matrix[2:nrow(total_matrix), 1:ncol(total_matrix)], drop = FALSE)


write.table(total_matrix_2, file = paste(str(patient_name), "filtered.fasta", sep = "_"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


#linux script to make blastdb
#in cga2
#makeblastdb -in iedb.fasta -input_type fasta -dbtype prot -title iedb -parse_seqids -out iedb.fasta

#convert netMHC to input txt format.
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


write.table(final_input_2, file = paste(str("neoantigens"), paste(patient_name, ".txt", sep = ""), sep = "_"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

