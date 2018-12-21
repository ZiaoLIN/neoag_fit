#analyze outputs from LAST (version 9).
args=commandArgs(TRUE)
last_path = args[1]
neoantigen_path = args[2]
patient_name = args[3]

last_output <- read.table(last_path, sep = "\t", header =FALSE, fill = TRUE, skip=32)

#"""#0. check whether each WT peptide maps to the proteome exactly, i.e. a unique single locus.
#subject_id <- as.character(last_output$V2)
#wt_list <- c()
#for (i in 1:nrow(last_output)){
#  if (grepl("WT", subject_id[i])){
#    wt_list <- c(wt_list, i)
#  }
#}

#last_output_select <- last_output[wt_list, ]
#subject_id_select <- as.character(last_output_select$V2)
#wt_unique_mapping <- names(table(subject_id_select))[table(subject_id_select) ==1]

#summarize the names of wt and mut.
#mut_unique_mapping <- c()
#for (i in 1:length(wt_unique_mapping)){
#  current_wt <- wt_unique_mapping[i]
#  split_wt <- unlist(strsplit(current_wt, "[|]"))
#  current_mut <- paste(split_wt[1], "MUT", split_wt[3], split_wt[4], split_wt[5], sep = "|")
#  mut_unique_mapping <- c(mut_unique_mapping, current_mut)
#}

#total_unique_mapping <- c(wt_unique_mapping, mut_unique_mapping)
#subset the dataframe with selected wt and mut mapping rows.
#unique_mapping <- last_output[as.character(last_output$V2) %in% total_unique_mapping, ]"""


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

process_1_mut <- last_output[as.character(last_output$V2) %in% mut_exact_id, ]
process_1_table <- rbind(process_1, process_1_mut)

#"""#2. check whether MUT and WT map to the same locus.
#same_locus_id <- c()
#for (i in 1:length(wt_exact_id)){
#  current_wt <- wt_exact_id[i]
#  current_mut <- mut_exact_id[i]
#  if (unique_mapping_1[as.character(unique_mapping_1$V2) == current_wt, "V1"] == unique_mapping_1[as.character(unique_mapping_1$V2) == current_mut, "V1"]){
#    same_locus_id <- c(same_locus_id, i)
#  }
#}
#
#wt_same_id <- wt_exact_id[same_locus_id]
#mut_same_id <- mut_exact_id[same_locus_id]
#
#total_mapping_2 <- c(wt_same_id, mut_same_id)
#unique_mapping_2 <- unique_mapping_1[as.character(unique_mapping_1$V2) %in% total_mapping_2, ]"""

#2. check whether MUT is only 1 mismatch.

process_1_mut <- process_1_table[as.character(process_1_table$V2) %in% mut_exact_id, ]
process_2_mut_exclude <- process_1_mut[process_1_mut$V5 != 1, ]
mut_exclude_id <- as.character(process_2_mut_exclude$V2)
mut_final_id <- as.character(process_1_mut$V2)[!as.character(process_1_mut$V2) %in% mut_exclude_id]

wt_final_id <- c()
for (i in 1:length(mut_final_id)){
  current_mut <- mut_final_id[i]
  split_mut <- unlist(strsplit(current_mut, "[|]"))
  current_wt <- paste(split_mut[1], "WT", split_mut[3], split_mut[4], split_mut[5], sep = "|")
  wt_final_id <- c(wt_final_id, current_wt)
}

process_2_mut_table <- process_1_mut[process_1_mut$V2 %in% mut_final_id, ]
process_2_wt_table <- process_1_table[process_1_table$V2 %in% wt_final_id, ]
process_2_table <- rbind(process_2_wt_table, process_2_mut_table)
#total_final_id <- c(mut_final_id, wt_final_id)


#3. check whether MUT maps to the exact locus as WT.
mut_final_3 <- c()
for (i in 1:length(mut_final_id)){
    current_mut <- mut_final_id[i]
    current_mut_row <- process_2_table[process_2_table$V2 == current_mut, ,drop = FALSE]
    
    split_mut <- unlist(strsplit(current_mut, "[|]"))
    current_wt <- paste(split_mut[1], "WT", split_mut[3], split_mut[4], split_mut[5], sep = "|")
    current_wt_row <- process_2_table[process_2_table$V2 == current_wt, ,drop = FALSE]
    count = 0
    for (j in 1:nrow(current_mut_row)){
        if (as.character(current_mut_row$V1)[j] %in% as.character(current_wt_row$V1)){
            count = count + 1
        }
    }
    if (count == nrow(current_mut_row)){
        mut_final_3 <- c(mut_final_3, current_mut)
    }
    
}

wt_final_3 <- c()
for (i in 1:length(mut_final_3)){
    current_mut <- mut_final_3[i]
    split_mut <- unlist(strsplit(current_mut, "[|]"))
    current_wt <- paste(split_mut[1], "WT", split_mut[3], split_mut[4], split_mut[5], sep = "|")
    wt_final_3 <- c(wt_final_3, current_wt)
}

total_final_id <- c(mut_final_3, wt_final_3)


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


#write.table(total_matrix_2, file = "/Users/LinZiao/Desktop/Getz_lab/neoantigen/filter_by_rules_292T-Tumor-SM-BZ6K2.fasta", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


write.table(total_matrix_2, file = paste(str(patient_name), "filtered.fasta", sep = "_"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


