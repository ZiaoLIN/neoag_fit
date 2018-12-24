#The script is intended to map neoantigen_fitness output peptides back to the mutations that generate them.

args=commandArgs(TRUE)
neoantigen_fitness_output_path = args[1]
neoantigen_netMHC_output_path = args[2]
mutation_input_path = args[3]
patient_name = args[4]


#neoantigen_fitness_output <- read.table("/Users/LinZiao/Desktop/Getz_lab/neoantigen/SupplementaryDataFile7/output_new/neoantigen_fitness_292T.txt", sep = "\t", header =TRUE, fill = TRUE, stringsAsFactors = FALSE)
#neoantigen_netMHC_output <- read.table("/Users/LinZiao/Desktop/Getz_lab/neoantigen/292T-TP-NB-SM-BZ6K2-SM-AYGX7.combined.all.binders.txt.annot.txt.clean.txt", sep = "\t", header =TRUE, fill = TRUE, stringsAsFactors = FALSE)
#mutation_input <- read.table("/Users/LinZiao/Desktop/Getz_lab/neoantigen/292T-TP-NB-SM-BZ6K2-SM-AYGX7.validated.maf", sep = "\t", header =TRUE, fill = TRUE, quote = "", stringsAsFactors = FALSE)


neoantigen_fitness_output <- read.table(neoantigen_fitness_output_path, sep = "\t", header =TRUE, fill = TRUE, stringsAsFactors = FALSE)
neoantigen_netMHC_output <- read.table(neoantigen_netMHC_output_path, sep = "\t", header =TRUE, fill = TRUE, stringsAsFactors = FALSE)
mutation_input <- read.table(mutation_input_path, sep = "\t", header =TRUE, fill = TRUE, quote = "", stringsAsFactors = FALSE)


gene_mutation_list <- c()
for (i in 1:nrow(neoantigen_fitness_output)){
  split_mutation <- unlist(strsplit(neoantigen_fitness_output[i, ]$Mutation, "_"))
  gene_mutation <- paste(split_mutation[1], split_mutation[2], sep = "_")
  gene_mutation_list <- c(gene_mutation_list, gene_mutation)
}

neoantigen_fitness_output$gene_mutation <- gene_mutation_list
neoantigen_netMHC_output$gene_protein_change <- paste(neoantigen_netMHC_output$gene, neoantigen_netMHC_output$protein_change, sep = "_")

#match back to netMHC_output to get $cdna_change.
cdna_change_list <- c()
for (i in 1:length(gene_mutation_list)){
  gene_mutation <- gene_mutation_list[i]
  cdna_change <- unique(neoantigen_netMHC_output[neoantigen_netMHC_output$gene_protein_change == gene_mutation, ]$cdna_change)
  cdna_array <- unlist(strsplit(cdna_change, ">"))
  last_character <- substr(cdna_array[1], nchar(cdna_array[1]), nchar(cdna_array[1]))
  new_cdna <- paste(last_character, cdna_array[2], sep = ">")
  cdna_change_list <- c(cdna_change_list, new_cdna)
}


neoantigen_fitness_output$cdna_change <- cdna_change_list

gene_list <- c()
for (i in 1:nrow(neoantigen_fitness_output)){
  split_mutation <- unlist(strsplit(neoantigen_fitness_output[i, ]$Mutation, "_"))
  gene <- split_mutation[1]
  gene_list <- c(gene_list, gene)
}
neoantigen_fitness_output$gene <- gene_list

neoantigen_fitness_output$gene_cdna_change <- paste(neoantigen_fitness_output$gene, neoantigen_fitness_output$cdna_change, sep = "_")


mutation_cdna_change <- c()
for (i in 1:nrow(mutation_input)){
  mutation_cdna <- mutation_input$cDNA_Change[i]
  cdna_array <- unlist(strsplit(mutation_cdna, ">"))
  last_character <- substr(cdna_array[1], nchar(cdna_array[1]), nchar(cdna_array[1]))
  new_cdna <- paste(last_character, cdna_array[2], sep = ">")
  mutation_cdna_change <- c(mutation_cdna_change, new_cdna)
}

mutation_input$new_cdna_change <- mutation_cdna_change
mutation_input$gene_cdna_change <- paste(mutation_input$Hugo_Symbol, mutation_input$new_cdna_change, sep = "_")

chrom_list <- c()
start_position_list <- c()
end_position_list <- c()
ref_allele_list <- c()
tumor_seq_1_list <- c()
tumor_seq_2_list <- c()
cDNA_Change_list <- c()
#add columns of chrom/position/ref/alt back to neoantigen_fitness_output table by cdna column.
for (i in 1:nrow(neoantigen_fitness_output)){
  current_gene_mutation <- neoantigen_fitness_output$gene_cdna_change[i]
  if (sum(mutation_input$gene_cdna_change == current_gene_mutation) == 0){
    print (i)
    chrom <- " "
    Start_position <- " "
    end_position <- " "
    ref_allele <- " "
    tumor_seq_allele <- " "
    tumor_seq_allele_2 <- " "
    cDNA_Change <- " "
    chrom_list <- c(chrom_list, chrom)
    start_position_list <- c(start_position_list, Start_position)
    end_position_list <- c(end_position_list, end_position)
    ref_allele_list <- c(ref_allele_list, ref_allele)
    tumor_seq_1_list <- c(tumor_seq_1_list, tumor_seq_allele)
    tumor_seq_2_list <- c(tumor_seq_2_list, tumor_seq_allele_2)
    cDNA_Change_list <- c(cDNA_Change_list, cDNA_Change)
  }
  else {
  chrom <- mutation_input[mutation_input$gene_cdna_change == current_gene_mutation, ]$Chromosome
  Start_position  <- mutation_input[mutation_input$gene_cdna_change == current_gene_mutation, ]$Start_position
  end_position  <- mutation_input[mutation_input$gene_cdna_change == current_gene_mutation, ]$End_position
  ref_allele  <- mutation_input[mutation_input$gene_cdna_change == current_gene_mutation, ]$Reference_Allele
  tumor_seq_allele  <- mutation_input[mutation_input$gene_cdna_change == current_gene_mutation, ]$Tumor_Seq_Allele1
  tumor_seq_allele_2  <- mutation_input[mutation_input$gene_cdna_change == current_gene_mutation, ]$Tumor_Seq_Allele2
  #cDNA_Change <- mutation_input[mutation_input$gene_cdna_change == current_gene_mutation, ]$cDNA_Change
  chrom_list <- c(chrom_list, chrom)
  start_position_list <- c(start_position_list, Start_position)
  end_position_list <- c(end_position_list, end_position)
  ref_allele_list <- c(ref_allele_list, ref_allele)
  tumor_seq_1_list <- c(tumor_seq_1_list, tumor_seq_allele)
  tumor_seq_2_list <- c(tumor_seq_2_list, tumor_seq_allele_2)
  #cDNA_Change_list <- c(cDNA_Change_list, cDNA_Change)
}
}

neoantigen_fitness_output$chrom <- chrom_list
neoantigen_fitness_output$start_position <- start_position_list
neoantigen_fitness_output$end_position <- end_position_list
neoantigen_fitness_output$ref_allele <- ref_allele_list
neoantigen_fitness_output$tumor_seq_1 <- tumor_seq_1_list
neoantigen_fitness_output$tumor_seq_2 <- tumor_seq_2_list


write.table(neoantigen_fitness_output, file = paste(as.character("neoantigens_map_mutations"), paste(as.character(patient_name), ".txt", sep = ""), sep = "_"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)






