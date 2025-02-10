## remove all existing variables
rm(list=ls())

## load input parameter
## args 1: model_dir
## args 2: model_family_cutoff
## args 3: model_genus_cutoff
## args 4: model_species_cutoff
## args 5: input_plasmid_matrix
## args 6: output_dir

input_args <- commandArgs(trailingOnly = TRUE)

## load model directory
model_f_dir <- paste(input_args[1], "model_Family/", sep = "")
model_g_dir <- paste(input_args[1], "model_Genus/", sep = "")
model_s_dir <- paste(input_args[1], "model_Species/", sep = "")

## load model cutoff
input_family <- read.table(input_args[2], sep = ",", header = TRUE, quote = "")
input_genus <- read.table(input_args[3], sep = ",", header = TRUE, quote = "")
input_species <- read.table(input_args[4], sep = ",", header = TRUE, quote = "")

## load input plasmid reference distance matrix
input_word <- read.table(input_args[5], sep = "\t", header = FALSE, quote = "")
colnames(input_word) <- c("ID", rep(paste("P", seq(1, 923), sep = "")))
rownames(input_word) <- input_word$ID
input_word <- input_word[, -1]

## create result score data frame
prob_family <- data.frame(matrix(rep(0, nrow(input_word) * nrow(input_family)), nrow = nrow(input_word), ncol = nrow(input_family), byrow = TRUE))
prob_genus <- data.frame(matrix(rep(0, nrow(input_word) * nrow(input_genus)), nrow = nrow(input_word), ncol = nrow(input_genus), byrow = TRUE))
prob_species <- data.frame(matrix(rep(0, nrow(input_word) * nrow(input_species)), nrow = nrow(input_word), ncol = nrow(input_species), byrow = TRUE))

## define function: return predict result
library(e1071)
model_prob_get <- function(svm_model, input_word){
  pred <- predict(svm_model, input_word, decision.values = TRUE)
  result <- data.frame(ID = rownames(input_word), distance = rep(0, nrow(input_word)))
  result$distance <- as.numeric(attr(pred, "decision.values"))
  return(result)
}

## define function: return final predict result
model_prob_predict <- function(input_word, model_dir, cutoff_list, result_class){
  for (rown in 1:nrow(cutoff_list)) {
    unit_choose <- cutoff_list[rown, 1]
    unit_choose <- gsub('[ ]', '_', unit_choose)
    unit_dir <- paste(model_dir, unit_choose, sep = "")
    unit_model <- paste(unit_dir, "_model.rds", sep = "")
    s_model <- readRDS(unit_model)
    s_result <- model_prob_get(s_model, input_word)
    result_class[, rown] <- s_result$distance
    colnames(result_class)[rown] <- unit_choose
    rm(rown, unit_choose, unit_dir, s_model, s_result)
  }
  rownames(result_class) <- rownames(input_word)
  return(result_class)
}

## get probabilities
HR_Prob_Family <- model_prob_predict(input_word, model_f_dir, input_family, prob_family)
HR_Prob_Genus <- model_prob_predict(input_word, model_g_dir, input_genus, prob_genus)
HR_Prob_Species <- model_prob_predict(input_word, model_s_dir, input_species, prob_species)

## output prob result
write.table(HR_Prob_Family, file = paste(input_args[6], "prob_result_family.tsv", sep = ""), sep = ",", col.names = TRUE, row.names = TRUE, quote = FALSE)
write.table(HR_Prob_Genus, file = paste(input_args[6], "prob_result_genus.tsv", sep = ""), sep = ",", col.names = TRUE, row.names = TRUE, quote = FALSE)
write.table(HR_Prob_Species, file = paste(input_args[6], "prob_result_species.tsv", sep = ""), sep = ",", col.names = TRUE, row.names = TRUE, quote = FALSE)


## create result class data frame
result_family <- data.frame(matrix(rep("No", nrow(input_word) * nrow(input_family)), nrow = nrow(input_word), ncol = nrow(input_family), byrow = TRUE))
result_genus <- data.frame(matrix(rep("No", nrow(input_word) * nrow(input_genus)), nrow = nrow(input_word), ncol = nrow(input_genus), byrow = TRUE))
result_species <- data.frame(matrix(rep("No", nrow(input_word) * nrow(input_species)), nrow = nrow(input_word), ncol = nrow(input_species), byrow = TRUE))

## define function: return class
model_result_get <- function(unit_choose, input_prob, cutoff){
  result <- data.frame(ID = rownames(input_prob), distance = rep(0, nrow(input_prob)), class = rep("No", nrow(input_prob)))
  result$distance <- input_prob[, c(unit_choose)]
  result$class[result$distance > cutoff] <- "Yes"
  return(result)
}

## define function: return final predict result
model_predict <- function(input_prob, cutoff_list, result_class){
  for (rown in 1:nrow(cutoff_list)) {
    unit_choose <- cutoff_list[rown, 1]
    unit_choose <- gsub('[ ]', '_', unit_choose)
    s_result <- model_result_get(unit_choose, input_prob, cutoff_list[rown, 2])
    result_class[, rown] <- s_result$class
    colnames(result_class)[rown] <- unit_choose
    rm(rown, unit_choose, s_result)
  }
  result_class$ID <- rownames(input_prob)
  result_class <- result_class[, c(ncol(result_class), 1:(ncol(result_class)-1)) ]
  return(result_class)
}

## load R package and predict
HR_Family <- model_predict(HR_Prob_Family, input_family, result_family)
HR_Genus <- model_predict(HR_Prob_Genus, input_genus, result_genus)
HR_Species <- model_predict(HR_Prob_Species, input_species, result_species)

## output result
write.table(HR_Family, file = paste(input_args[6], "result_family.tsv", sep = ""), sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(HR_Genus, file = paste(input_args[6], "result_genus.tsv", sep = ""), sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(HR_Species, file = paste(input_args[6], "result_species.tsv", sep = ""), sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
