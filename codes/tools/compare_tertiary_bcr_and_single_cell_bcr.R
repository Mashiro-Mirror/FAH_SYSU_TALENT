compare_tertiary_bcr_and_single_cell_bcr = function(sc.bcell.meta, ter_bcr_data,
                                                    nt_or_aa = "nt", 
                                                    similarity_threshold = .8){
  ##Read before formal running:
  ##There is a major difference between the tertiary bcr data and the single cell bcr data(10X)
  ##That is, the length of the CDR3 sequence is different.
  ##Regarding the nucleotide seqs, the tertiary bcr data is shorter than the single cell bcr data.
  ##To be more specific, in the tertiary bcr data, 3 bps are missing in the 3' end of the CDR3nt seqs 
  ##and the same as the head of the CDR3 seqs.
  ##In regards to the amino acid seqs, the CDR3aa is one aa shorter than the single cell bcr data and the head of the CDR3aa
  ##And one aa shorter than the single cell bcr data in the end of the CDR3aa.
  ##So, when comparing the tertiary bcr data and the single cell bcr data, we need to take the above difference into consideration.
  ##The following code is to compare the tertiary bcr data and the single cell bcr data
  ##And return the bins of the tertiary bcr data that shares the same CDR3aa with the single cell bcr data(based on some criteria).
  
  ##Trim the 10X bcr data to the same length as the tertiary bcr data
  if (nt_or_aa == "aa"){
    sc.bcell.meta$cdr3 = substr(sc.bcell.meta$cdr3, 2, nchar(sc.bcell.meta$cdr3)-1)
  } else if (nt_or_aa == "nt"){
    sc.bcell.meta$cdr3_nt = substr(sc.bcell.meta$cdr3_nt, 4, nchar(sc.bcell.meta$cdr3_nt)-3)
  }
  ##Check the tertiary bcr data, make sure that all records are "IGH"
  ter_bcr_data = ter_bcr_data[ter_bcr_data$locus == "IGH",]
  ##Define a function to calculate the similarity between two strings
  ##The similarity is defined as the number of two equal length strings divided by the length of the string
  string_similarity = function(str1, str2){
    # Convert strings to sets of characters
    set1 = strsplit(str1, split = "")[[1]]
    set2 = strsplit(str2, split = "")[[1]]
    # Calculate the number of characters in common, if the strings are of equal length
    # If the strings are of unequal length, return 0
    if (length(set1) != length(set2)) {
      return(0)
    }
    for (i in 1:length(set1)) {
      if (set1[i] == set2[i]) {
        set1[i] = 1
      } else {
        set1[i] = 0
      }
    }
    # Convert character vector to numeric
    set1 <- as.numeric(set1)
    return(sum(set1)/length(set1))
  }
  ##Since there are multiple records in the tertiary bcr data, every records should generate a similarity score
  ##And the highest score will be used to determine the similarity between the tertiary bcr data and the single cell bcr data
  if (nt_or_aa == "aa"){
    bcr_highest_similarity = c()
    for (ter_bcr in ter_bcr_data$junction_aa) {
      similarity = sapply(sc.bcell.meta$cdr3, string_similarity, str2 = ter_bcr)
      bcr_highest_similarity = c(bcr_highest_similarity, max(similarity))
    }
  } else if (nt_or_aa == "nt"){
    bcr_highest_similarity = c()
    for (ter_bcr in ter_bcr_data$cdr3) {
      similarity = sapply(sc.bcell.meta$cdr3_nt, string_similarity, str2 = ter_bcr)
      bcr_highest_similarity = c(bcr_highest_similarity, max(similarity))
    }
  }
  new_ter_bcr_data = ter_bcr_data
  new_ter_bcr_data$highest_similarity = bcr_highest_similarity
  ##Filter the new_ter_bcr_data based on the similarity score
  new_ter_bcr_data = new_ter_bcr_data[new_ter_bcr_data$highest_similarity >= similarity_threshold,]
  return(new_ter_bcr_data)
}
