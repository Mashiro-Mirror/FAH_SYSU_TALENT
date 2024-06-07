cjy_single_cell_proportion_test = function(single_cell_freq_input, column_one, column_two) {
  ## Force all columns into the same names
  names(single_cell_freq_input) = c(column_one, column_two, "Freq")
  ## Start proportion calculating
  nrows = nrow(single_cell_freq_input)
  sample_number = length(unique(single_cell_freq_input[,column_two]))
  for (i in seq(1, nrows, by = nrows/sample_number)) {
    sum = sum(single_cell_freq_input$Freq[i:(i+((nrows/sample_number)-1))])
    single_cell_freq_input$Proportion[i:(i+((nrows/sample_number)-1))] = single_cell_freq_input$Freq[i:(i+((nrows/sample_number)-1))] / sum
  }
  ## Return the result
  return(single_cell_freq_input)
}