intramutate = function(df, orig_search_list, replacement_list){
  
  for (i in 1:length(replacement_list)){
    
    df = case_when(
      df == orig_search_list[i] ~ replacement_list[i],
      df != orig_search_list[i] ~ df
    )
  }
  
  return(df)
}
