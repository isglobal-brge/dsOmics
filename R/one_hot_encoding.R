#' Create dummy variables
#' 
#' @details Function adapted from https://stackoverflow.com/a/66817515
#'
#' @param df Input data frame
#' @param columns Columns to create dummies
#'
#' @return
#'
#' @examples
.one_hot_encoding = function(df, columns="season"){
  # create a copy of the original data.frame for not modifying the original
  df = cbind(df)
  # convert the columns to vector in case it is a string
  columns = c(columns)
  # for each variable perform the One hot encoding
  for (column in columns){
    unique_values = sort(unique(df[column])[,column])
    non_reference_values  = unique_values#[c(-1)] # the first element is going 
    # to be the reference by default
    for (value in non_reference_values){
      # the new dummy column name
      new_col_name = paste0(column,value)
      # create new dummy column for each value of the non_reference_values
      df[new_col_name] <- with(df, ifelse(df[,column] == value, 1, 0))
    }
    # delete the one hot encoded column
    df[column] = NULL
    
  }
  return(df)
}
