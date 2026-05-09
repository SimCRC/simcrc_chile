# Identify strategies that we should remove because there
# Is another strategy that is effectively the same
# (e.g. COL4575q10 & COL4580q10, screening happens at 45,55,65,75, so we will remove one of them)
# Claudia 
# November 7, 2019

identify_strategies_to_remove <- function(modality = "COL",
                                          start_ages = c(45,50,55),
                                          stop_ages  = c(70,75,80,85),
                                          intervals  = c(5,10,15) 
){
  table <- NULL
  for(int in intervals){
    for(start in start_ages){
      for(stop in stop_ages){
        
        # Create a string with ages in which screening happens (full sequence)
        sequence <- paste(seq(start, stop, by = int), collapse = "_")
        
        # Store information 
        new_row <- c("start" = start, "stop" = stop, "interval" = int, "sequence" = sequence)
        table <- as.data.frame(rbind(table, new_row), stringsAsFactors = FALSE)
        
      }
    }
  }
  
  
  # Print the duplicated sequences
  df_dup_strategies <- table[duplicated(table$sequence) | duplicated(table$sequence, fromLast = TRUE),]
  
  # Keep one, and remove all other duplicated strategies
  df_keep <- df_dup_strategies[!duplicated(df_dup_strategies$sequence),]
  
  # List all the duplicate strategies to remove
  df_remove <- df_dup_strategies[duplicated(df_dup_strategies$sequence),]
  
  # For those strategies we are removing, store a strategy label
  df_remove$strategy_label <- paste0(modality, df_remove$start, "", df_remove$stop, "q", df_remove$interval)
  
  # Return the list of strategies
  l_strategies_to_remove <- df_remove$strategy_label 
  return(l_strategies_to_remove)
}



