# pcb_tools.r

#' convert_to_coelution_reporting_format
#' 
#' PCB-100, with flag C93 implies that PCB-100 coelutes with C93, we want a single concentration for  PCB-93+100 to avoid double counting
#'    
#'    This function converts the Test America Flag format to the classic pattern PCB-10+14#' input must contain:
#'    txtConstituent
#'    txtQual with Coelution format as "C([0-9]+)"
#' 
#' @param df 
#'
#' @return output with a new names"PCB" that includes the coeluting congeners with plus signs, like PCB-4+10 )
convert_to_coelution_reporting_format <- function(df){
  # creates a data frame that contains the txtConstituent Name, and coelutant if applicable
  m <- as.data.frame(cbind(df$txtConstituent, 
                           stringr::str_match(df$txtQual, pattern = "C([0-9]+)")), 
                     stringsAsFactors = F)
  # converts the colutant to a number
  m$V3 <- as.numeric(m$V3)
  # removes the txtconstituent prefix and converts to a number
  m$V4 <- as.numeric(stringr::str_replace(m$V1, pattern = "PCB-", replacement = ""))
  m$ORIGINAL <- df$txtConstituent
  m <- m[order(m$V4),]
  for(i in 1:dim(m)[1]){
    coelution <- m$V3[i]
    if (!is.na(coelution)){
      ind = which(m$V4 == coelution)
      m$V1[ind] <- paste(m$V1[ind], m$V4[i], sep = "+")
    }
  }
  table <- data.frame(ORIGINAL = m$ORIGINAL, PCB= m$V1 , stringsAsFactors = F)
  output <-  left_join(df, table, by = c("txtConstituent" = "ORIGINAL"))
  return(output)
}

#' process_pcb_df
#' this is a wrapper function for convert_to_coelution_reporting_format(), 
#  it also sorts congener in ascending order
#
#' @param df data.frame containing test america format PCB congener results 
#'
#' @return df4 data.frame
#' @example 
#'  list1 <- split(PCB, f =  PCB$txtSampleID)
#'  list2 <- lapply(list1, function(i) process(i))
#'  dPCB <- do.call(rbind, list2)
#' 
#' 
#' 
process_individual_df <- function(df){
  # call convert_to_coelution_reporting_format
  # inodrer to append PCB column wiht PCB-Z+X+Y
  df2 <- convert_to_coelution_reporting_format(df) # Produces new PCB label PCB-Z+X+Y
  df3 <- df2[!grepl(df2$txtQual, pattern = "C[0-9]+"),] # Drops all samples that were reported as a coelution with Z
  df3$SORTBY <- as.numeric(stringr::str_replace(df3$txtConstituent, pattern = "PCB-", replacement = "")) # extracts a numeric congener values to sort by
  df4 <- arrange(df3, SORTBY) # sorts so PCB-20 comes before PCB-100 
  return(df4)
}

process_pcb_df <- function(df){
  list1 <- split(df, f =  df$txtSampleID)
  list2 <- lapply(list1, function(i) process_individual_df(i))
  dPCB <- do.call(rbind, list2)
  return(dPCB)
}
