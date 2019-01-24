# BACKGROUND

# Laboratories report concentrations of individual chemicals; however, 
# environmental  regulations are often based on the total aggregate concentration of multiple 
# analytes of the same class 
# (e.g., Total PCBs, Total cPAHs, Total DDT + homologs)

# Rules and reporting requirements for summed parameters, particularly 
# when some results are below analytical method detection limits, 
# can be agency, task, or site-specific. 

# Summary concentrations can also depend on equivalent toxicity corrections.
# For instance, PAHs that are known to be carcinogenic (cPAHs) 
# are often summed with a correction factor to report cPAH concentration 
# equivalent to the potency of cPAHs—A benzo(a)pyrene (BaP).

# Total cPAHs—A benzo(a)pyrene (BaP) equivalent (BaPEq) concentration of a sample 
# is calculated by multiplying the carcinogenic PAHs (cPAHs) by their respective 
# potency equivalency factors (PEFs) and summing the resulting concentrations. 

# The functions in this library provide a basis for a general a reproducible
# workflow for producing  Electronic Data Deliverables (EDDs) 
# for these summary results


#' safe add
#' 
#' add NAs such that they are treated as zeros
#' 
#' @param a numeric or numeric vector 
#' @param b numeric or numeric vector
#'
#' @return aplusb the addition of a and b where NA is treated as zero
safe_add<-function(a,b){
  a[is.na(a)] <- 0
  b[is.na(b)] <- 0
  aplusb = a+b
  return(aplusb)
}

#' get_unique_chars
#' 
#' given any string return an sorted list of unique single charcters in alphabetical order
#'
#' @param x a string
#'
#' @return a list of unique single characters making up the input string
get_unique_chars         <- function(x) {sort(unique(strsplit(x, "")[[1]]))}


#' collapse_list_of_chars 
#' 
#' collapse a vector of strings toa single string
#' 
#' @param x a vector of charcters or numeric values
#'
#' @return a string of the characters in the input vector
collapse_list_of_chars   <- function(x) {paste(x, collapse = "")}

#' append_special_flag
#'
#' append a special character, such as a T flag, to a string
#' 
#' @param x a string
#' @param special a string (defaults to T for "Total") to add to the input string x
#'
#' @return a string with the appended special character
append_special_flag      <- function(x, special = "T") {paste(x,special, sep = "")}

#' remove_special_flag
#' 
#' removes all instances of a designated string from the input string
#'
#' @param x a string
#' @param special a substring to remove from the primary input string x
#'
#' @return a string with the offending special character
remove_special_flag      <- function(x, special = "NA") {
  stringr::str_replace_all(x, pattern =special, replacement = "")
}

#' return all all unique flags
#' 
#' from a vecotr of characters, return the unique character asa string, with instances of "U" and "NA" removed and "T" appended
#'
#' @param x a vector of character flags (ie. c("U","K","B"))
#'
#' @return a string containing all the unique character in the input string, with instances of "U" and "NA" removed and "T" appended
return_all_unique_flags <- function(x) {
  result <- collapse_list_of_chars(sort(unique(x))) %>%
    remove_special_flag("NA") %>%
    remove_special_flag("U") %>% 
    get_unique_chars() %>%
    collapse_list_of_chars() %>%
    append_special_flag("T")
  return(result) 
}

#' general_analytical_summary
#'
#' General purpose summary function for analytical results.
#' 
#' @param df data.frame object must contain columns:
#' CLASS, txtConstituent, txtSampleID, dblResult, dblLimit, txtQual
#' @param PDI boolean defaults to TRUE, if PDI is true "Calculated totals are 
#' the sum of all detected results and one half of the highest reporting
#' detection limit for the results reported as not detected (ND) in a sample"
#' (see page 3. | Appendix E: 2018 PDI Surface Sediment Database Description) 
#'
#' @return df.sums data.frame object containing the following outputs as columns:
#'  CLASS, txtSampleID - grouping variables
#'  sumResult_0ND      - summary treating all non-detect (ND) as a concentration of zero 
#'  sumResult_halfND   - summary treating all ND as a concentration of half the DL
#'  sumResult_fullND   - summary treating all ND as a concentration as the DL
#'  maxLimit           - the highest detection limit in the analyte set
#'  n                  - the number of analytes in the summary
#'  nd                 - the number of analytes with a non-detect (ND) result
#'  txtQual            - the summary flag (depends on all_unique_flags() function) 
#'  ND_zero            - if all results are ND returns the highest dection limit, othewise sumResult_0ND 
#'  ND_half            - if all results are ND returns the highest dection limit, otherwise sumResult_halfND 
#'  ND_full            - if all results are ND returns the highest dection limit, otherwise or sumResult_fullND 

# this is a version of the general analytical summary meant to comply strictly with 
# page 3. | Appendix E: 2018 PDI Surface Sediment Database Description
# "Calculated totals are the sum of all detected results and one half of the highest reporting
# detection limit for the results reported as not detected (ND) in a sample."

# Some datasets will include Potency Equivalent Factors. 
# If PEFs are included, then check that every analyte has a PEF
# Otherwise assign all analytes a PEF of 1

general_analytical_summary <- function(df, PDI = T){

  if (!"PEF" %in% names(df)){
    df$PEF <- 1
  }else{
    if(sum(is.na(df$PEF)) > 0){
      stop("Data is compatabile with general_analytical_summary method\n
           PEFs are defined for some analytes and not others.")
    } 
  }
  
  if (!"CLASS" %in% names(df)){
    df$CLASS <- "ALL"
  }else{
    if(sum(is.na(df$CLASS)) > 0){
      stop("Data is compatabile with general_analytical_summary method\n
           not all analytes are assigned a class.")
    } 
  }
  
  
  
  df.maxLimits <- df %>% select(CLASS,
                                txtConstituent,
                                txtSampleID,
                                dblLimit,
                                dblResult,
                                PEF) %>%
    mutate(DL = is.na(dblResult)) %>%                      # This labels ND samples as DL = TRUE
    group_by(CLASS, txtSampleID, DL) %>%                   # This groups all the ND samples together per CLASS,txtSampleID
    summarise(dblmaxLimitSample = max(dblLimit*PEF, na.rm= T), # maxLimit calculates highest limits among NDs, 
              totalNDs = n()) %>% 
    mutate(dblmaxLimitSample = ifelse(dblmaxLimitSample == -Inf, 0, dblmaxLimitSample)) %>%
    filter(DL == TRUE) %>% ungroup()  # this filters so that we only have maxLimit from the samples with atleaast one ND
  
  df.sums <- df %>% select(CLASS,
                           txtConstituent,
                           txtSampleID, 
                           dblResult, 
                           dblLimit, 
                           txtQual,
                           PEF) %>% 
    left_join(df.maxLimits) %>%
    # page 3. | Appendix E: 2018 PDI Surface Sediment Database Description
    # "Calculated totals are the sum of all detected results and one half of the highest reporting
    # detection limit for the results reported as not detected (ND) in a sample."
    mutate(dblResultzeroND = ifelse(is.na(dblResult), 0  *dblLimit, dblResult*PEF))  %>%
    mutate(dblResulthalfND = ifelse(is.na(dblResult), 0.5*dblLimit, dblResult*PEF))  %>%
    mutate(dblResultfullND = ifelse(is.na(dblResult), 1  *dblLimit, dblResult*PEF))  %>%
    group_by(CLASS, txtSampleID) %>%
    summarise(sumResult_0ND    = sum(dblResultzeroND, na.rm=TRUE), 
              sumResult_halfND = sum(dblResulthalfND, na.rm=TRUE),
              sumResult_fullND = sum(dblResultfullND, na.rm=TRUE),
              maxLimit         = max(dblLimit),
              maxNDLimit       = max(dblmaxLimitSample),
              n                = n(),
              nd               = sum(is.na(dblResult)),
              txtQual       = return_all_unique_flags(txtQual) ) %>%
    mutate( sumResult_halfPDI = safe_add(sumResult_0ND, 0.5*maxNDLimit)) %>%
    mutate(`ND_zero` = ifelse(sumResult_0ND == 0, maxLimit, sumResult_0ND), 
           `ND_half` = ifelse(sumResult_0ND == 0, maxLimit, sumResult_halfND),
           `ND_full` = ifelse(sumResult_0ND == 0, maxLimit, sumResult_fullND),
           `ND_PDI`  = ifelse(sumResult_0ND == 0, maxLimit, sumResult_halfPDI),
           `txtQual` = ifelse(n == nd, paste("U",txtQual, sep = ""),txtQual))
  return(df.sums)
  
  }



#' output_EDD_from_sums_and_ref
#'
#' @param x.sum a data.frame with column names
#'   "txtSampleID", 
#'   "sumResult_0ND", 
#'   "sumResult_halfND", 
#'   "sumResult_fullND", 
#'   "maxLimit", 
#'   "n", 
#'   "nd", 
#'   "ND_zero", 
#'   "ND_half", 
#'   "ND_full",
#'   "ND_PDI", 
#'   "txtQual", 
#'   "CLASS"
#' @param x.ref a data.frame with other sample information which will be joined. It must contain the following headers to avoid an error
#'    txtProjectName
#'    txtEvent
#'    txtLabID
#'    dtmSampleDate
#'    txtAnalysisMethod
#'    txtConstiuent
#'    txtPCode
#'    txtSampleID
#'    dtmExtractedDate
#'    dtmAnalyzedDate
#'    intDilution
#'    dblResult
#'    txtUnits
#'    txtQual
#'    maxLimit
#'    txtQA
#'    txtComments
#'
#' @return x3 - data.frame with summation results in EDD format
output_EDD_from_sums_and_ref<- function(x.sum, x.ref, my_signif = 2){
  
  ref_info <- x.ref %>% group_by(txtSampleID, CLASS) %>% slice(1) %>% 
    select(-txtConstituent, -txtCASID, -txtQual, -dblLimit, -dblResult)
  

  x1 <- x.sum %>% select(CLASS, txtSampleID, ND_zero, ND_half, ND_full, ND_PDI, txtQual, maxLimit) %>%
    gather(summation, dblResult, -CLASS, -txtSampleID, -txtQual, -maxLimit) %>% 
    mutate(txtConstiuent = paste0(CLASS, "_", summation)) %>% ungroup()
  
  x2 <- x1 %>% left_join(ref_info)
  
  x3 <- x2 %>% transmute(PRJ_NAME  = txtProjectName,
                         EVENT     = txtEvent,
                         LAB_ID    = txtLabID,
                         SAMPLE_DT = dtmSampleDate,
                         ANALYSIS_METH  = txtAnalysisMethod,
                         CONSTITUENT    = txtConstiuent,
                         CAS_ID         = NA,
                         P_CODE         = txtPCode,
                         SAMPLE_ID      = txtSampleID,
                         EXTRACTED_DATE = dtmExtractedDate,
                         ANALYZED_DT    = dtmAnalyzedDate,
                         DILUTION       = intDilution,
                         RESULT         = dblResult,
                         UNITS          = txtUnits,
                         QUAL           = txtQual,
                         LIMIT          = maxLimit,
                         QA             = txtQA,
                         COMMENTS       = txtComments)
  # apply significant figures to result
  x3$RESULT <- signif(x3$RESULT, my_signif)
  #change posix date format to MM/DD/YYY hh:mm character
  x3$SAMPLE_DT<-as.character(x3$SAMPLE_DT, format='%m/%d/%y %H:%M')
  x3$EXTRACTED_DATE<-as.character(x3$EXTRACTED_DATE, format='%m/%d/%y %H:%M')
  x3$ANALYZED_DT<-as.character(x3$ANALYZED_DT, format='%m/%d/%y %H:%M')
  #change NAs and NA characters to blanks
  x3[is.na(x3)] <- ""
  return(x3)
}




# Define selectors and assigners
selectorPCB <- function(x){ 
  y = filter(x, txtAnalysisMethod == "E1668A", grepl("PCB-*", txtConstituent), !grepl("T", txtQual))
  return(y)
}

selectorAroclor <- function(x){ 
  y = filter(x,txtAnalysisMethod == "SW8082A", !grepl(txtQual, pattern = "T"))
  return(y)
}

selectorDDx <- function(x){
  y = filter(x, txtConstituent %in% c("o,p'-DDD", 
                                      "4,4'-DDD", 
                                      "o,p-DDE",
                                      "4,4'-DDE", 
                                      "o,p'-DDT",  
                                      "4,4'-DDT"),
             !grepl(txtQual, pattern = "T"))
  return(y)
}


selectorCPAH <- function(x){
  y = filter(x, txtConstituent %in%  c("Benzo(a)anthracene",
                                       "Benzo(a)pyrene",
                                       "Benzo(b)fluoranthene",
                                       "Benzo(k)fluoranthene",
                                       "Chrysene",
                                       "Dibenzo(a,h)anthracene",
                                       "Indeno(1,2,3-cd)pyrene"), 
             !grepl("T", txtQual))
  return(y)
}

selectorLPAH <- function(x){
  y = filter(x, txtConstituent %in%  c("2-Methylnaphthalene", 
                                       "Acenaphthene", 
                                       "Acenaphthylene",
                                       "Anthracene", 
                                       "Fluorene", 
                                       "Naphthalene", 
                                       "Phenanthrene"), 
             !grepl("T", txtQual))
  return(y)
}

selectorHPAH <- function(x){
  y = filter(x, txtConstituent %in%  c("Fluoranthene", 
                                       "Pyrene", 
                                       "Benzo(a)anthracene", 
                                       "Chrysene", 
                                       "Benzofluoranthene", 
                                       "Benzo(a)pyrene" ,
                                       "Indeno(1,2,3-cd)pyrene", 
                                       "Dibenzo(a,h)anthracene", 
                                       "Benzo(g,h,i)perylene"), 
             !grepl("T", txtQual))
  return(y)
}

selectorPAH <- function(x){
  y = filter(x, txtConstituent %in%  c("2-Methylnaphthalene", 
                                       "Acenaphthene", 
                                       "Acenaphthylene",
                                       "Anthracene", 
                                       "Fluorene", 
                                       "Naphthalene", 
                                       "Phenanthrene",
                                       "Fluoranthene", 
                                       "Pyrene", 
                                       "Benzo(a)anthracene", 
                                       "Chrysene", 
                                       "Benzofluoranthene", 
                                       "Benzo(a)pyrene" ,
                                       "Indeno(1,2,3-cd)pyrene", 
                                       "Dibenzo(a,h)anthracene", 
                                       "Benzo(g,h,i)perylene"), 
             !grepl("T", txtQual))
  return(y)
}

selectorTCDD <- function(x){
  y = filter(x, txtConstituent %in%  c("2,3,7,8-TCDD",
                                       "1,2,3,7,8-PeCDD",
                                       "2,3,4,7,8-PeCDF",
                                       "2,3,7,8-TCDF",
                                       "1,2,3,4,7,8-HxCDD",
                                       "1,2,3,4,7,8-HxCDF",
                                       "1,2,3,7,8,9-HxCDF",
                                       "1,2,3,7,8,9-HxCDD",
                                       "1,2,3,6,7,8-HxCDD",
                                       "1,2,3,6,7,8-HxCDF",
                                       "2,3,4,6,7,8-HxCDF",
                                       "1,2,3,7,8-PeCDF",
                                       "1,2,3,4,7,8,9-HpCDF",
                                       "1,2,3,4,6,7,8-HpCDD",
                                       "1,2,3,4,6,7,8-HpCDF",
                                       "OCDF",
                                       "OCDD"),
             !grepl("T", txtQual))
  return(y)
}



assignClassPAH <- function(x){
  y <- x %>% mutate(CLASS = 
                      ifelse(txtConstituent %in%  c("Fluoranthene", 
                                                    "Pyrene", 
                                                    "Benzo(a)anthracene", 
                                                    "Chrysene", 
                                                    "Benzofluoranthene", 
                                                    "Benzo(a)pyrene" ,
                                                    "Indeno(1,2,3-cd)pyrene", 
                                                    "Dibenzo(a,h)anthracene", 
                                                    "Benzo(g,h,i)perylene"), "HPAH", 
                             ifelse(txtConstituent %in% c("2-Methylnaphthalene", 
                                                          "Acenaphthene", 
                                                          "Acenaphthylene",
                                                          "Anthracene", 
                                                          "Fluorene", 
                                                          "Naphthalene", 
                                                          "Phenanthrene"), "LPAH", NA)))
  
}

assignClassDDx <- function(x){
  y <- x %>% mutate(CLASS = 
                      ifelse(txtConstituent %in% c("o,p'-DDD", "4,4'-DDD"), "Total DDD", 
                             ifelse(txtConstituent %in% c("o,p'-DDE", "4,4'-DDE"), "Total DDE", 
                                    ifelse(txtConstituent %in% c("o,p'-DDT", "4,4'-DDT"), "Total DDT", NA)
                             )
                      )
  )
  
}




assignPEFCPAH <- function(x){
  CPAH <- data.frame( txtConstituent = c("Benzo(a)anthracene",
                                         "Benzo(a)pyrene",
                                         "Benzo(b)fluoranthene",
                                         "Benzo(k)fluoranthene",
                                         "Chrysene",
                                         "Dibenzo(a,h)anthracene",
                                         "Indeno(1,2,3-cd)pyrene"),
                      PEF = c(0.1, 1, 0.1, 0.01,0.001,1,0.1),
                      stringsAsFactors=FALSE
  )
  y = left_join(x,CPAH)
  return(y)
}

assignPEFTCCD <- function(x){
  TCDD <- data.frame( txtConstituent = c("2,3,7,8-TCDD",
                                         "1,2,3,7,8-PeCDD",
                                         "2,3,4,7,8-PeCDF",
                                         "2,3,7,8-TCDF",
                                         "1,2,3,4,7,8-HxCDD",
                                         "1,2,3,4,7,8-HxCDF",
                                         "1,2,3,7,8,9-HxCDF",
                                         "1,2,3,7,8,9-HxCDD",
                                         "1,2,3,6,7,8-HxCDD",
                                         "1,2,3,6,7,8-HxCDF",
                                         "2,3,4,6,7,8-HxCDF",
                                         "1,2,3,7,8-PeCDF",
                                         "1,2,3,4,7,8,9-HpCDF",
                                         "1,2,3,4,6,7,8-HpCDD",
                                         "1,2,3,4,6,7,8-HpCDF",
                                         "OCDF",
                                         "OCDD"),
                      PEF = c(1,
                              1,
                              0.3,
                              0.1,
                              0.1,
                              0.1,
                              0.1,
                              0.1,
                              0.1,
                              0.1,
                              0.1,
                              0.03,
                              0.01,
                              0.01,
                              0.01,
                              0.0003,
                              0.0003), 
                      stringsAsFactors = F)
  y = left_join(x,TCDD)
  return(y)
  
  
}



check_for_all_txtConstituents<- function(x, txtConstituents){

  check_df <- x %>% group_by(CLASS,txtConstituent) %>% summarise(n=n())
  if(!identical(sort(check_df$txtConstituent), sort(txtConstituents))){
    stop("ERROR: Not all analytes are present, summary could be incomplete\n
         Try using group_by(CLASS,txtConstituent) %>% summarise(n=n()) to diagnose the issue")
  }
  check_df <- x %>% group_by(CLASS,txtSampleID) %>% summarise(n=n())
  if( sum(check_df$n != length(txtConstituents)) > 0) {
    print(check_df)
    stop(sprintf("ERROR: Not all sampleIDs have %i analytes, summary could be incomplete\n
                 Try using group_by(CLASS,txtSampleID) %>% summarise(n=n()) to diagnose the issue", length(txtConstituents)))
  }
  return(TRUE)
}


check_for_all_TCDD <- function(x){
  check <- check_for_all_txtConstituents(x, c("2,3,7,8-TCDD",
                                     "1,2,3,7,8-PeCDD",
                                     "2,3,4,7,8-PeCDF",
                                     "2,3,7,8-TCDF",
                                     "1,2,3,4,7,8-HxCDD",
                                     "1,2,3,4,7,8-HxCDF",
                                     "1,2,3,7,8,9-HxCDF",
                                     "1,2,3,7,8,9-HxCDD",
                                     "1,2,3,6,7,8-HxCDD",
                                     "1,2,3,6,7,8-HxCDF",
                                     "2,3,4,6,7,8-HxCDF",
                                     "1,2,3,7,8-PeCDF",
                                     "1,2,3,4,7,8,9-HpCDF",
                                     "1,2,3,4,6,7,8-HpCDD",
                                     "1,2,3,4,6,7,8-HpCDF",
                                     "OCDF",
                                     "OCDD"))
  return(check)
}

check_for_all_CPAH <- function(x) {
  check <- check_for_all_txtConstituents(x,  c("Benzo(a)anthracene",
                                               "Benzo(a)pyrene",
                                               "Benzo(b)fluoranthene",
                                               "Benzo(k)fluoranthene",
                                               "Chrysene",
                                               "Dibenzo(a,h)anthracene",
                                               "Indeno(1,2,3-cd)pyrene"))
  return(check)
}

check_for_all_LPAH <- function(x) {
  check <- check_for_all_txtConstituents(x, c("2-Methylnaphthalene", 
                                              "Acenaphthene", 
                                              "Acenaphthylene",
                                              "Anthracene", 
                                              "Fluorene", 
                                              "Naphthalene", 
                                              "Phenanthrene"))
  return(check)
}

check_for_all_HPAH <- function(x) {
  check <- check_for_all_txtConstituents(x, c("Fluoranthene", 
                                              "Pyrene", 
                                              "Benzo(a)anthracene", 
                                              "Chrysene", 
                                              "Benzofluoranthene", 
                                              "Benzo(a)pyrene" ,
                                              "Indeno(1,2,3-cd)pyrene", 
                                              "Dibenzo(a,h)anthracene", 
                                              "Benzo(g,h,i)perylene"))
  return(check)
}

check_for_all_PAH <- function(x) {
  check <- check_for_all_txtConstituents(x,c("2-Methylnaphthalene", 
                                             "Acenaphthene", 
                                             "Acenaphthylene",
                                             "Anthracene", 
                                             "Fluorene", 
                                             "Naphthalene", 
                                             "Phenanthrene",
                                             "Fluoranthene", 
                                             "Pyrene", 
                                             "Benzo(a)anthracene", 
                                             "Chrysene", 
                                             "Benzofluoranthene", 
                                             "Benzo(a)pyrene" ,
                                             "Indeno(1,2,3-cd)pyrene", 
                                             "Dibenzo(a,h)anthracene", 
                                             "Benzo(g,h,i)perylene") )
  return(check)
}


