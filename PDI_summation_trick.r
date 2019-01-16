general_analytical_summary_PDI <- function(df){
  df.maxLimits <- df %>% select(CLASS,
                              txtConstituent,
                              txtSampleID,
                              dblLimit,
                              dblResult) %>%
    mutate(DL = is.na(dblResult)) %>%                      # This labels ND samples as DL = TRUE
    group_by(CLASS, txtSampleID, DL) %>%                   # This groups all the ND samples together per CLASS,txtSampleID
    summarise(dblmaxLimitSample = max(dblLimit, na.rm= T), # maxLimit calculates highest limits among NDs, 
              totalNDs = n()) %>%     
    filter(DL == TRUE) %>% ungroup()  # this filters so that we only have maxLimit from the samples with atleaast one ND
  
  df.sums <- df %>% select(CLASS,
                           txtConstituent,
                           txtSampleID, 
                           dblResult, 
                           dblLimit, 
                           txtQual) %>% 
    left_join(df.maxLimits) %>%
    # page 3. | Appendix E: 2018 PDI Surface Sediment Database Description
    # "Calculated totals are the sum of all detected results and one half of the highest reporting
    # detection limit for the results reported as not detected (ND) in a sample."
    mutate(dblResultzeroND = ifelse(is.na(dblResult), 0*dblmaxLimitSample,   dblResult)) %>%
    mutate(dblResulthalfND = ifelse(is.na(dblResult), 0.5*dblmaxLimitSample, dblResult)) %>%
    mutate(dblResultfullND = ifelse(is.na(dblResult), 1*dblmaxLimitSample,   dblResult)) %>% 
    group_by(CLASS, txtSampleID) %>%
    summarise(sumResult_0ND    = sum(dblResultzeroND, na.rm=TRUE), 
              sumResult_halfND = sum(dblResulthalfND, na.rm=TRUE),
              sumResult_fullND = sum(dblResultfullND, na.rm=TRUE),
              maxLimit         = max(dblLimit), 
              n                = n(),
              nd               = sum(is.na(dblResult)),
              txtQual       = return_all_unique_flags(txtQual) ) %>%
    mutate(`ND_zero` = ifelse(sumResult_0ND == 0, maxLimit, sumResult_0ND), 
           `ND_half` = ifelse(sumResult_0ND == 0, maxLimit, sumResult_halfND),
           `ND_full` = ifelse(sumResult_0ND == 0, maxLimit, sumResult_fullND))
  return(df.sums)
  
}

# HERE IS A TEST:
A <- data %>% filter(txtAnalysisMethod == "SW8082A", !grepl(txtQual, pattern = "T")) 
A$CLASS <- "AROCLORS"
Atest <- A
Atest$dblLimit[1] <- 10000#head(Atest)
#J5-SC-40to60-102218




# HERE IS A TEST WHERE I CHANGED THE MAX DETECTION LIMIT OF ONLY ONE ANALYTE IN 
general_analytical_summary_PDI(Atest) %>% filter(txtSampleID== "J5-SC-40to60-102218")
# vs. old method
general_analytical_summary(Atest) %>% filter(txtSampleID== "J5-SC-40to60-102218")
general_analytical_summary(A) %>% filter(txtSampleID== "J5-SC-40to60-102218")
general_analytical_summary_PDI(A) %>% filter(txtSampleID== "J5-SC-40to60-102218")




