
# Test - Show that PCB, Aroclors, CPAH, LPAH, HPAH, DDx can all be done with one regular work flow

# January 21, 2019
require(dplyr)
require(tidyr)
require(stringr)
require(magrittr)

source("pcb_tools.R")
source("summation_tools.r")

# Load data cleanly Load
path = "L:/DTNA/"
db <- "PGG_SedDB_2018_SQL.mdb"
query_string = "SELECT * FROM ExR_2bSUM_All"


# Load Data
con2 <- odbcConnectAccess(paste0(path,db))  ###Make sure dbase is closed####
data <-sqlQuery(con2, query_string)
data <- data %>% mutate_if(is.factor, as.character)
close(con2)


class_name = "PDI_CALC_TOTAL_PCB"
output_suffix = ".PDI_CALC_TOTAL_PCB_CONGENERS.tsv"
# This is the entire process to produce PCB Congener results
R <- data %>% 
  selectorPCB() %>%    #!#
  process_pcb_df() %>% #!#
  mutate(CLASS = class_name) %T>% 
  assign(x = "R", value = ., pos = 1) %>% 
  general_analytical_summary(PDI = T) %>%
  output_EDD_from_sums_and_ref(x.sum = ., x.ref = R, my_signif = 2) %>%
  write.table(paste0(db, output_suffix), sep= "|", row.names = F, quote = F)


class_name = "PDI_CALC_TOTAL_AROCLORS"
output_suffix = ".PDI_CALC_TOTAL_PCB_AROCLOR.tsv"
# This is the entire process to produce Aroclor Results
data %>% 
  selectorAroclor() %>% 
  mutate(CLASS = class_name) %T>% 
  assign(x = "R", value = ., pos = 1) %>% 
  general_analytical_summary(PDI = T) %>%
  output_EDD_from_sums_and_ref(x.sum = ., x.ref = R, my_signif = 2) %>%
  write.table(paste0(db, output_suffix), sep= "|", row.names = F, quote = F)


# This is the entire process for DDx results
class_name = "PDI_CALC_TOTAL_DDX"
output_suffix = ".PDI_CALC_TOTAL_DDX.tsv"
paste0(db, output_suffix) = paste0(db, output_suffix)
# This is the entire process to produce Aroclor Results
data %>% 
  selectorDDx() %>% 
  assignClassDDx() %T>% 
  assign(x = "R", value = ., pos = 1) %>% 
  general_analytical_summary(PDI = T) %>%
  output_EDD_from_sums_and_ref(x.sum = ., x.ref = R, my_signif = 2) %>%
  write.table(paste0(db, output_suffix), sep= "|", row.names = F, quote = F)


class_name = "PDI_CALC_TOTAL_CPAH"
output_suffix = ".PDI_CALC_TOTAL_CPAH.tsv"
paste0(db, output_suffix) = paste0(db, output_suffix)
# This is the entire proces for cPAH results
data %>% 
  selectorCPAH() %>% 
  assignPEFCPAH() %>%
  mutate(CLASS = class_name) %T>%  
  assign(x = "R", value = ., pos = 1) %>% 
  general_analytical_summary(PDI = T) %>%
  output_EDD_from_sums_and_ref(x.sum = ., x.ref = R, my_signif = 2) %>%
  write.table(paste0(db, output_suffix), sep= "|", row.names = F, quote = F)

class_name = "PDI_CALC_TOTAL_LPAH"
output_suffix = ".PDI_CALC_TOTAL_LPAH.tsv"
paste0(db, output_suffix) = paste0(db, output_suffix)
# This is the entire proces for LPAH results
data %>% 
  selectorLPAH() %>% 
  mutate(CLASS = class_name) %T>%  
  assign(x = "R", value = ., pos = 1) %>% 
  general_analytical_summary(PDI = T) %>%
  output_EDD_from_sums_and_ref(x.sum = ., x.ref = R, my_signif = 2) %>%
  write.table(paste0(db, output_suffix), sep= "|", row.names = F, quote = F)

class_name = "PDI_CALC_TOTAL_HPAH"
output_suffix = ".PDI_CALC_TOTAL_HPAH.tsv"
paste0(db, output_suffix) = paste0(db, output_suffix)
# This is the entire proces for LPAH results
data %>% 
  selectorHPAH() %>% 
  mutate(CLASS = class_name) %T>%  
  assign(x = "R", value = ., pos = 1) %>% 
  general_analytical_summary(PDI = T) %>%
  output_EDD_from_sums_and_ref(x.sum = ., x.ref = R, my_signif = 2) %>%
  write.table(paste0(db, output_suffix), sep= "|", row.names = F, quote = F)


class_name = "PDI_CALC_TOTAL_PAH"
output_suffix = ".PDI_CALC_TOTAL_PAH.tsv"
# This is the entire proces for LPAH results
data %>% 
  selectorPAH() %>% 
  mutate(CLASS = class_name) %T>%  
  assign(x = "R", value = ., pos = 1) %>% 
  general_analytical_summary(PDI = T) %>%
  output_EDD_from_sums_and_ref(x.sum = ., x.ref = R, my_signif = 2) %>%
  write.table(paste0(db, output_suffix), sep= "|", row.names = F, quote = F)



# alternative for PAHs if you want to do LPAH and HPAH in the same file
# output_suffix = ".PDI_CALC_TOTAL_LplusHPAH.tsv"
# paste0(db, output_suffix) = paste0(db, output_suffix)
# # This is the entire proces for LPAH results
# data %>% 
#   selectorPAH() %>% 
#   assignClassPAH() %T>%  
#   assign(x = "R", value = ., pos = 1) %>% 
#   general_analytical_summary(PDI = T) %>%
#   output_EDD_from_sums_and_ref(x.sum = ., x.ref = R, my_signif = 2) %>%
#   write.table(paste0(db, output_suffix), sep= "|", row.names = F, quote = F)











