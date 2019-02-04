# Forensic Tools.r

#' consider_nd_constituent
#' 
#' A preview function, to show the ratio of ND results of each constituent.
#' Caution should be used when running PCA, when variables are ND for most samples
#'
#' @param df data.frame with PGG standard variable names (eg. txtConstituent, txtQual)
#'
#' @return
#' @export
#'
#' @examples
consider_nd_constituent<- function(df){
  df %>% mutate(nd = ifelse(grepl(pattern ="U", x= txtQual), 1, 0)) %>%
    group_by(txtConstituent) %>%
    summarise(n = n(), nds = sum(nd), s = nds/n )
}


#' consider_varience_in_nd 
#'
#
#' @param df 
#'
#' @return
consider_varience_in_nd  <- function(df){
  df %>% mutate(nd = ifelse(grepl(pattern ="U", x= txtQual), 1, 0)) %>%
    group_by(txtConstituent, nd) %>%
    summarise(nds = sum(nd), min(dblLimit, na.rm = T), 
              max(dblLimit, na.rm =T),
              median(dblLimit, na.rm = T),
              var(dblLimit, na.rm = T)) %>%
    filter(nd ==1)
}


#' change_non_detect_to_zero
#' 
#' Replace NA with zero if sample was below detection limit
#'
#' @param df data.frame with PGG standard variable names (eg. txtConstituent, txtQual)
#'
#' @return data.frame with variable name dblResultPCA used in downstream functions
change_non_detect_to_zero <- function(df){
  df %>% mutate(dblResultPCA = ifelse(grepl(pattern ="U", x= txtQual),
                                      0, dblResult))
}

#' check_counts 
#'
#' @param df data.frame with PGG standard variable names (eg. txtSampleID, txtConstituent, txtQual)
#'
#' @return returns true if all samples have the same number of results, stops otherwise
check_counts <- function(df){
  df %>% group_by(txtSampleID) %>% 
    summarise(n =n()) -> y
  if (max(y$n, na.rm = T) ==  min(y$n,na.rm =T)){
    return(TRUE)
  }else{
      stop("Some samples have more results than other, normalization may be effected")
    }
}

#' normalize_by_sample
#'
#' @param df data.frame with PGG standard variable names, plus crucially dblResultPCA
#'
#' @return data.frame with results normalized as fraction of total mass across samples
#' @export
#'
#' @examples
normalize_by_sample <- function(df){
  totals <- df %>% group_by(txtSampleID) %>%
    summarise(s = sum(dblResultPCA,na.rm =T)) %>%
    select(txtSampleID, s)
  df %>% left_join(totals) %>% 
    mutate(dblResultPCAnorm = dblResultPCA/s) -> y
  return(y)
  
}


#' matrixify
#'
#' converts data.frame with character first column and numeric additional columns into a matrix
#' 
#' @param df data.frame with analyte names in the first column and numerics only elsewhere 
#'
#' @return matrix with all but the first row, first row is used as rownames
matrixify <- function(df){
  m = as.matrix(df[,2:ncol(df)])
  rownames(m) <- df[,1]
  return(m)
}

# see page 227 of Env. Forensics for Miesch calculated CD
miesch_cd <- function(obs,low_rank_approx){
  (var(obs) - var(obs-low_rank_approx)) / var(obs)
}


#' remove_zero_variance_cols
#'
#' PCA can't run on a matrix with variable with zero variance
#'
#' @param m matrix with numeric entries only
#'
#' @return matrix with columns reoved that have zero variance
remove_zero_variance_cols <- function(m){
  ind <- (which(!apply(m, 2, var) == 0))
  return(m[,ind])
}

#' remove_zero_variance_row 
#'
#' @param m matrix with numeric entries 
#'
#' @return matrix with row removed that have zero variance
remove_zero_variance_row <- function(m){
  ind <- (which(!apply(m, 1, var) == 0))
  return(m[ind,])
}


#' remove columns
#'
#' @param x matrix or data.frame
#' @param remove_ind column index used to remove colums, all other indices are kept
#'
#' @return
remove_columns <- function(x, remove_ind = c()){
  return(x[,-remove_ind])
}


#' better_biplot
#' 
#' Makes a simple PCA plot with aesthetics that are an improvement over
#' built in prcomp biplot() function 
#'
#' @param pca precomp object 
#' @param data.matrix data.matrix used as t(data.matrix) in prcomp principal compenents analysis
#' @param f factor defaulting at 10 for drawing biplot vectors
#' @param dim1 first PC score to plot on x-axis (default is 1)
#' @param dim2 second PC scores to plot on y-axis (default is 2)
#' @param label_as_number boolean deafult is true - labels are column position in data.matrix, 
#' if false, labels are colnames in data.matrix 
#'
#' @return plot command to make PCA biplot
better_biplot<- function(pca, data.matrix, f = 10, dim1 =1, dim2=2, label_as_number = T){
  
  
  var_percent <- round(100*pca$sdev^2/sum(pca$sdev^2),1)
  plot(pca$x[,dim1 ], 
       pca$x[,dim2], 
       pch = 20, 
       col = "#00000050", 
       xlab = paste0("PC",dim1, " (", var_percent[dim1], "% of Variance)"),
       ylab = paste0("PC",dim2, " (",var_percent[dim2], "% of Variance)"))
  if (label_as_number){
    text(pca$x[,dim1 ], 
         pca$x[,dim2], 
         labels= 1:ncol(data.matrix), 
         cex = .5, 
         pos =4)
  } else {
    text(pca$x[,dim1 ], 
         pca$x[,dim2], 
         labels= colnames(data.matrix), 
         cex = .5, 
         pos = 4)
  }
  arrows(0,0,f*pca$rotation[,dim1 ],f*pca$rotation[,dim2], col = "#0000FF50", length = .05)
  text(f*pca$rotation[,dim1 ],
       f*pca$rotation[,dim2], 
       label = rownames(data.matrix),
       col = "#0000FF50",
       cex = .5) 
}

# 
# # Full Work Flow
# require(dplyr)
# require(tidyr)
# require(stringr)
# require(magrittr)
# require(DBI)
# require(RODBC)
# 
# # Load data cleanly Load
# path = "L:/DTNA/"
# db <- "PGG_SedDB_2018_SQL.mdb"
# query_string = "SELECT * FROM ExR_2bSUM_All"
# source("pcb_tools.R")
# source("summation_tools.r")
# # Load Data
# con2 <- odbcConnectAccess(paste0(path,db))  ###Make sure dbase is closed####
# data <-sqlQuery(con2, query_string)
# data <- data %>% mutate_if(is.factor, as.character)
# 
# testdata <- data %>%
#   selectorTCDD() #!# filters to 17 Dioxin Furans
# 
# data.matrix <- testdata %>% change_non_detect_to_zero() %>%
#   normalize_by_sample() %>%
#   select(txtSampleID, txtConstituent, dblResultPCAnorm) %>%
#   spread(key = txtSampleID, value = dblResultPCAnorm) %>%
#   matrixify() %>%
#   remove_zero_variance_cols()
# class(data.matrix)
# pca <- prcomp(t(data.matrix), scale = T, center = T)
# better_biplot(pca, data.matrix ,f = 10)
# better_biplot(pca, data.matrix ,f = 10, label_as_number = F)
# better_biplot(pca, data.matrix, f = 10, 3,2)
# dim(remove_columns(data.matrix, c(1,2,3)))
# # variable-by-variable goodness
# 
# # recomose the matrix
# M = (t(pca$x[] %*% t(pca$rotation[])) * pca$scale  + pca$center)
# # recompose matrix using only < n > compoents
# n = 3
# M_approx = (t(pca$x[,1:n] %*% t(pca$rotation[,1:n])) * pca$scale  + pca$center)
# 
# par(mfrow = c(4,4), mar = c(2,2,2,2))
# for(i in 2:nrow(M)) {
#   cd = round(miesch_cd(M[i,], M_approx[i,]),2)
#   plot(M[i,],
#        M_approx[i,],
#        pch = 20,
#        cex = .5,
#        col = "#00000050",
#        main = paste(rownames(M)[i], "(rank:" ,n ,")") ,
#        ylab = "low-rank prediction",
#        xlab = "measured-value",
#        ylim = c(0,max(c( M[i,],
#                        M_approx[i,]) ) ),
#        xlim = c(0,max(c( M[i,],
#                          M_approx[i,]) ) ) )
#   abline(0,1, col = "gray")
#   text(M[i,],
#        M_approx[i,],
#        labels= 1:ncol(M),
#        cex = .75,
#        pos = 2,
#        col = "#00000050")
#   text(0,0.95*max(M_approx[i,], M[i,]), label = paste("CD =", cd), pos =4, cex = 1 )
# }