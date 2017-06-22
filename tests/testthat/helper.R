#####
# Useful scripts to save and load .RDS files for testing

#' @importFrom utils head object.size
#' @importFrom stringr str_match
save_to_test <- function(obj, file_name, tolerance = sqrt(.Machine$double.eps)){
  if(!interactive())
    stop("save_to_test called not in interactive mode. Likely an error")
  
  cat("Largest sizes:\n")
  if(is.list(obj))
    print(head(sort(unlist(lapply(obj, object.size)), decreasing = T))) else
      print(object.size(obj))
  
  out_file <- paste0(
    str_match(getwd(), ".+pre"), "/tests/testthat/previous_results/", file_name, ".RDS")
  saveRDS(obj, compress = T, out_file)
  
  cat("RDS file size is ", file.size(out_file) / 1000, "KB\n", sep = "")
  
  str_tol <- if(any_numeric(obj)) 
    paste0(", tolerance = ", signif(tolerance, 4)) else ""
  
  cat("Call 'expect_equal(", deparse(substitute(obj)), ", read_to_test(\"",  file_name, "\")",
      str_tol, ")' to test\n", sep = "")
}

any_numeric <- function(x){
  if(!is.recursive(x))
    return(is.numeric(x))
  
  for(i in x){
    out <- any_numeric(i)
    if(out)
      return(TRUE)
  }
  
  return(FALSE)
}

# Sanity check
# any_numeric(c(1, 2))
# any_numeric(c("a", "b"))
# any_numeric(
#   list(a = "a", b = list(a = "a", b = 1)))
# any_numeric(
#   list(a = "a", b = list(a = "a", b = "c")))

#' @importFrom stringr str_match
read_to_test <- function(file_name){
  path <- if(!interactive()) "./previous_results/" else
    paste0(stringr::str_match(getwd(), ".+pre"), "/tests/testthat/previous_results/")
  
  readRDS(paste0(path, file_name, ".RDS"))
}

#####
# Load libraries if interactive
if(interactive()){
  library(pre)
  library(stringr)
  library(testthat)
  library(partykit)
  
  list.rules <- environment(pre)$list.rules
  list.all_rules_wo_complements <- 
    environment(pre)$list.all_rules_wo_complements
}