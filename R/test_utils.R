#' @importFrom utils head object.size
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
  
  cat("Call 'expect_equal(", deparse(substitute(obj)), ", read_to_test(\"",  file_name, "\"), tolerance = ",
      tolerance, ")' to test\n", sep = "")
}

#' @importFrom stringr str_match
read_to_test <- function(file_name){
  path <- if(!interactive()) "./previous_results/" else
    paste0(str_match(getwd(), ".+pre"), "/tests/testthat/previous_results/")
  
  readRDS(paste0(path, file_name, ".RDS"))
}