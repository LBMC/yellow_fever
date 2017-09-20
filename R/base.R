#' Select ERCC from gene list
#'
#' @param data scdata object
#' @param minus (default: false) if true return non ERCC
#' @return return a vector of ERCC gene names
#' @examples
#' \dontrun{
#' genes = ERCC(data, minus = true)
#' data = data$select(genes = genes)
#' }
#' @export ERCC
ERCC <- function(data, minus = false) {
  return(data$getgenes[grep('ERCC\\.', data$getgenes)])
}

get_files <- function(counts, regexp) {
  file_list <- list.files(path = counts, full.names = TRUE, recursive = TRUE)
  file_list <- file_list[grepl(regexp, file_list, perl = T)]
  return(file_list)
}
