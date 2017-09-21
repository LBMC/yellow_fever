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
ERCC <- function(data, minus = FALSE) {
  return(data$getgenes[grep('ERCC\\.', data$getgenes)])
}

get_files <- function(path, regexp) {
  file_list <- base::list.files(path = path, full.names = TRUE, recursive = TRUE)
  file_list <- file_list[grepl(regexp, file_list, perl = T)]
  return(file_list)
}
