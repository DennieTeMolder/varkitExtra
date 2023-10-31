##' @export
compute_vcf_index <- function(file) {
  file <- path.expand(file)
  if (!file.exists(file))
    stop("Could not locate the provided file: '", file, "'")
  tibble::as_tibble(.index_vcfC(file))
}

##' @export
##' @importFrom varkit varkit_index
index_vcf_fast <- function(file, ...) {
  cli::cli_progress_step("Indexing VCF")
  index <- compute_vcf_index(file)
  if (nrow(index) < 2) {
    stop("The file does not appear to be a VCF or tabular data of any sort!")
  }

  index$linum <- cumsum(c(1L, utils::head(index$n_lines, n = -1)))
  index <- index[c(1, 3, 2)]

  cli::cli_progress_step("Validating Index")
  h_idx <- 2 # Header line row-index
  if (index$id[1] != "##") {
    if (index$n_lines[1] != 1)
      stop("The file is missing a header or column names!")
    warning("No meta lines detected, using the first line as header.")
    h_idx <- 1
  }

  if (!stringr::str_detect(tolower(index$id[h_idx]), "^#?chrom$"))
    stop("At line ", index$linum[h_idx], ": Expected a header line with CHROM as the first field!")
  if (index$n_lines[h_idx] != 1L)
    stop("At line ", index$linum[h_idx], ": Expected a single a header line, but found ", index$n_lines[h_idx], " instead!")

  # Remove meta and header lines from index (these are implicitly stored in linum)
  index <- utils::tail(index, n = -h_idx)

  # Check body
  if (any(nchar(index$id) == 0))
    warning("Some CHROM fields appear to be empty. These might be skipped when reading, invalidating the index.")
  if (any(substr(index$id, 1, 1) == "#"))
    warning("The VCF contains comments. These might be skipped when reading, invalidating the index.")
  if (anyDuplicated(index$id) != 0)
    warning("The VCF appears to be unsorted, this could cause errors later.")

  varkit_index(index = index, file = file, ...)
}
