#' check gz files
#'
#' Runs `gzip -t` on a set of files to verify there is nothing wrong with
#' them before running pump and unify.
#' `gzip` must be installed and available on the system.
#' Because this is potentially making **a lot** of `system2` calls via `purrr`,
#' it is not easy to interrupt once begun.
#'
#' @param fasta the directory to check, or a list of files directly.
#'
#' @export
#' @importFrom purrr map list_rbind
rc3_check_gz = function(fasta, regexp = "fq.gz|fastq.gz") {
  is_dir = fs::dir_exists(fasta)
  if (any(is_dir)) {
    all_files = fs::dir_ls(fasta[is_dir], regexp = regexp)
  } else {
    is_file = fs::file_exists(fasta)
    all_files = fasta[is_file]
  }

  if (length(all_files) == 0) {
    cli::cli_abort(
      message = c('No files found in {.file fasta}.')
    )
  }

  gzip_res = purrr::map(
    all_files,
    \(in_file) {
      tmp = suppressWarnings(system2(
        "gzip",
        args = c("-t", in_file),
        stderr = TRUE
      ))
      data.frame(file = fs::path_file(in_file), gzip_status = tmp[2])
    },
    .progress = TRUE
  ) |>
    purrr::list_rbind()
  gzip_res
}
