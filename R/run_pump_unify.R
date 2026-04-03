#' run recount-pump and recount-unify
#'
#' Given the parameters, attempts to create and run the commands to carry out the
#' runs of recount-pump and recount-unify on a set of samples.
#'
#' @param fasta the directory with the input samples to be run
#' @param outputs where to keep outputs, assumed to be the parent of fasta
#' @param monorail where is the monorail repository
#' @param recount_pump where is the recount-pump sif image
#' @param recount_unify where is the recount-unify sif image
#' @param studyid the study id
#' @param shortid short id
#' @param ncore how many cores to use
#'
#' @importFrom rlang caller_arg caller_env
#' @importFrom withr local_dir
#' @export
#' @return unify output directory
rc3_run_pump_unify = function(
  fasta = getwd(),
  outputs = fs::path_dir(fs::path_expand(fasta)),
  monorail = "monorail-external",
  recount_pump = "recount-pump_1.1.3.sif",
  recount_unify = "recount-pump_1.1.3.sif",
  reference_path = ".",
  reference = "hg38",
  studyid = "other1",
  shortid = "test",
  ncore = 1
) {
  withr::local_dir(outputs)
  arg = rlang::caller_arg(fasta)
  sample_list = check_samples(fasta, arg = rlang::caller_arg(fasta))

  monorail_paths = check_monorail(monorail, arg = rlang::caller_arg(monorail))
  recount_pump = check_exists(
    recount_pump,
    type = "file",
    arg = rlang::caller_arg(recount_pump)
  )
  recount_unify = check_exists(
    recount_unify,
    type = "file",
    arg = rlang::caller_arg(recount_unify)
  )
  reference_path = check_exists(
    reference_path,
    type = "dir",
    arg = rlang::caller_arg(reference_path)
  )

  reference_main = check_exists(
    fs::path(reference_path, reference),
    type = "dir",
    arg = arg
  )
  reference_unify = check_exists(
    fs::path(reference_path, paste0(reference, "_unify")),
    type = "dir",
    arg = arg
  )

  # pump steps
  pump_dir = fs::path(outputs, "pump_output")
  if (fs::dir_exists(pump_dir)) {
    cli::cli_inform("Deleting the directory {.file {pump_dir}}.")
    fs::dir_delete(pump_dir)
  }
  fs::dir_create(pump_dir)
  setwd(pump_dir)

  for (isample in names(sample_list)) {
    run_sample = glue::glue(
      "{monorail[1]} {recount_pump} {isample} local {reference} {ncore} {reference_path} {sample_list[[isample]][1]} {sample_list[[isample]][2]} {studyid}"
    )
    system2("/bin/bash", args = run_sample)
  }

  check_pump_outputs(unique_samples, outputs)

  # unify steps
  unify_dir = fs::path(outputs, "unify_output")
  if (fs::dir_exists(unify_dir)) {
    cli::cli_inform("Deleting the directory {.file {unify_dir}}.")
    fs::dir_delete(unify_dir)
  }
  fs::dir_create(unify_dir)
  setwd(unify_dir)
  sample_table = tibble::tibble(
    study_id = studyid,
    sample_id = fs::path_file(unique_samples)
  )
  sample_metadata_path = fs::path(outputs, "sample_metadata.tsv")
  write.table(
    sample_table,
    file = sample_metadata_path,
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t",
    quote = FALSE
  )
  run_unify = glue::glue(
    "{monorail[2]} {recount_unify} {reference} {reference_path} {unify_dir} {pump_dir} {sample_metadata_path} {ncore} {shortid}:101"
  )
  system2("/bin/bash", args = run_unify)

  return(unify_dir)
}

check_pump_outputs = function(unique_samples, outputs) {
  just_samples = fs::path_file(unique_samples)
  pump_outputs = fs::dir_ls(fs::path(outputs, "output"))
  sample_not_pump = purrr::map_lgl(just_samples, \(x) {
    !grepl(x, pump_outputs)
  })

  if (any(sample_not_pump)) {
    potential_samples = just_samples[sample_not_pump]
    names(potential_samples) = rep("*", length(potential_samples))
    cli::cli_abort(
      message = c(
        'The following samples did not get through the pump stage:',
        potential_samples,
        'i' = 'Maybe check them using {.fun rc3_check_gzip} before rerunning.'
      )
    )
  }
}

check_monorail = function(
  monorail,
  arg = rlang::caller_arg(monorail),
  call = rlang::caller_env()
) {
  mono_dir_exists = fs::dir_exists(fs::path_dir(monorail[1]))
  if (!mono_dir_exists) {
    cli::cli_abort(
      message = c(
        'Directory given for {.var monorail}, {.arg {arg}} does not exist.'
      ),
      call = call
    )
  }

  if (length(monorail) == 1) {
    out_monorail = c(
      fs::path(monorail, "singularity", "run_recount_pump.sh"),
      fs::path(monorail, "singularity", "run_recount_unify.sh")
    )
  } else if (length(monorail) == 2) {
    out_monorail = monorail
  } else {
    cli::cli_abort(
      message = c(
        '{.arg {arg}} {.code {length({arg})} returns {length(monorail)}.'
      )
    )
  }
  out_exists = fs::file_exists(out_monorail)
  if (any(!out_exists)) {
    cli::cli_abort(
      message = c(
        'One or more files that should be present in {.var monorail}, {.arg {arg}} dont exist.',
        'i' = 'Check for {.file singularity/run_recount_pump.sh} and {.file singularity/run_recount_unify.sh}'
      ),
      call = call
    )
  }

  return(out_monorail)
}

check_exists = function(
  in_path,
  type = "file",
  arg = rlang::caller_arg(in_path),
  call = rlang::caller_env()
) {
  if (type %in% "file") {
    exist_path = fs::file_exists(in_path)
  } else if (type %in% "dir") {
    exist_path = fs::dir_exists(in_path)
  }

  if (!exist_path) {
    cli::cli_abort(message = c('{.file {arg}} doesnt exist.'))
  }
  in_path
}


check_samples = function(
  fasta,
  arg = rlang::caller_arg(fasta),
  call = rlang::caller_env()
) {
  all_fasta = fs::dir_ls(fasta, regexp = "fq.gz$|fastq.gz$")
  n_files = length(all_fasta)
  if (n_files == 0) {
    cli::cli_abort(
      message = c('No fasta files found in {.arg {arg}}.'),
      call = call
    )
  }
  if ((n_files %% 2) != 0 && (n_files > 0)) {
    cli::cli_abort(
      message = c(
        '{.arg {arg}} must be a directory with pairs of fq.gz or fastq.gz files.'
      ),
      call = call
    )
  }
  unique_fasta = gsub(
    "_1.fq.gz|_2.fq.gz|_1.fastq.gz|_2.fastq.gz",
    "",
    all_fasta
  ) |>
    unique()
  n_unique = length(unique_fasta)
  if (n_files / n_unique != 2) {
    cli::cli_abort(
      message = c(
        'The directory {.arg {arg}} has an unexpected number of fq.gz files.',
        'i' = 'Do all sample files end with "_1.fq.gz" and "_2.fq.gz"?'
      ),
      call = call
    )
  }

  last_two = stringr::str_sub(unique_fasta, -2)
  has_punct = grepl("[[:punct:]]| ", last_two)

  if (any(has_punct)) {
    cli::cli_abort(
      message = c(
        'One or more file names contain punctuation or spaces in the sample ID portion.',
        'i' = 'Check all filenames for punctuation or spaces in {.file {arg}}.'
      )
    )
  }

  fasta_list = purrr::map(fs::path_file(unique_fasta), \(in_fasta) {
    tmp = grep(in_fasta, all_fasta, value = TRUE) |> sort()
    names(tmp) = NULL
    if (length(tmp) == 0) {
      return(NULL)
    }
    tmp
  })
  names(fasta_list) = fs::path_file(unique_fasta)
  notnull_list = !purrr::map_lgl(fasta_list, is.null)
  fasta_list = fasta_list[notnull_list]

  return(fasta_list)
}
