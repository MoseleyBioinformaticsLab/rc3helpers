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
#' @param run which parts of the workflow to run? (both, pump, unify)
#' @param delete_oldpump should the old recount-pump results be deleted and regenerated?
#' @param run_or_show do you want the commands shown, or actually run?
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
  ncore = 1,
  run = "both",
  delete_oldpump = "yes",
  run_or_show = "run"
) {
  withr::local_dir(outputs)

  sample_list = check_samples(
    fasta,
    fasta_arg = rlang::caller_arg(fasta)
  )

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
  if (run %in% c("both", "pump")) {
    pump_dir = fs::path(outputs, "pump_output")

    # if delete and rerun, nuke the old results.
    # otherwise, check the outputs, delete only what is necessary
    # and rerun them before running unify
    if (delete_oldpump %in% "yes") {
      if (fs::dir_exists(pump_dir)) {
        cli::cli_inform("Deleting the directory {.file {pump_dir}}.")
        fs::dir_delete(pump_dir)
      }
      fs::dir_create(pump_dir)
      run_list = sample_list
    } else {
      notdone_samples = check_pump_outputs(
        names(sample_list),
        pump_dir,
        die = "no"
      )
      run_list = sample_list[notdone_samples]
    }

    setwd(pump_dir)

    cli::cli_inform("Running recount-pump on samples.")
    for (isample in names(run_list)) {
      run_sample = glue::glue(
        "{monorail_paths[1]} {recount_pump} {isample} local {reference} {ncore} {reference_path} {sample_list[[isample]][1]} {sample_list[[isample]][2]} {studyid}"
      )
      if (run_or_show %in% "run") {
        system2("/bin/bash", args = run_sample, stdout = TRUE, stderr = "")
      } else {
        cli::cli_inform(
          message = c(
            "To run from the command line:",
            'i' = paste0("{.field /bin/bash ", run_sample, "}")
          )
        )
      }
    }
  }

  # unify steps
  if (run %in% c("both", "unify")) {
    check_pump_outputs(
      names(sample_list),
      pump_dir,
      die = "yes"
    )
    unify_dir = fs::path(outputs, "unify_output")
    if (fs::dir_exists(unify_dir)) {
      cli::cli_inform("Deleting the directory {.file {unify_dir}}.")
      fs::dir_delete(unify_dir)
    }
    fs::dir_create(unify_dir)
    setwd(unify_dir)
    sample_table = data.frame(
      study_id = studyid,
      sample_id = names(sample_list)
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
    cli::cli_inform("Running recount-unify.")
    run_unify = glue::glue(
      "{monorail_paths[2]} {recount_unify} {reference} {reference_path} {unify_dir} {pump_dir} {sample_metadata_path} {ncore} {shortid}:101"
    )
    if (run_or_show %in% "run") {
      system2("/bin/bash", args = run_unify, stdout = TRUE, stderr = "")
    } else {
      cli::cli_inform(
        message = c(
          "To run from the command line: ",
          'i' = paste0("{.field /bin/bash ", run_unify, "}")
        )
      )
    }

    cli::cli_inform("Done!")
  }

  return(unify_dir)
}

check_pump_outputs = function(unique_samples, pump_dir, die = "yes") {
  # logic:
  #   get the sample directories for pump;
  #   generate path to std.out (which it will have it ran at all)
  #   check for the done message in std.out
  #   filter the std.out paths for those that did complete
  #   check sample ids that they have directories with completion
  #   if any don't, either die with the error, or return the missing ones
  pump_outputs = fs::dir_ls(fs::path(pump_dir, "output"))
  std_out_locs = fs::path(pump_outputs, "std.out")

  job_complete = purrr::map_lgl(std_out_locs, \(in_file) {
    snakemake_log = readLines(in_file)
    has_done = which(grepl("19 of 19 steps \\(100\\%\\) done", snakemake_log))
    if ((length(has_done) > 0) && (has_done == (length(snakemake_log) - 1))) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  })

  std_out_complete = std_out_locs[job_complete]

  sample_not_complete = purrr::map_lgl(unique_samples, \(x) {
    all(!grepl(x, std_out_complete))
  })

  if (any(sample_not_complete)) {
    potential_samples = unique_samples[sample_not_complete]
    if (die %in% "yes") {
      names(potential_samples) = rep("*", length(potential_samples))
      cli::cli_abort(
        message = c(
          'The following samples did not get through the pump stage:',
          potential_samples,
          'i' = 'Maybe check them using {.fun rc3_check_gzip} before rerunning.'
        )
      )
    } else {
      return(potential_samples)
    }
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
  fasta_arg = rlang::caller_arg(fasta),
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
        '{.arg {fasta_arg}} must be a directory with pairs of fq.gz or fastq.gz files.'
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
