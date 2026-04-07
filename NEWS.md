# rc3helpers 0.0.4

* Added better checking that a run properly completed by checking contents of **std.out**.
* Made it possible to run just the **pump**, **unify**, or **both** pump & unify steps (see argument `run`  of `rc3_run_pump_unify`).
* Can keep previous **pump** outputs, and only re-run the failed samples (see argument `delete_oldpump` of `rc3_run_pump_unify`).
* Can generate and **show** the list of commands instead of running them (see argument `run_or_show` of `rc3_run_pump_unify`).

# rc3helpers 0.0.3

* Suppresses output that is going to *stdout*.
* Handles both fq.gz and fastq.gz files.
