#!/usr/bin/env Rscript

# Call chr15 DMRs between hap1 and hap2 independently for each sample.
#
# Input BED columns (1-based column numbers):
#   1 chromosome, 2 BED start, 6 total coverage, 7 methylated reads
# The BED start is converted from 0-based to the 1-based coordinate DSS uses.

suppressPackageStartupMessages({
  library(data.table)
  library(DSS)
})

defaults <- list(
  input_dir = "/home/rare/arlen/pws-as-hifi-cis-methylation/results/06_methylation",
  output_dir = "/home/rare/arlen/pws-as-hifi-cis-methylation/results/06_methylation/DSS",
  chromosome = "chr15",
  samples = NULL,
  ncores = 1L,
  smoothing_span = 500L,
  delta = 0,
  p_threshold = 1e-5,
  minlen = 50L,
  min_cg = 3L,
  merge_distance = 100L,
  pct_sig = 0.5,
  overwrite = FALSE,
  dry_run = FALSE
)

usage <- function(status = 0L) {
  cat(
    "Usage: dss_chr15_hap1_vs_hap2.R [options]\n\n",
    "Calls DSS DMRs for hap1 versus hap2 separately in every sample.\n\n",
    "Options:\n",
    "  --input-dir DIR          Input root [results/06_methylation]\n",
    "  --output-dir DIR         Output directory [results/06_methylation/DSS]\n",
    "  --chromosome CHR         Chromosome name [chr15]\n",
    "  --samples ID1,ID2        Only process these sample IDs [all]\n",
    "  --ncores INT             Cores used by each DSS test [1]\n",
    "  --smoothing-span INT     DSS smoothing span in bp [500]\n",
    "  --delta NUM              Minimum methylation difference [0]\n",
    "  --p-threshold NUM        DSS DMR p-value threshold [1e-5]\n",
    "  --minlen INT             Minimum DMR length [50]\n",
    "  --min-cg INT             Minimum CpGs per DMR [3]\n",
    "  --merge-distance INT     Merge DMRs within this distance [100]\n",
    "  --pct-sig NUM            Required fraction significant CpGs [0.5]\n",
    "  --overwrite              Replace existing per-sample outputs\n",
    "  --dry-run                Validate inputs without running DSS\n",
    "  -h, --help               Show this help\n",
    sep = ""
  )
  quit(save = "no", status = status)
}

parse_args <- function(args) {
  opts <- defaults
  value_options <- c(
    "--input-dir", "--output-dir", "--chromosome", "--samples", "--ncores",
    "--smoothing-span", "--delta", "--p-threshold", "--minlen", "--min-cg",
    "--merge-distance", "--pct-sig"
  )
  key_map <- c(
    "--input-dir" = "input_dir", "--output-dir" = "output_dir",
    "--chromosome" = "chromosome", "--samples" = "samples",
    "--ncores" = "ncores", "--smoothing-span" = "smoothing_span",
    "--delta" = "delta", "--p-threshold" = "p_threshold",
    "--minlen" = "minlen", "--min-cg" = "min_cg",
    "--merge-distance" = "merge_distance", "--pct-sig" = "pct_sig"
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("-h", "--help")) usage(0L)
    if (arg == "--overwrite") {
      opts$overwrite <- TRUE
      i <- i + 1L
      next
    }
    if (arg == "--dry-run") {
      opts$dry_run <- TRUE
      i <- i + 1L
      next
    }
    if (!arg %in% value_options) {
      stop("Unknown option: ", arg, call. = FALSE)
    }
    if (i == length(args)) {
      stop("Missing value for ", arg, call. = FALSE)
    }
    opts[[key_map[[arg]]]] <- args[[i + 1L]]
    i <- i + 2L
  }

  integer_keys <- c("ncores", "smoothing_span", "minlen", "min_cg", "merge_distance")
  numeric_keys <- c("delta", "p_threshold", "pct_sig")
  for (key in integer_keys) opts[[key]] <- as.integer(opts[[key]])
  for (key in numeric_keys) opts[[key]] <- as.numeric(opts[[key]])
  if (!is.null(opts$samples)) {
    opts$samples <- unique(trimws(strsplit(opts$samples, ",", fixed = TRUE)[[1L]]))
    opts$samples <- opts$samples[nzchar(opts$samples)]
  }

  if (is.na(opts$ncores) || opts$ncores < 1L) stop("--ncores must be >= 1", call. = FALSE)
  if (is.na(opts$smoothing_span) || opts$smoothing_span < 1L) stop("--smoothing-span must be >= 1", call. = FALSE)
  if (is.na(opts$delta) || opts$delta < 0 || opts$delta > 1) stop("--delta must be between 0 and 1", call. = FALSE)
  if (is.na(opts$p_threshold) || opts$p_threshold <= 0 || opts$p_threshold > 1) stop("--p-threshold must be in (0, 1]", call. = FALSE)
  if (is.na(opts$pct_sig) || opts$pct_sig <= 0 || opts$pct_sig > 1) stop("--pct-sig must be in (0, 1]", call. = FALSE)
  if (!dir.exists(opts$input_dir)) stop("Input directory does not exist: ", opts$input_dir, call. = FALSE)
  opts
}

read_haplotype <- function(path, chromosome) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)

  dat <- fread(
    path,
    header = FALSE,
    select = c(1L, 2L, 6L, 7L),
    showProgress = FALSE
  )
  setnames(dat, c("chr", "start0", "N", "X"))
  dat <- dat[chr == chromosome]
  if (nrow(dat) == 0L) stop("No ", chromosome, " records in ", path, call. = FALSE)

  dat[, `:=`(
    start0 = as.numeric(start0),
    N = as.numeric(N),
    X = as.numeric(X)
  )]
  if (anyNA(dat)) stop("Missing or non-numeric chr15 count fields in ", path, call. = FALSE)
  if (any(dat$start0 < 0) || any(dat$N <= 0) || any(dat$X < 0) || any(dat$X > dat$N)) {
    stop("Invalid coordinates or methylation counts in ", path, call. = FALSE)
  }

  # Combining duplicate positions preserves their total methylated/read counts.
  dat[, pos := start0 + 1]
  dat <- dat[, .(N = sum(N), X = sum(X)), by = .(chr, pos)]
  setorder(dat, chr, pos)
  dat[, `:=`(pos = as.integer(pos), N = as.integer(N), X = as.integer(X))]
  as.data.frame(dat[, .(chr, pos, N, X)])
}

write_dmr_outputs <- function(dmrs, sample_id, output_dir, chromosome) {
  prefix <- file.path(output_dir, paste0(sample_id, ".", chromosome, ".hap1_vs_hap2"))
  result <- as.data.table(dmrs)
  result[, `:=`(sample = sample_id, comparison = "hap1_vs_hap2")]
  setcolorder(result, c("sample", "comparison", setdiff(names(result), c("sample", "comparison"))))
  fwrite(result, paste0(prefix, ".DMR.tsv"), sep = "\t", quote = FALSE, na = "NA")

  if (nrow(result) == 0L) {
    bed <- data.table(
      chr = character(), start = integer(), end = integer(), name = character(),
      score = integer(), strand = character(), diff.Methy = numeric(),
      meanMethy1 = numeric(), meanMethy2 = numeric(), nCG = integer(),
      areaStat = numeric()
    )
  } else {
    bed <- result[, .(
      chr,
      start = pmax(0L, as.integer(start) - 1L),
      end = as.integer(end),
      name = sprintf("%s_DMR_%d", sample_id, seq_len(.N)),
      score = pmin(1000L, as.integer(round(abs(diff.Methy) * 1000))),
      strand = ".",
      diff.Methy,
      meanMethy1,
      meanMethy2,
      nCG,
      areaStat
    )]
  }
  fwrite(bed, paste0(prefix, ".DMR.bed"), sep = "\t", col.names = FALSE, quote = FALSE)
  invisible(prefix)
}

opts <- tryCatch(parse_args(commandArgs(trailingOnly = TRUE)), error = function(e) {
  message("ERROR: ", conditionMessage(e))
  usage(2L)
})

sample_dirs <- list.dirs(opts$input_dir, full.names = FALSE, recursive = FALSE)
sample_dirs <- sort(sample_dirs[sample_dirs != "DSS"])
samples <- if (is.null(opts$samples)) sample_dirs else opts$samples
if (length(samples) == 0L) stop("No samples selected", call. = FALSE)

unknown_samples <- setdiff(samples, sample_dirs)
if (length(unknown_samples)) {
  stop("Sample directories not found: ", paste(unknown_samples, collapse = ", "), call. = FALSE)
}

input_pairs <- lapply(samples, function(sample_id) {
  list(
    sample = sample_id,
    hap1 = file.path(opts$input_dir, sample_id, paste0(sample_id, ".hap1.bed")),
    hap2 = file.path(opts$input_dir, sample_id, paste0(sample_id, ".hap2.bed"))
  )
})
missing_inputs <- unlist(lapply(input_pairs, function(x) c(x$hap1, x$hap2)[!file.exists(c(x$hap1, x$hap2))]))
if (length(missing_inputs)) stop("Missing input file(s):\n", paste(missing_inputs, collapse = "\n"), call. = FALSE)

dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)

parameter_table <- data.table(
  parameter = c(
    "input_dir", "output_dir", "chromosome", "ncores", "smoothing",
    "smoothing_span", "equal_disp", "delta", "p_threshold", "minlen",
    "min_cg", "merge_distance", "pct_sig", "DSS_version"
  ),
  value = c(
    opts$input_dir, opts$output_dir, opts$chromosome, opts$ncores, TRUE,
    opts$smoothing_span, FALSE, opts$delta, opts$p_threshold, opts$minlen,
    opts$min_cg, opts$merge_distance, opts$pct_sig, as.character(packageVersion("DSS"))
  )
)
fwrite(parameter_table, file.path(opts$output_dir, "run_parameters.tsv"), sep = "\t")

summary_rows <- vector("list", length(input_pairs))
had_error <- FALSE

for (i in seq_along(input_pairs)) {
  pair <- input_pairs[[i]]
  sample_id <- pair$sample
  prefix <- file.path(opts$output_dir, paste0(sample_id, ".", opts$chromosome, ".hap1_vs_hap2"))
  output_tsv <- paste0(prefix, ".DMR.tsv")
  started <- Sys.time()
  message(sprintf("[%d/%d] %s", i, length(input_pairs), sample_id))

  if (file.exists(output_tsv) && !opts$overwrite && !opts$dry_run) {
    message("  Skipping existing output (use --overwrite to replace it)")
    existing_dmr_count <- nrow(fread(output_tsv, select = 1L, showProgress = FALSE))
    summary_rows[[i]] <- data.table(
      sample = sample_id, status = "skipped_existing", hap1_chr_cpgs = NA_integer_,
      hap2_chr_cpgs = NA_integer_, tested_cpgs = NA_integer_, dmrs = existing_dmr_count,
      elapsed_seconds = 0, message = ""
    )
    fwrite(rbindlist(summary_rows[seq_len(i)], fill = TRUE), file.path(opts$output_dir, "run_summary.tsv"), sep = "\t")
    next
  }

  # local() releases the large BSseq/DML objects before the next sample starts.
  result <- tryCatch(local({
    hap1 <- read_haplotype(pair$hap1, opts$chromosome)
    hap2 <- read_haplotype(pair$hap2, opts$chromosome)
    message(sprintf("  Loaded %s hap1 CpGs and %s hap2 CpGs", nrow(hap1), nrow(hap2)))

    if (opts$dry_run) {
      list(status = "validated", n1 = nrow(hap1), n2 = nrow(hap2), tested = NA_integer_, ndmr = NA_integer_)
    } else {
      bsobj <- makeBSseqData(list(hap1, hap2), c("hap1", "hap2"))
      dml <- DMLtest(
        bsobj,
        group1 = "hap1",
        group2 = "hap2",
        equal.disp = FALSE,
        smoothing = TRUE,
        smoothing.span = opts$smoothing_span,
        ncores = opts$ncores
      )
      dmrs <- callDMR(
        dml,
        delta = opts$delta,
        p.threshold = opts$p_threshold,
        minlen = opts$minlen,
        minCG = opts$min_cg,
        dis.merge = opts$merge_distance,
        pct.sig = opts$pct_sig
      )
      write_dmr_outputs(dmrs, sample_id, opts$output_dir, opts$chromosome)
      message(sprintf("  Wrote %d DMRs", nrow(dmrs)))
      list(status = "complete", n1 = nrow(hap1), n2 = nrow(hap2), tested = sum(!is.na(dml$pval)), ndmr = nrow(dmrs))
    }
  }), error = function(e) {
    had_error <<- TRUE
    message("  ERROR: ", conditionMessage(e))
    list(status = "failed", n1 = NA_integer_, n2 = NA_integer_, tested = NA_integer_, ndmr = NA_integer_, error = conditionMessage(e))
  })

  summary_rows[[i]] <- data.table(
    sample = sample_id,
    status = result$status,
    hap1_chr_cpgs = result$n1,
    hap2_chr_cpgs = result$n2,
    tested_cpgs = result$tested,
    dmrs = result$ndmr,
    elapsed_seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 2),
    message = if (is.null(result$error)) "" else result$error
  )
  fwrite(rbindlist(summary_rows[seq_len(i)], fill = TRUE), file.path(opts$output_dir, "run_summary.tsv"), sep = "\t")
  rm(result)
  invisible(gc())
}

capture.output(sessionInfo(), file = file.path(opts$output_dir, "sessionInfo.txt"))
if (had_error) quit(save = "no", status = 1L)
message("Finished. Results: ", opts$output_dir)
