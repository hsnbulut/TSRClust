args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("Usage: Rscript scripts/01_run_worker.R manifest.csv worker_id worker_count")
manifest_path <- args[1L]
worker_id <- as.integer(args[2L])
worker_count <- as.integer(args[3L])
source("config/config.R")
source("R/utils.R")
source("R/core_tsr.R")
source("R/data_generation.R")
source("R/evaluation_metrics.R")
source("R/competitors.R")
source("v2/config_v2.R")
source("v2/R/competitors_v2.R")
if (!requireNamespace("tclust", quietly = TRUE)) stop("Package 'tclust' is required.")
if (!requireNamespace("dbscan", quietly = TRUE)) stop("Package 'dbscan' is required.")
if (!requireNamespace("otrimle", quietly = TRUE)) stop("Package 'otrimle' is required.")
if (!requireNamespace("teigen", quietly = TRUE)) stop("Package 'teigen' is required.")
manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
selected <- manifest[((manifest$task_id - 1L) %% worker_count) == worker_id, , drop = FALSE]
if (nrow(selected) == 0L) quit(save = "no", status = 0L)

base_result_row <- function(row, rep_id, seed, method, fit_run, metrics = NULL,
                            selected_k = NA_integer_, selected_alpha = NA_real_,
                            selected_c = NA_real_, eps = NA_real_, micro_p = NA_real_,
                            micro_score = NA_real_) {
  base <- data.frame(
    study = as.character(row$study), cell_id = as_int(row$cell_id),
    task_id = as_int(row$task_id), replication = rep_id, seed = seed,
    scenario_code = as.character(row$scenario_code), method = method,
    status = if (fit_run$ok) "OK" else "ERROR",
    error = fit_run$error, warning = fit_run$warning,
    elapsed_seconds = fit_run$elapsed,
    preprocess_seconds = NA_real_,
    macro_seconds = NA_real_,
    micro_seconds = NA_real_,
    fit_seconds = fit_run$elapsed,
    end_to_end_seconds = NA_real_,
    selected_g = NA_integer_,
    otrimle_code = NA_integer_,
    estimated_noise_proportion = NA_real_,
    fit_object_bytes = if (fit_run$ok) as.numeric(utils::object.size(fit_run$value)) else NA_real_,
    n = as_int(row$n), p = as_int(row$p), main_k = as_int(row$main_k),
    separation = as_num(row$separation), main_structure = as.character(row$main_structure),
    main_distribution = as.character(row$main_distribution),
    main_weights = as.character(row$main_weights), residual_total = as_num(row$residual_total),
    minority_ratio = as_num(row$minority_ratio), q = as_int(row$q),
    minority_shape = as.character(row$minority_shape), minority_sd = as_num(row$minority_sd),
    minority_proximity = as.character(row$minority_proximity), noise_type = as.character(row$noise_type),
    factor_changed = as.character(row$factor_changed), auto_macro = as.character(row$auto_macro),
    requested_k = as_int(row$k_fit), requested_alpha = as_num(row$alpha_fit),
    requested_c = as_num(row$restr_fact), nstart = as_int(row$nstart),
    null_reps = as_int(row$null_reps),
    fitted_k = selected_k, fitted_alpha = selected_alpha, fitted_c = selected_c,
    eps = eps, micro_p_value = micro_p, micro_score = micro_score,
    stringsAsFactors = FALSE
  )
  if (is.null(metrics)) {
    metric_names <- names(evaluate_partition(c(1L, 2L), c(1L, 2L), 2L))
    for (nm in metric_names) base[[nm]] <- NA_real_
    base
  } else cbind(base, metrics)
}

run_standard_methods <- function(row, rep_id, data_obj, x, data_seed) {
  cfg <- TSR_CONFIG
  k_fit <- as_int(row$k_fit, data_obj$main_k)
  alpha_fit <- as_num(row$alpha_fit, data_obj$residual_total)
  c_fit <- as_num(row$restr_fact, 20)
  nstart <- as_int(row$nstart, cfg$nstart_main)
  method_set <- as.character(row$method_set)
  out <- list()
  macro_seed <- seed_for(cfg$base_seed, row$study, row$cell_id, rep_id, 1000L)
  set.seed(macro_seed)
  macro_run <- run_safe(fit_tclust_fixed(x, k_fit, alpha_fit, c_fit, nstart, cfg))
  if (macro_run$ok) {
    macro_labels <- as.integer(macro_run$value$cluster)
    met <- evaluate_partition(data_obj$truth, macro_labels, data_obj$main_k)
    out[[length(out) + 1L]] <- base_result_row(row, rep_id, data_seed, "TCLUST", macro_run, met,
                                               k_fit, alpha_fit, c_fit)
    if (method_set %in% c("core", "macro_focus", "runtime")) {
      set.seed(seed_for(cfg$base_seed, row$study, row$cell_id, rep_id, 2000L))
      uv_run <- run_safe(fit_tsr_from_macro(
        x, macro_labels, k_fit, validated = FALSE,
        eps_quantile = cfg$eps_quantile, eps_grid_size = cfg$eps_grid_size
      ))
      if (uv_run$ok) {
        met <- evaluate_partition(data_obj$truth, uv_run$value$cluster, data_obj$main_k)
        out[[length(out) + 1L]] <- base_result_row(
          row, rep_id, data_seed, "TSRClust_unvalidated", uv_run, met,
          k_fit, alpha_fit, c_fit, uv_run$value$micro$eps,
          uv_run$value$micro$p_value, uv_run$value$micro$score
        )
      } else out[[length(out) + 1L]] <- base_result_row(row, rep_id, data_seed, "TSRClust_unvalidated", uv_run)
      set.seed(seed_for(cfg$base_seed, row$study, row$cell_id, rep_id, 3000L))
      val_run <- run_safe(fit_tsr_from_macro(
        x, macro_labels, k_fit, validated = TRUE,
        null_reps = as_int(row$null_reps, cfg$null_reps_main),
        structure_alpha = cfg$structure_alpha,
        min_density_ratio = cfg$min_density_ratio,
        novelty_prob = cfg$novelty_prob,
        eps_quantile = cfg$eps_quantile,
        eps_grid_size = cfg$eps_grid_size
      ))
      if (val_run$ok) {
        met <- evaluate_partition(data_obj$truth, val_run$value$cluster, data_obj$main_k)
        out[[length(out) + 1L]] <- base_result_row(
          row, rep_id, data_seed, "TSRClust_validated", val_run, met,
          k_fit, alpha_fit, c_fit, val_run$value$micro$eps,
          val_run$value$micro$p_value, val_run$value$micro$score
        )
      } else out[[length(out) + 1L]] <- base_result_row(row, rep_id, data_seed, "TSRClust_validated", val_run)
    }
  } else {
    out[[length(out) + 1L]] <- base_result_row(row, rep_id, data_seed, "TCLUST", macro_run)
    if (method_set %in% c("core", "macro_focus", "runtime")) {
      out[[length(out) + 1L]] <- base_result_row(row, rep_id, data_seed, "TSRClust_unvalidated", macro_run)
      out[[length(out) + 1L]] <- base_result_row(row, rep_id, data_seed, "TSRClust_validated", macro_run)
    }
  }
  if (method_set %in% c("core", "runtime")) {
    set.seed(seed_for(cfg$base_seed, row$study, row$cell_id, rep_id, 4000L))
    tk_run <- run_safe(fit_tkmeans_fixed(x, k_fit, alpha_fit, nstart, cfg))
    if (tk_run$ok) {
      pred <- normalise_tkmeans_labels(tk_run$value$cluster, k_fit)
      out[[length(out) + 1L]] <- base_result_row(row, rep_id, data_seed, "Trimmed_kmeans", tk_run,
                                                 evaluate_partition(data_obj$truth, pred, data_obj$main_k),
                                                 k_fit, alpha_fit, NA_real_)
    } else out[[length(out) + 1L]] <- base_result_row(row, rep_id, data_seed, "Trimmed_kmeans", tk_run)
    db_run <- run_safe(fit_full_dbscan(x, cfg$eps_quantile))
    if (db_run$ok) {
      out[[length(out) + 1L]] <- base_result_row(row, rep_id, data_seed, "DBSCAN_full", db_run,
                                                 evaluate_partition(data_obj$truth, db_run$value$cluster, data_obj$main_k),
                                                 NA_integer_, NA_real_, NA_real_, db_run$value$eps)
    } else out[[length(out) + 1L]] <- base_result_row(row, rep_id, data_seed, "DBSCAN_full", db_run)
    hd_run <- run_safe(fit_full_hdbscan(x))
    if (hd_run$ok) {
      out[[length(out) + 1L]] <- base_result_row(row, rep_id, data_seed, "HDBSCAN_full", hd_run,
                                                 evaluate_partition(data_obj$truth, hd_run$value$cluster, data_obj$main_k))
    } else out[[length(out) + 1L]] <- base_result_row(row, rep_id, data_seed, "HDBSCAN_full", hd_run)
    if (method_set == "core" && data_obj$q > 0L) {
      set.seed(seed_for(cfg$base_seed, row$study, row$cell_id, rep_id, 5000L))
      oracle_run <- run_safe(fit_tclust_fixed(
        x, data_obj$main_k + data_obj$q, data_obj$noise_ratio, c_fit, nstart, cfg
      ))
      if (oracle_run$ok) {
        out[[length(out) + 1L]] <- base_result_row(
          row, rep_id, data_seed, "Oracle_TCLUST", oracle_run,
          evaluate_partition(data_obj$truth, oracle_run$value$cluster, data_obj$main_k),
          data_obj$main_k + data_obj$q, data_obj$noise_ratio, c_fit
        )
      } else out[[length(out) + 1L]] <- base_result_row(row, rep_id, data_seed, "Oracle_TCLUST", oracle_run)
    }
  }
  do.call(rbind, out)
}

run_nstart_pilot <- function(row, rep_id, data_obj, x, data_seed) {
  cfg <- TSR_CONFIG
  values <- split_ints(row$nstart_values)
  fits <- vector("list", length(values))
  for (i in seq_along(values)) {
    set.seed(seed_for(cfg$base_seed, row$study, row$cell_id, rep_id, 6000L))
    fits[[i]] <- run_safe(fit_tclust_fixed(x, as_int(row$k_fit), as_num(row$alpha_fit),
                                          as_num(row$restr_fact), values[i], cfg))
  }
  ref_idx <- which(values == max(values))[1L]
  ref <- fits[[ref_idx]]
  out <- vector("list", length(values))
  for (i in seq_along(values)) {
    ari_ref <- NA_real_
    rel_gap <- NA_real_
    obj <- NA_real_
    if (fits[[i]]$ok) obj <- fits[[i]]$value$obj
    if (fits[[i]]$ok && ref$ok) {
      ari_ref <- adjusted_rand_index(ref$value$cluster, fits[[i]]$value$cluster)
      rel_gap <- max(0, (ref$value$obj - fits[[i]]$value$obj) / max(abs(ref$value$obj), 1))
    }
    out[[i]] <- data.frame(
      study = as.character(row$study), cell_id = as_int(row$cell_id), task_id = as_int(row$task_id),
      replication = rep_id, seed = data_seed, scenario_code = as.character(row$scenario_code),
      nstart = values[i], status = if (fits[[i]]$ok) "OK" else "ERROR",
      error = fits[[i]]$error, warning = fits[[i]]$warning,
      elapsed_seconds = fits[[i]]$elapsed, objective = obj,
      ARI_vs_reference = ari_ref, RelativeObjectiveGap = rel_gap,
      n = as_int(row$n), p = as_int(row$p), stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

run_null_pilot <- function(row, rep_id, data_obj, x, data_seed) {
  cfg <- TSR_CONFIG
  values <- split_ints(row$null_reps_values)
  set.seed(seed_for(cfg$base_seed, row$study, row$cell_id, rep_id, 7000L))
  macro <- run_safe(fit_tclust_fixed(x, as_int(row$k_fit), as_num(row$alpha_fit),
                                     as_num(row$restr_fact), as_int(row$nstart, 200L), cfg))
  out <- list()
  if (!macro$ok) {
    for (v in values) out[[length(out) + 1L]] <- data.frame(
      study = as.character(row$study), cell_id = as_int(row$cell_id), task_id = as_int(row$task_id),
      replication = rep_id, seed = data_seed, scenario_code = as.character(row$scenario_code),
      null_reps = v, status = "ERROR", error = macro$error, warning = macro$warning,
      elapsed_seconds = macro$elapsed, stringsAsFactors = FALSE
    )
    return(do.call(rbind, out))
  }
  for (v in values) {
    set.seed(seed_for(cfg$base_seed, row$study, row$cell_id, rep_id, 8000L))
    fit <- run_safe(fit_tsr_from_macro(
      x, macro$value$cluster, as_int(row$k_fit), validated = TRUE, null_reps = v,
      structure_alpha = cfg$structure_alpha, min_density_ratio = cfg$min_density_ratio,
      novelty_prob = cfg$novelty_prob, eps_quantile = cfg$eps_quantile,
      eps_grid_size = cfg$eps_grid_size
    ))
    if (fit$ok) {
      met <- evaluate_partition(data_obj$truth, fit$value$cluster, data_obj$main_k)
      out[[length(out) + 1L]] <- cbind(data.frame(
        study = as.character(row$study), cell_id = as_int(row$cell_id), task_id = as_int(row$task_id),
        replication = rep_id, seed = data_seed, scenario_code = as.character(row$scenario_code),
        null_reps = v, status = "OK", error = "", warning = fit$warning,
        elapsed_seconds = fit$elapsed, micro_p_value = fit$value$micro$p_value,
        stringsAsFactors = FALSE
      ), met)
    } else {
      out[[length(out) + 1L]] <- data.frame(
        study = as.character(row$study), cell_id = as_int(row$cell_id), task_id = as_int(row$task_id),
        replication = rep_id, seed = data_seed, scenario_code = as.character(row$scenario_code),
        null_reps = v, status = "ERROR", error = fit$error, warning = fit$warning,
        elapsed_seconds = fit$elapsed, micro_p_value = NA_real_, stringsAsFactors = FALSE
      )
    }
  }
  all_names <- unique(unlist(lapply(out, names)))
  out <- lapply(out, function(z) { for (nm in setdiff(all_names, names(z))) z[[nm]] <- NA; z[all_names] })
  do.call(rbind, out)
}

run_automacro <- function(row, rep_id, data_obj, x, data_seed) {
  cfg <- TSR_CONFIG
  out <- list()
  set.seed(seed_for(cfg$base_seed, row$study, row$cell_id, rep_id, 9000L))
  auto <- run_safe(fit_auto_tclust(x, cfg, as_int(row$nstart, cfg$nstart_main)))
  if (auto$ok) {
    selected <- auto$value$selected
    labels <- as.integer(auto$value$model$cluster)
    out[[1L]] <- base_result_row(row, rep_id, data_seed, "TCLUST_auto", auto,
                                 evaluate_partition(data_obj$truth, labels, data_obj$main_k),
                                 selected$k, selected$alpha, selected$restr_fact)
    set.seed(seed_for(cfg$base_seed, row$study, row$cell_id, rep_id, 10000L))
    tsr <- run_safe(fit_tsr_from_macro(
      x, labels, selected$k, validated = TRUE,
      null_reps = as_int(row$null_reps, cfg$null_reps_main),
      structure_alpha = cfg$structure_alpha, min_density_ratio = cfg$min_density_ratio,
      novelty_prob = cfg$novelty_prob, eps_quantile = cfg$eps_quantile,
      eps_grid_size = cfg$eps_grid_size
    ))
    if (tsr$ok) out[[2L]] <- base_result_row(
      row, rep_id, data_seed, "TSRClust_auto", tsr,
      evaluate_partition(data_obj$truth, tsr$value$cluster, data_obj$main_k),
      selected$k, selected$alpha, selected$restr_fact,
      tsr$value$micro$eps, tsr$value$micro$p_value, tsr$value$micro$score
    ) else out[[2L]] <- base_result_row(row, rep_id, data_seed, "TSRClust_auto", tsr)
  } else {
    out[[1L]] <- base_result_row(row, rep_id, data_seed, "TCLUST_auto", auto)
    out[[2L]] <- base_result_row(row, rep_id, data_seed, "TSRClust_auto", auto)
  }
  do.call(rbind, out)
}

for (rr in seq_len(nrow(selected))) {
  row <- selected[rr, , drop = FALSE]
  study_dir <- file.path("v2", "results", "final_sim_raw", as.character(row$study))
  safe_dir_create(study_dir)
  outfile <- file.path(study_dir, sprintf("task_%06d.rds", as_int(row$task_id)))
  if (file.exists(outfile)) next
  result_parts <- list()
  for (rep_id in seq.int(as_int(row$rep_start), as_int(row$rep_end))) {
    data_seed <- seed_for(TSR_CONFIG$base_seed, row$study, row$cell_id, rep_id, 0L)
    data_obj <- simulate_protocol_data(row, data_seed)
    prep_start <- proc.time()[[3L]]

    x <- robust_standardize(data_obj$x)

    prep_elapsed <- proc.time()[[3L]] - prep_start
    job_type <- as.character(row$job_type)
    if (job_type == "nstart_pilot") {
      result_parts[[length(result_parts) + 1L]] <- run_nstart_pilot(row, rep_id, data_obj, x, data_seed)
    } else if (job_type == "null_pilot") {
      result_parts[[length(result_parts) + 1L]] <- run_null_pilot(row, rep_id, data_obj, x, data_seed)
    } else if (job_type == "automacro") {
      result_parts[[length(result_parts) + 1L]] <- run_automacro(row, rep_id, data_obj, x, data_seed)
    } else {
      z <- run_standard_methods(
        row, rep_id, data_obj, x, data_seed
      )

      # Correct end-to-end timing
      macro_elapsed <- z$fit_seconds[z$method == "TCLUST"][1L]

      z$preprocess_seconds <- prep_elapsed
      z$end_to_end_seconds <- prep_elapsed + z$fit_seconds

      if (length(macro_elapsed) && is.finite(macro_elapsed)) {

        ii <- which(z$method == "TCLUST")
        if (length(ii)) {
          z$macro_seconds[ii] <- macro_elapsed
          z$micro_seconds[ii] <- 0
          z$end_to_end_seconds[ii] <- prep_elapsed + macro_elapsed
        }

        for (mm in c(
          "TSRClust_unvalidated",
          "TSRClust_validated"
        )) {
          ii <- which(z$method == mm)

          if (length(ii)) {
            z$macro_seconds[ii] <- macro_elapsed
            z$micro_seconds[ii] <- z$fit_seconds[ii]
            z$end_to_end_seconds[ii] <-
              prep_elapsed + macro_elapsed + z$fit_seconds[ii]
          }
        }
      }

      z$elapsed_seconds <- z$end_to_end_seconds

      # OTRIMLE fixed K
      set.seed(seed_for(
        TSR_CONFIG$base_seed,
        row$study,
        row$cell_id,
        rep_id,
        6000L
      ))

      oo <- run_safe(
        fit_otrimle_fixed(
          x,
          G = as_int(row$k_fit, data_obj$main_k),
          erc = TSR_V2$otrimle_erc,
          npr.max = TSR_V2$otrimle_npr_max
        )
      )

      if (oo$ok) {

        oo_row <- base_result_row(
          row,
          rep_id,
          data_seed,
          "OTRIMLE_fixedK",
          oo,
          evaluate_partition(
            data_obj$truth,
            oo$value$cluster,
            data_obj$main_k
          ),
          selected_k = as_int(
            row$k_fit,
            data_obj$main_k
          )
        )

        oo_row$selected_g <-
          as_int(row$k_fit, data_obj$main_k)

        oo_row$otrimle_code <- oo$value$code
        oo_row$estimated_noise_proportion <-
          oo$value$noise_proportion

      } else {

        oo_row <- base_result_row(
          row,
          rep_id,
          data_seed,
          "OTRIMLE_fixedK",
          oo
        )
      }

      oo_row$preprocess_seconds <- prep_elapsed
      oo_row$fit_seconds <- oo$elapsed
      oo_row$end_to_end_seconds <-
        prep_elapsed + oo$elapsed
      oo_row$elapsed_seconds <-
        oo_row$end_to_end_seconds

      z <- rbind(z, oo_row)

      # tEIGEN BIC
      set.seed(seed_for(
        TSR_CONFIG$base_seed,
        row$study,
        row$cell_id,
        rep_id,
        7000L
      ))

      gmax <- min(
        6L,
        data_obj$main_k +
          max(data_obj$q, 1L) +
          TSR_V2$teigen_g_extra
      )

      tg <- run_safe(
        fit_teigen_bic(
          x,
          Gs = seq_len(gmax)
        )
      )

      if (tg$ok) {

        tg_row <- base_result_row(
          row,
          rep_id,
          data_seed,
          "tEIGEN_BIC",
          tg,
          evaluate_partition(
            data_obj$truth,
            tg$value$cluster,
            data_obj$main_k
          )
        )

        tg_row$selected_g <- tg$value$selected_G

      } else {

        tg_row <- base_result_row(
          row,
          rep_id,
          data_seed,
          "tEIGEN_BIC",
          tg
        )
      }

      tg_row$preprocess_seconds <- prep_elapsed
      tg_row$fit_seconds <- tg$elapsed
      tg_row$end_to_end_seconds <-
        prep_elapsed + tg$elapsed
      tg_row$elapsed_seconds <-
        tg_row$end_to_end_seconds

      z <- rbind(z, tg_row)

      result_parts[[length(result_parts) + 1L]] <- z
    }
  }
  saveRDS(
    do.call(rbind, result_parts),
    outfile,
    compress = "xz"
  )
  gc(FALSE)
}
