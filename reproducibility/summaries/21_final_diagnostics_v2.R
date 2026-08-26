studies <- c("architecture", "stress", "negative_control")
outdir <- "v2/results/final_sim_summary"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

metrics <- c(
  "MajorARI",
  "MinorityDetectionF1",
  "MinorityClusterF1",
  "AnyMicroDiscovery",
  "AbsQError",
  "SpuriousExtraClusters",
  "end_to_end_seconds"
)

reference <- "TSRClust_validated"

for (study in studies) {

  x <- readRDS(
    file.path(
      "v2/results/final_sim_combined",
      paste0(study, "_replications.rds")
    )
  )

  x <- x[x$status == "OK", , drop = FALSE]

  # --------------------------------------------------------
  # CELL-LEVEL SUMMARY
  # --------------------------------------------------------
  key <- interaction(x$cell_id, x$method, drop = TRUE)

  cell_summary <- do.call(
    rbind,
    lapply(split(x, key), function(z) {

      out <- data.frame(
        study = study,
        cell_id = z$cell_id[1],
        scenario = z$scenario_code[1],
        factor_changed = z$factor_changed[1],
        method = z$method[1],
        n = z$n[1],
        p = z$p[1],
        main_k = z$main_k[1],
        separation = z$separation[1],
        main_structure = z$main_structure[1],
        main_distribution = z$main_distribution[1],
        minority_ratio = z$minority_ratio[1],
        q = z$q[1],
        minority_sd = z$minority_sd[1],
        minority_proximity = z$minority_proximity[1],
        stringsAsFactors = FALSE
      )

      for (m in metrics) {
        v <- z[[m]]
        v <- v[is.finite(v)]
        out[[m]] <- if (length(v)) mean(v) else NA_real_
      }

      out
    })
  )

  write.csv(
    cell_summary,
    file.path(outdir, paste0(study, "_cell_summary.csv")),
    row.names = FALSE
  )

  # --------------------------------------------------------
  # PAIRED DIFFERENCES AGAINST VALIDATED TSRClust
  # Positive difference = TSRClust validated is better.
  # For AnyMicroDiscovery/AbsQError/Spurious clusters/time,
  # smaller is better, so sign is reversed.
  # --------------------------------------------------------
  ref <- x[
    x$method == reference,
    c("cell_id", "replication", metrics),
    drop = FALSE
  ]

  names(ref)[-(1:2)] <- paste0(metrics, "_ref")

  comparators <- setdiff(unique(x$method), reference)
  rows <- list()

  for (method in comparators) {

    z <- x[
      x$method == method,
      c("cell_id", "replication", metrics),
      drop = FALSE
    ]

    m <- merge(
      z,
      ref,
      by = c("cell_id", "replication"),
      all = FALSE
    )

    if (!nrow(m)) next

    for (cid in sort(unique(m$cell_id))) {

      zz <- m[m$cell_id == cid, , drop = FALSE]

      out <- data.frame(
        study = study,
        cell_id = cid,
        comparator = method,
        NPairs = nrow(zz),
        stringsAsFactors = FALSE
      )

      for (metric in metrics) {

        a <- zz[[paste0(metric, "_ref")]]
        b <- zz[[metric]]

        if (metric %in% c(
          "AnyMicroDiscovery",
          "AbsQError",
          "SpuriousExtraClusters",
          "end_to_end_seconds"
        )) {
          d <- b - a
        } else {
          d <- a - b
        }

        d <- d[is.finite(d)]

        out[[paste0(metric, "_MeanAdvantage")]] <-
          if (length(d)) mean(d) else NA_real_

        out[[paste0(metric, "_MCSE")]] <-
          if (length(d) > 1)
            sd(d) / sqrt(length(d))
          else NA_real_
      }

      rows[[length(rows) + 1L]] <- out
    }
  }

  paired <- do.call(rbind, rows)

  write.csv(
    paired,
    file.path(outdir, paste0(study, "_paired_vs_TSRvalidated.csv")),
    row.names = FALSE
  )
}

cat("FINAL_V2_DIAGNOSTICS_OK\n")
