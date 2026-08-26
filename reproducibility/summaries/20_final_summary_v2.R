studies <- c(
  "architecture",
  "stress",
  "negative_control"
)

outdir <- "v2/results/final_sim_summary"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else median(x)
}

summarise_group <- function(z) {
  data.frame(
    N = nrow(z),
    FailureRate = mean(z$status != "OK"),
    MajorARI = safe_mean(z$MajorARI),
    MinorityDetectionF1 = safe_mean(z$MinorityDetectionF1),
    MinorityClusterF1 = safe_mean(z$MinorityClusterF1),
    AnyMicroDiscovery = safe_mean(z$AnyMicroDiscovery),
    AbsQError = safe_mean(z$AbsQError),
    SpuriousExtraClusters = safe_mean(z$SpuriousExtraClusters),
    EndToEndMedian = safe_median(z$end_to_end_seconds),
    EndToEndMean = safe_mean(z$end_to_end_seconds),
    SelectedG = safe_mean(z$selected_g)
  )
}

for (study in studies) {

  x <- readRDS(
    file.path(
      "v2/results/final_sim_combined",
      paste0(study, "_replications.rds")
    )
  )

  key <- interaction(
    x$scenario_code,
    x$method,
    drop = TRUE
  )

  sm <- do.call(
    rbind,
    lapply(split(x, key), function(z) {

      ans <- summarise_group(z)

      cbind(
        data.frame(
          study = study,
          scenario = z$scenario_code[1],
          method = z$method[1],
          stringsAsFactors = FALSE
        ),
        ans
      )
    })
  )

  sm <- sm[order(sm$scenario, sm$method), ]

  write.csv(
    sm,
    file.path(
      outdir,
      paste0(study, "_headline.csv")
    ),
    row.names = FALSE
  )
}

# ----------------------------------------------------------
# Architecture: pooled M1-M4
# ----------------------------------------------------------

a <- readRDS(
  "v2/results/final_sim_combined/architecture_replications.rds"
)

a_pos <- a[a$scenario_code %in% c("M1","M2","M3","M4"), ]

arch_pos <- do.call(
  rbind,
  lapply(split(a_pos, a_pos$method), function(z) {

    ans <- summarise_group(z)

    cbind(
      data.frame(
        method = z$method[1],
        stringsAsFactors = FALSE
      ),
      ans
    )
  })
)

arch_pos <- arch_pos[
  order(-arch_pos$MajorARI),
]

write.csv(
  arch_pos,
  file.path(
    outdir,
    "architecture_M1_M4_pooled.csv"
  ),
  row.names = FALSE
)

# ----------------------------------------------------------
# Architecture N0
# ----------------------------------------------------------

a0 <- a[a$scenario_code == "N0", ]

arch_n0 <- do.call(
  rbind,
  lapply(split(a0, a0$method), function(z) {

    ans <- summarise_group(z)

    cbind(
      data.frame(
        method = z$method[1],
        stringsAsFactors = FALSE
      ),
      ans
    )
  })
)

write.csv(
  arch_n0,
  file.path(
    outdir,
    "architecture_N0.csv"
  ),
  row.names = FALSE
)

# ----------------------------------------------------------
# Stress pooled
# ----------------------------------------------------------

s <- readRDS(
  "v2/results/final_sim_combined/stress_replications.rds"
)

stress_pool <- do.call(
  rbind,
  lapply(split(s, s$method), function(z) {

    ans <- summarise_group(z)

    cbind(
      data.frame(
        method = z$method[1],
        stringsAsFactors = FALSE
      ),
      ans
    )
  })
)

write.csv(
  stress_pool,
  file.path(
    outdir,
    "stress_pooled.csv"
  ),
  row.names = FALSE
)

# ----------------------------------------------------------
# Expanded negative controls
# ----------------------------------------------------------

n0 <- readRDS(
  "v2/results/final_sim_combined/negative_control_replications.rds"
)

neg_pool <- do.call(
  rbind,
  lapply(split(n0, n0$method), function(z) {

    ans <- summarise_group(z)

    cbind(
      data.frame(
        method = z$method[1],
        stringsAsFactors = FALSE
      ),
      ans
    )
  })
)

neg_pool <- neg_pool[
  order(neg_pool$AnyMicroDiscovery),
]

write.csv(
  neg_pool,
  file.path(
    outdir,
    "negative_control_pooled.csv"
  ),
  row.names = FALSE
)

cat("\n===== ARCHITECTURE: M1-M4 POOLED =====\n")
print(
  arch_pos[, c(
    "method",
    "MajorARI",
    "MinorityDetectionF1",
    "MinorityClusterF1",
    "AnyMicroDiscovery",
    "AbsQError",
    "EndToEndMedian"
  )],
  row.names = FALSE
)

cat("\n===== ARCHITECTURE: N0 =====\n")
print(
  arch_n0[, c(
    "method",
    "MajorARI",
    "AnyMicroDiscovery",
    "SpuriousExtraClusters",
    "EndToEndMedian"
  )],
  row.names = FALSE
)

cat("\n===== STRESS: POOLED =====\n")
print(
  stress_pool[, c(
    "method",
    "MajorARI",
    "MinorityDetectionF1",
    "MinorityClusterF1",
    "AnyMicroDiscovery",
    "AbsQError",
    "EndToEndMedian"
  )],
  row.names = FALSE
)

cat("\n===== NEGATIVE CONTROL: POOLED =====\n")
print(
  neg_pool[, c(
    "method",
    "MajorARI",
    "AnyMicroDiscovery",
    "SpuriousExtraClusters",
    "EndToEndMedian"
  )],
  row.names = FALSE
)

cat("\nFINAL_V2_SUMMARY_OK\n")
