options(stringsAsFactors = FALSE)

root <- normalizePath(file.path(getwd(), ".."), mustWork = TRUE)
out <- normalizePath(getwd(), mustWork = TRUE)

dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "supplementary"), recursive = TRUE, showWarnings = FALSE)

read_csv <- function(...) read.csv(file.path(root, ...), check.names = FALSE)
fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "--", formatC(x, format = "f", digits = digits))
}
latex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([_&%$#{}])", "\\\\\\1", x, perl = TRUE)
  x
}
write_table <- function(path, tabular, caption, label, notes = NULL) {
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\caption{" %+% caption %+% "}",
    "\\label{" %+% label %+% "}",
    "\\centering",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{2pt}",
    tabular
  )
  if (!is.null(notes)) lines <- c(lines, "\\begin{flushleft}\\footnotesize " %+% notes %+% "\\end{flushleft}")
  lines <- c(lines, "\\end{table}")
  writeLines(lines, file.path(out, path))
}
fit_width <- function(tabular) {
  z <- sub("\\\\begin\\{tabular\\}", "\\\\resizebox{\\\\textwidth}{!}{%\n\\\\begin{tabular}", tabular)
  sub("\\\\end\\{tabular\\}", "\\\\end{tabular}%\n}", z)
}
`%+%` <- function(a, b) paste0(a, b)

manifest_summary <- do.call(rbind, lapply(c("architecture", "stress", "negative_control"), function(study) {
  m <- read_csv("github_results", "simulation_design", paste0(study, ".csv"))
  r <- read.csv(gzfile(file.path(root, "github_results", "raw_results", paste0(study, "_replications.csv.gz"))), check.names = FALSE)
  data.frame(
    study = study,
    cells = length(unique(m$cell_id)),
    generated_datasets = sum(m$rep_end - m$rep_start + 1),
    method_evaluations = nrow(r),
    successful_evaluations = sum(r$status == "OK"),
    methods = length(unique(r$method))
  )
}))
write.csv(manifest_summary, file.path(out, "tables", "simulation_scale.csv"), row.names = FALSE)

design_tab <- data.frame(
  Component = c("Architecture", "Stress", "Negative control"),
  Purpose = c(
    "Canonical architecture grid",
    "One-factor stress departures",
    "No-minority residual controls"
  ),
  Cells = manifest_summary$cells,
  Datasets = manifest_summary$generated_datasets,
  Evaluations = manifest_summary$method_evaluations,
  MainFactors = c(
    "$n$, $p$, separation, covariance, N0/M1--M4",
    "$n$, $p$, $K$, separation, tails, imbalance, minority geometry, noise",
    "$p$, residual fraction, tails, diffuse-noise type, $q=0$"
  )
)
tab <- paste(c(
  "\\begin{tabular}{p{0.18\\linewidth}p{0.28\\linewidth}rrrp{0.25\\linewidth}}",
  "\\hline",
  "Study & Purpose & Cells & Datasets & Evaluations & Main factors \\\\",
  "\\hline",
  apply(design_tab, 1, function(z) paste(latex_escape(z[1]), latex_escape(z[2]), z[3], z[4], z[5], z[6], sep = " & ") %+% " \\\\"),
  "\\hline",
  "\\end{tabular}"
), collapse = "\n")
write_table("tables/table1_simulation_design.tex", tab, "Final simulation design and scale.", "tab:simulation-design")

arch <- read_csv("github_results", "summary_results", "architecture_M1_M4_pooled.csv")
stress <- read_csv("github_results", "summary_results", "stress_pooled.csv")
neg <- read_csv("github_results", "summary_results", "negative_control_pooled.csv")
methods <- c("Oracle_TCLUST", "TSRClust_validated", "TCLUST", "TSRClust_unvalidated",
             "tEIGEN_BIC", "OTRIMLE_fixedK", "HDBSCAN_full", "DBSCAN_full", "Trimmed_kmeans")
headline_rows <- rbind(
  cbind(Study = "Architecture", arch[match(intersect(methods, arch$method), arch$method), c("method","MajorARI","MinorityClusterF1","AnyMicroDiscovery","EndToEndMedian")]),
  cbind(Study = "Stress", stress[match(intersect(methods, stress$method), stress$method), c("method","MajorARI","MinorityClusterF1","AnyMicroDiscovery","EndToEndMedian")]),
  cbind(Study = "Negative control", neg[match(intersect(methods, neg$method), neg$method), c("method","MajorARI","MinorityClusterF1","AnyMicroDiscovery","EndToEndMedian")])
)
names(headline_rows) <- c("Study","Method","MajorARI","MinorityClusterF1","AnyMicroDiscovery","EndToEndMedian")
write.csv(headline_rows, file.path(out, "tables", "headline_results.csv"), row.names = FALSE)
tab <- paste(c(
  "\\begin{tabular}{llrrrr}",
  "\\hline",
  "Study & Method & MajorARI & Minority F1 & Any/false micro & Median s \\\\",
  "\\hline",
  apply(headline_rows, 1, function(z) paste(latex_escape(z[1]), latex_escape(z[2]), fmt(as.numeric(z[3])), fmt(as.numeric(z[4])), fmt(as.numeric(z[5])), fmt(as.numeric(z[6])), sep = " & ") %+% " \\\\"),
  "\\hline",
  "\\end{tabular}"
), collapse = "\n")
write_table(
  "tables/table2_headline_results.tex", tab,
  "Headline simulation results from final replication-level exports.",
  "tab:headline-results",
  "Architecture excludes the N0 negative-control scenario and pools M1--M4; Negative control has $q=0$, so the Any/false micro column is a false micro-discovery probability. Oracle TCLUST uses true macro-plus-minority cluster count and true noise proportion and is not a practical unsupervised competitor. Runtime is end-to-end median seconds."
)

fit <- read_csv("github_results", "real_data_results", "dataset_fit_summary.csv")
tab <- paste(c(
  "\\begin{tabular}{lrrrrrrr}",
  "\\hline",
  "Dataset & $n$ & $p$ & $K$ & $\\alpha$ & $c$ & $|R_0|$ & Microclusters \\\\",
  "\\hline",
  apply(fit, 1, function(z) paste(latex_escape(gsub("_", " ", z[1])), z[2], z[3], z[4], fmt(as.numeric(z[5]), 2), z[6], z[7], z[8], sep = " & ") %+% " \\\\"),
  "\\hline",
  "\\end{tabular}"
), collapse = "\n")
write_table("tables/table3_realdata_settings.tex", tab, "Real-data dimensions and final TSRClust settings.", "tab:realdata-settings")

ari <- read_csv("github_results", "real_data_results", "external_macro_ari.csv")
resid <- read_csv("github_results", "real_data_results", "residual_external_eval.csv")
time <- read_csv("github_results", "real_data_results", "timings.csv")
real <- merge(merge(ari, resid, by = c("dataset", "method"), all = TRUE), time[, c("dataset", "method", "end_to_end_seconds")], by = c("dataset", "method"), all = TRUE)
keep <- real$method %in% c("TCLUST", "TSRClust_validated", "TSRClust_unvalidated", "DBSCAN_full", "HDBSCAN_full", "OTRIMLE_fixedK", "tEIGEN_BIC")
real <- real[keep, ]
tab <- paste(c(
  "\\begin{tabular}{llrrrr}",
  "\\hline",
  "Dataset & Method & Macro ARI & Residual ARI & Residual assigned & Seconds \\\\",
  "\\hline",
  apply(real, 1, function(z) paste(latex_escape(gsub("_", " ", z["dataset"])), latex_escape(z["method"]), fmt(as.numeric(z["MacroLabelARI"])), fmt(as.numeric(z["ResidualLocalARI"])), fmt(as.numeric(z["ResidualAssignedFraction"])), fmt(as.numeric(z["end_to_end_seconds"])), sep = " & ") %+% " \\\\"),
  "\\hline",
  "\\end{tabular}"
), collapse = "\n")
write_table("tables/table4_realdata_results.tex", tab, "Real-data external agreement, residual agreement, residual assignment, and runtime.", "tab:realdata-results")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  sim_plot <- rbind(
    cbind(Study = "Architecture M1--M4", arch),
    cbind(Study = "Stress", stress)
  )
  sim_plot <- sim_plot[sim_plot$method %in% methods, ]
  sim_plot$Study <- factor(sim_plot$Study, levels = c("Architecture M1--M4", "Stress"))
  label_map <- c(
    Oracle_TCLUST = "Oracle",
    TSRClust_validated = "TSR valid",
    TCLUST = "TCLUST",
    TSRClust_unvalidated = "TSR unval",
    tEIGEN_BIC = "tEIGEN",
    OTRIMLE_fixedK = "OTRIMLE",
    HDBSCAN_full = "HDBSCAN",
    DBSCAN_full = "DBSCAN",
    Trimmed_kmeans = "Trim. kmeans"
  )
  sim_plot$method_label <- unname(label_map[sim_plot$method])
  sim_plot$dx <- 0
  sim_plot$dy <- 0
  sim_plot$hjust <- 0.5
  sim_plot$vjust <- -0.7
  sim_plot$dx[sim_plot$method == "TCLUST"] <- -0.02
  sim_plot$hjust[sim_plot$method == "TCLUST"] <- 1
  sim_plot$dx[sim_plot$method == "TSRClust_validated"] <- 0.02
  sim_plot$hjust[sim_plot$method == "TSRClust_validated"] <- 0
  sim_plot$dy[sim_plot$method == "TSRClust_unvalidated"] <- -0.045
  sim_plot$vjust[sim_plot$method == "TSRClust_unvalidated"] <- 1
  sim_plot$dy[sim_plot$method == "Oracle_TCLUST"] <- 0.035
  sim_plot$dy[sim_plot$method == "tEIGEN_BIC"] <- -0.04
  sim_plot$vjust[sim_plot$method == "tEIGEN_BIC"] <- 1
  sim_plot$dy[sim_plot$method == "HDBSCAN_full"] <- -0.04
  sim_plot$vjust[sim_plot$method == "HDBSCAN_full"] <- 1
  p <- ggplot(sim_plot, aes(MajorARI, MinorityClusterF1, colour = Study)) +
    geom_point(size = 2.6, stroke = 0.8) +
    geom_text(aes(x = MajorARI + dx, y = MinorityClusterF1 + dy,
                  label = method_label, hjust = hjust, vjust = vjust),
              size = 2.35, colour = "black", lineheight = 0.9) +
    facet_wrap(~ Study, nrow = 1) +
    scale_x_continuous(limits = c(-0.08, 1.12)) +
    scale_y_continuous(limits = c(-0.04, 1.04)) +
    scale_colour_manual(values = c("Architecture M1--M4" = "#1B6CA8", "Stress" = "#B35C1E")) +
    theme_bw(base_size = 9) +
    theme(legend.position = "none",
          strip.background = element_rect(fill = "grey92", colour = "grey60")) +
    labs(x = "Major-cluster ARI", y = "Minority-cluster F1")
  ggsave(file.path(out, "figures", "fig2_macro_minor_tradeoff.pdf"), p, width = 7.2, height = 3.8)

  stress_cell <- read_csv("github_results", "summary_results", "stress_cell_summary.csv")
  st_methods <- c("TSRClust_validated", "TCLUST", "tEIGEN_BIC", "HDBSCAN_full", "OTRIMLE_fixedK")
  st <- stress_cell[stress_cell$method %in% st_methods & stress_cell$factor_changed != "baseline", ]
  st <- aggregate(cbind(MajorARI, MinorityClusterF1) ~ factor_changed + method, st, mean, na.rm = TRUE)
  st$method_label <- unname(label_map[st$method])
  st$method_label <- factor(st$method_label, levels = unname(label_map[st_methods]))
  factor_order <- aggregate(MinorityClusterF1 ~ factor_changed, st[st$method == "TSRClust_validated", ], mean)
  factor_order <- factor_order$factor_changed[order(factor_order$MinorityClusterF1)]
  st$factor_changed <- factor(st$factor_changed, levels = factor_order)
  p <- ggplot(st, aes(MinorityClusterF1, factor_changed, colour = method_label, shape = method_label)) +
    geom_point(size = 2.0, alpha = 0.92, position = position_dodge(width = 0.55)) +
    theme_bw(base_size = 8) +
    theme(legend.position = "bottom") +
    labs(x = "Mean minority-cluster F1", y = "Stress factor", colour = "Method", shape = "Method")
  ggsave(file.path(out, "figures", "fig3_stress_factors.pdf"), p, width = 6.8, height = 4.8)

  negp <- neg[neg$method %in% c("TSRClust_validated", "TSRClust_unvalidated", "tEIGEN_BIC", "OTRIMLE_fixedK", "HDBSCAN_full", "DBSCAN_full", "TCLUST"), ]
  p <- ggplot(negp, aes(reorder(method, AnyMicroDiscovery), AnyMicroDiscovery)) +
    geom_col(fill = "#3B6EA8", width = 0.7) +
    coord_flip() +
    theme_bw(base_size = 9) +
    labs(x = NULL, y = "False micro-discovery probability")
  ggsave(file.path(out, "figures", "fig4_negative_control.pdf"), p, width = 6.2, height = 3.5)

  runp <- rbind(cbind(Study = "Architecture", arch), cbind(Study = "Stress", stress), cbind(Study = "Negative control", neg))
  runp <- runp[runp$method %in% c("TCLUST", "TSRClust_validated", "tEIGEN_BIC", "OTRIMLE_fixedK", "DBSCAN_full", "HDBSCAN_full"), ]
  p <- ggplot(runp, aes(method, EndToEndMedian, fill = Study)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_y_log10() +
    coord_flip() +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom") +
    labs(x = NULL, y = "Median end-to-end seconds (log scale)")
  ggsave(file.path(out, "figures", "fig5_runtime.pdf"), p, width = 6.8, height = 4.0)

  olive <- read_csv("github_results", "real_data_results", "olive_microcluster_composition.csv")
  shuttle <- read_csv("github_results", "real_data_results", "shuttle_microcluster_composition.csv")
  comp <- rbind(olive[, c("dataset", "cluster", "size", "dominant_class", "purity")], shuttle[, c("dataset", "cluster", "size", "dominant_class", "purity")])
  p <- ggplot(comp, aes(factor(cluster), size, fill = dominant_class)) +
    geom_col(width = 0.7) +
    facet_wrap(~ dataset, scales = "free_x") +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom") +
    labs(x = "Accepted TSRClust microcluster", y = "Size")
  ggsave(file.path(out, "figures", "fig6_realdata_microclusters.pdf"), p, width = 6.8, height = 3.8)
}

neg_raw <- read.csv(gzfile(file.path(root, "github_results", "raw_results", "negative_control_replications.csv.gz")), check.names = FALSE)
neg_tsr <- neg_raw[neg_raw$method == "TSRClust_validated" & neg_raw$status == "OK", ]
neg_false <- sum(neg_tsr$AnyMicroDiscovery > 0, na.rm = TRUE)
neg_total <- nrow(neg_tsr)
neg_ci <- binom.test(neg_false, neg_total)$conf.int
mc_uncertainty <- data.frame(
  method = "TSRClust_validated",
  event = "false_micro_discovery",
  count = neg_false,
  denominator = neg_total,
  estimate = neg_false / neg_total,
  exact95_low = neg_ci[1],
  exact95_high = neg_ci[2]
)
write.csv(mc_uncertainty, file.path(out, "tables", "negative_control_ci.csv"), row.names = FALSE)

supp <- c(
  "\\documentclass[pdflatex,sn-nature]{sn-jnl}",
  "\\usepackage{amsmath,amssymb,amsfonts}",
  "\\usepackage{booktabs}",
  "\\usepackage{graphicx}",
  "\\begin{document}",
  "\\title{Supplementary Information for TSRClust}",
  "\\author{Hasan Bulut}",
  "\\maketitle",
  "\\section*{Supplementary Methods}",
  "The reproducibility files contain the complete simulation manifests and replication-level result exports. The architecture study used 80 cells and 35,200 generated datasets; the stress study used 30 cells and 9,000 generated datasets; and the expanded negative-control study used 10 cells and 10,000 generated datasets.",
  "\\section*{S1. Complete simulation design}",
  "Architecture, stress, and negative-control cell-summary CSV files provide compact cell-level summaries. Full replication-level results are provided with the reproducibility package rather than reproduced as large tables here.",
  "\\section*{S2. Paired comparisons}",
  "Architecture, stress, and negative-control paired-comparison CSV files contain method differences against validated TSRClust with Monte Carlo standard errors.",
  "\\section*{S3. Monte Carlo uncertainty}",
  sprintf("For expanded negative controls, validated TSRClust produced %d false micro-discoveries in %d successful evaluations, giving an estimated probability of %.4f with an exact 95\\%% binomial confidence interval of [%.5f, %.5f].", neg_false, neg_total, neg_false / neg_total, neg_ci[1], neg_ci[2]),
  "\\section*{S4. OTRIMLE sensitivity and runtime}",
  "OTRIMLE sensitivity outputs and detailed runtime summaries are included in the reproducibility results folders. The main text reports the fixed-\\(K\\) OTRIMLE results used in the headline comparison.",
  "\\section*{S5. Real-data details}",
  "Real-data supplementary CSV files provide class-specific trimming rates and accepted residual micro-cluster compositions for the olive-oil and Shuttle analyses.",
  "\\begin{figure}[htbp]",
  "\\centering",
  "\\includegraphics[width=\\linewidth]{../figures/fig6_realdata_microclusters.pdf}",
  "\\caption{Accepted real-data TSRClust micro-clusters. Bar heights show accepted residual micro-cluster sizes; fill identifies the dominant external class.}",
  "\\label{fig:supp-realdata-microclusters}",
  "\\end{figure}",
  "\\end{document}"
)
writeLines(supp, file.path(out, "supplementary", "supplementary_information.tex"))

for (f in c("architecture_cell_summary.csv", "architecture_paired_vs_TSRvalidated.csv", "stress_cell_summary.csv", "stress_paired_vs_TSRvalidated.csv", "negative_control_cell_summary.csv", "negative_control_paired_vs_TSRvalidated.csv")) {
  file.copy(file.path(root, "github_results", "summary_results", f), file.path(out, "supplementary", f), overwrite = TRUE)
}
for (f in c("olive_microcluster_composition.csv", "shuttle_microcluster_composition.csv", "class_trimming_rates.csv")) {
  file.copy(file.path(root, "github_results", "real_data_results", f), file.path(out, "supplementary", f), overwrite = TRUE)
}

write.csv(data.frame(
  item = c("architecture_cells", "architecture_datasets", "stress_cells", "stress_datasets", "negative_control_cells", "negative_control_datasets", "total_method_evaluations", "negative_control_false_micro_count", "negative_control_false_micro_denominator", "negative_control_false_micro_probability", "negative_control_false_micro_exact95_low", "negative_control_false_micro_exact95_high"),
  value = c(manifest_summary$cells[1], manifest_summary$generated_datasets[1],
            manifest_summary$cells[2], manifest_summary$generated_datasets[2],
            manifest_summary$cells[3], manifest_summary$generated_datasets[3],
            sum(manifest_summary$method_evaluations),
            neg_false, neg_total, neg_false / neg_total, neg_ci[1], neg_ci[2])
), file.path(out, "tables", "computed_claims.csv"), row.names = FALSE)

cat("ASSETS_OK\n")
