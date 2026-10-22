# Run from the repository root: Rscript scripts/render_validate.R [output-dir]
# Uses a fresh evaluation environment and validates scientific regression checks.
args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[[1]] else "output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
input_before <- tools::md5sum(c("HW3.Rmd", list.files("data", recursive = TRUE,
                                                   full.names = TRUE)))
run_env <- new.env(parent = globalenv())
pdf <- rmarkdown::render("HW3.Rmd", output_dir = output_dir, envir = run_env,
                        quiet = FALSE)
checks <- character()
check <- function(name, value) {
  if (!isTRUE(value)) stop("VALIDATION FAILED: ", name)
  checks <<- c(checks, name)
}
e <- run_env
check("50 aligned samples and 22189 finite expression features",
  ncol(e$expr) == 50 && nrow(e$expr) == 22189 && all(is.finite(e$expr)) &&
    identical(colnames(e$expr), e$ph$sample))
t2d <- e$cohort_summary[e$cohort_summary$group == "T2D", ]
check("T2D SI mean uses 8 observations and equals 1.815",
  t2d$si_n == 8 && abs(t2d$si_mean - 1.815) < 1e-10)
check("age counts retain numeric type after DAVID execution",
  is.numeric(e$age_sig_counts) &&
    identical(unname(e$age_sig_counts[c("pearson_q05", "spearman_q05", "adjusted_q05")]),
              c(0L, 1L, 1L)) && is.data.frame(e$david_sig_terms) &&
    !exists("sig", envir = e, inherits = FALSE))
check("full-sample 5 percent adjusted age hit is UNC13C",
  identical(e$age_full_significant, "UNC13C"))
check("46-sample 5 percent adjusted age hit is PLAG1",
  identical(e$age_sens_significant, "PLAG1"))
plag_q <- e$age_sensitivity$q[match("PLAG1", e$age_sensitivity$gene)]
unc_q <- e$age_sensitivity$q[match("UNC13C", e$age_sensitivity$gene)]
check("PLAG1 q is about 0.013; UNC13C no longer passes 5 percent",
  abs(plag_q - 0.013) < 0.001 && unc_q > 0.05 && unc_q < 0.1)
# Independent lm() comparisons exercise the QR coefficient and SE extraction.
genes <- c("UNC13C", "PLAG1", "MFF", "SRF")
for (g in genes) {
  fit <- lm(e$expr[g, ] ~ age + sex + diabetes, data = e$ph)
  observed <- summary(fit)$coefficients["age", ]
  actual <- e$age_tab[match(g, e$age_tab$gene), ]
  check(paste("OLS matches lm for", g),
    abs(observed[[1]] - actual$beta_age_adj) < 1e-10 &&
      abs(observed[[4]] - actual$p_adj) < 1e-10)
  expected_sp <- cor.test(e$expr[g, ], e$ph$age,
                          method = "spearman", exact = FALSE)$p.value
  check(paste("Spearman matches cor.test for", g),
    abs(expected_sp - actual$p_spearman) < 1e-10)
}
for (g in c("PLAG1", "UNC13C")) {
  fit <- lm(e$expr[g, e$keep_s] ~ age + sex + diabetes,
            data = e$ph[e$keep_s, ])
  actual <- e$age_sensitivity[match(g, e$age_sensitivity$gene), ]
  check(paste("sensitivity OLS matches lm for", g),
    abs(summary(fit)$coefficients["age", 4] - actual$p) < 1e-10)
}
pc <- e$pc_associations
check("four complete PCA test families contain 42 tests each",
  nrow(pc) == 168 && all(table(pc$analysis) == 42))
si <- subset(pc, analysis == "matched_46" & PC == "PC2" & phenotype == "si")
check("nominal PC2-SI association is not significant after 42-test BH",
  nrow(si) == 1 && abs(si$p - 0.0177882) < 1e-5 &&
    abs(si$q_within - 0.2037034) < 1e-5 && si$q_all > 0.05)
check("matched PCA pair fixes features and scaling rule",
  identical(rownames(e$pca_matched$rotation), rownames(e$pca_s$rotation)) &&
    is.numeric(e$pca_matched$scale) && is.numeric(e$pca_s$scale))
check("frozen DAVID lists match the actual current inputs",
  isTRUE(e$check_david_inputs(e$david_current_symbols, e$david_frozen_symbols)))
changed <- e$david_current_symbols
changed$age_down <- c(changed$age_down, "INTENTIONAL_MISMATCH_FOR_VALIDATION")
check("DAVID cache guard rejects altered membership",
  inherits(try(e$check_david_inputs(changed, e$david_frozen_symbols), silent = TRUE),
           "try-error"))
check("DAVID cache guard rejects missing or extra lists",
  inherits(try(e$check_david_inputs(changed[-1], e$david_frozen_symbols),
               silent = TRUE), "try-error"))
check("adaptive and ordinary BH identify the same eight saved terms",
  nrow(e$david_sig_terms) == 8 &&
    identical(which(e$dv$FDR < 0.05), which(e$dv$Benjamini < 0.05)))
check("source and supplied data were not changed by rendering",
  identical(input_before, tools::md5sum(names(input_before))))
check("a nonempty PDF was generated", file.exists(pdf) && file.info(pdf)$size > 10000)
dir.create(file.path(output_dir, "results"), showWarnings = FALSE)
files <- list.files("results", full.names = TRUE)
check("result tables copied beside the PDF",
  all(file.copy(files, file.path(output_dir, "results"), overwrite = TRUE)))
writeLines(c(paste("PASS", checks, sep = "\t"),
  paste("PLAG1_adjusted_q_46", format(plag_q, digits = 12), sep = "\t"),
  paste("UNC13C_adjusted_q_46", format(unc_q, digits = 12), sep = "\t"),
  paste("R_version", R.version.string, sep = "\t")),
  file.path(output_dir, "validation.txt"))
writeLines(sub("[[:space:]]+$", "", capture.output(sessionInfo())),
           file.path(output_dir, "sessionInfo.txt"))
cat("VALIDATION_PASS:", length(checks), "checks\n")
