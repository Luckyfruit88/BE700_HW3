# BE700 Homework 3: Insulin Resistance and Aging

Exploratory analysis of the course-provided GSE25462 human muscle expression and phenotype tables. `HW3.Rmd` contains the analysis and interpretation. Runtime session information is saved separately in `output/sessionInfo.txt`.

## Run and validate

From the repository root, in an R environment with `rmarkdown`, `knitr`, `pheatmap`, `ggplot2`, `ggrepel`, and `cluster`, plus Pandoc and LaTeX:

```bash
Rscript --vanilla scripts/render_validate.R
```

On BU SCC, load the modules first:

```bash
module load R/4.4.3 pandoc/2.5 texlive/2023
Rscript --vanilla scripts/render_validate.R
```

The default output is `output/HW3.pdf`; tables are copied into `output/results/`. `output/validation.txt` records numerical and cache-validity checks, and `output/sessionInfo.txt` records the actual runtime. An alternative output directory may be supplied as the first script argument. Preserve an existing run before overwriting its outputs.

In VS Code, use the R extension's **Run Chunk** or **Run Above** for interactive work, **Open Preview to the Side** for an HTML preview, and **Knit Rmd** for the document's PDF output. Previewing renders in a separate R process; its objects are not populated into the interactive console. Final verification should use the command above in a fresh process, rather than depend on objects left by interactive chunks.

GitHub Actions runs the same render-and-validation script, uploads the report, and commits validated outputs. Changes to the Rmd, data, validation script, or workflow trigger it. Output-only bot commits do not retrigger rendering.

## Interpretation and validation boundaries

- The 50-sample dataset contains 22,189 annotated expression features, including multi-gene annotations. Course instructions identify it as already log-normalized. Original GEO processing and the course matrix's downstream transformation are described separately in the report.
- The cohort SI summary uses available observations. The T2D mean is 1.815 from 8 observations.
- Full-sample Pearson testing has no 5%-FDR age hits; UNC13C passes Spearman and age/sex/diabetes-adjusted testing. Removing the four exploratory QC flags changes the adjusted-model hit to PLAG1 (q about 0.013); UNC13C no longer passes 5% FDR.
- PC--phenotype results include raw p-values, BH adjustment within each 42-test PCA family, and adjustment across all four analyses. Nominal associations are not described as FDR-significant.
- The matched PCA pair uses the same feature selection and scaling rule, changing only sample inclusion. QC flags do not establish mislabeling or contamination.
- These are exploratory, within-cohort analyses. Expression associations and pathway annotation do not establish causal aging effects, SRF activity, or independent replication.

The validation script checks the stated dataset contract, independent `lm()`/`cor.test()` agreement, the PLAG1/UNC13C sensitivity result, missing-value summaries, PCA multiplicity, separate age/DAVID objects, and rejection of stale DAVID inputs. Passing these checks establishes implementation consistency, not biological validation.

## DAVID input provenance

`data/david/DAVID_chart_results.tsv` and `DAVID_run_info.tsv` preserve the imported October 21, 2026 DAVID results. Five normalized symbol lists are frozen in `data/david/input_symbols/`; `input_manifest.tsv` records their SHA-256 hashes, symbol counts, and the source commit from which the associated original lists were recovered.

The Rmd splits `///`, trims whitespace, normalizes symbol case, and compares current input membership against these frozen lists before interpreting cached results. A changed threshold, cluster result, or gene selection can therefore stop a render with a stale-DAVID error. In that case, rerun DAVID with the new inputs and update the frozen lists, manifest, run metadata, and imported results together. Do not replace the frozen lists merely to suppress the check. Rerendering the Rmd does not submit a new DAVID analysis.

The saved run used OFFICIAL_GENE_SYMBOL, Homo sapiens, the Human Genome U133 Plus 2 Array background, KEGG_PATHWAY and GOTERM_BP_DIRECT, Count >= 2, and EASE p <= 0.1. DAVID's `FDR` column is adaptive FDR; `Benjamini` is ordinary BH. Both identify the same eight age-down terms in this saved run.

## Files

| Path | Contents |
| --- | --- |
| `HW3.Rmd` | Editable report and analysis |
| `data/` | Course expression/phenotype inputs and saved DAVID evidence |
| `docs/` | Assignment handout and Hallmarks of Aging paper |
| `scripts/render_validate.R` | Full render and regression validation |
| `output/HW3.pdf` | Latest validated PDF |
| `output/results/` | Full result tables and DAVID input exports |
| `output/validation.txt` | Validation checks and selected numeric receipts |
| `output/sessionInfo.txt` | Full session information |
