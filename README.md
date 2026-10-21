# BE700 Homework 3: Insulin Resistance and Aging (human muscle, GSE25462)

Analysis report written in R Markdown, rendered to PDF in the cloud. Nothing needs to run locally.

| Path | What |
| --- | --- |
| `HW3.Rmd` | Full analysis + report (QC, heatmaps, clustering, phenotype enrichment, age-correlated genes, PCA, SI / T2D, DAVID enrichment) |
| `data/` | Expression matrix (22,189 genes x 50 samples, log-normalized) and GEO phenotype table |
| `data/david/` | DAVID Functional Annotation Chart output (KEGG + GO-BP) and run info, read by `HW3.Rmd` |
| `docs/` | Assignment handout and Hallmarks of Aging paper |
| `output/HW3.pdf` | Rendered report (written by GitHub Actions) |
| `output/results/` | Gene tables and DAVID input lists |

## How to run

**Automatic:** every push that changes `HW3.Rmd` or `data/` triggers the *Render HW3 report* workflow, which knits the PDF in a `rocker/verse` container and commits it to `output/`. It can also be started by hand from the Actions tab (Run workflow).

**Codespaces:** Code -> Codespaces -> Create codespace. The dev container has R 4.4 + LaTeX. Then in the terminal:

```bash
Rscript -e 'rmarkdown::render("HW3.Rmd", output_dir = "output")'
```

## Part III (DAVID)

Done on 2026-10-21. The five lists in `output/results/DAVID_*.txt` (age_down, age_up, cluster1_up, cluster2_up, SI_top200) were run through the DAVID Functional Annotation Chart (OFFICIAL_GENE_SYMBOL, Homo sapiens, background = built-in *Human Genome U133 Plus 2 Array*, KEGG_PATHWAY + GOTERM_BP_DIRECT, Count >= 2, EASE <= 0.1). Results are in `data/david/DAVID_chart_results.tsv` and summarized in the DAVID section of the report. Only `age_down` has terms at FDR < 0.05 (carbon metabolism, TCA cycle, malate-aspartate shuttle, mitochondrial fission).
