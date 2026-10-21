# BE700 Homework 3: Insulin Resistance and Aging (human muscle, GSE25462)

Analysis report written in R Markdown, rendered to PDF in the cloud. Nothing needs to run locally.

| Path | What |
| --- | --- |
| `HW3.Rmd` | Full analysis + report (QC, heatmaps, clustering, phenotype enrichment, age-correlated genes, PCA, SI / T2D, DAVID gene lists) |
| `data/` | Expression matrix (22,189 genes x 50 samples, log-normalized) and GEO phenotype table |
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

Upload `output/results/DAVID_age_down.txt` (and the other `DAVID_*.txt` lists) to https://davidbioinformatics.nih.gov with `DAVID_background.txt` as background; add the top KEGG / GO terms to the DAVID section of `HW3.Rmd`.
