# case_study_bioconductor_deseq2

Case study of Bioconductor for CBIO7501, with a DESeq2 differential expression analysis of the airway dataset (dexamethasone-treated vs untreated human airway smooth muscle cells).

## Run

From the project folder, with R, DESeq2 and airway installed:

```
Rscript scripts/deseq2_demo.r
```

## Data

- Source: Bioconductor package `airway`. See `metadata/README.md`.
- Filtering: genes with count >= 10 in at least 4 samples.

## Folders

| Folder | Contents |
|---|---|
| `.github/workflows/` | CI workflows (none yet) |
| `config/` | Analysis parameters (none yet) |
| `containers/` | Container definitions (not used yet) |
| `data/raw/` | Read-only downloaded data (empty; data comes from the package) |
| `data/processed/` | Derived, rebuildable data (empty) |
| `docs/` | Extra documentation (empty) |
| `envs/` | Software environment files |
| `metadata/` | Data sources and abbreviations |
| `notebooks/` | Exploratory notebooks (not used yet) |
| `reports/` | Homework report |
| `results/` | Tables and figures made by the script |
| `scripts/` | Analysis code |
| `workflow/` | Pipeline definition (not used yet) |
