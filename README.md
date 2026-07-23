# Single-cell transcriptomic analysis of acute chikungunya virus infection

This repository contains the R analysis scripts associated with the manuscript
describing peripheral immune responses in patients with acute chikungunya virus
(CHIKV) infection. The study profiles peripheral blood mononuclear cells from
29 patients with acute CHIKV infection and 10 healthy controls and includes a
comparison with a public dengue virus (DENV) single-cell dataset.

## Repository structure

| Directory | Analysis |
| --- | --- |
| `01_initial_clustering_and_roe` | Quality control, Seurat integration, major-cell annotation, and observed/expected cell-type ratios |
| `02_ifna_and_inflammation_scores` | IFN-alpha response and acute inflammatory-response gene-set scoring |
| `03_innate_immune_cells` | Innate immune-cell subclustering, differential expression, enrichment, migration, and apoptosis analyses |
| `04_innate_inflammatory_response` | Innate-cell cytokine and inflammatory-response scoring |
| `05_b_cell_subclustering` | B-cell subclustering, differential expression, immunoglobulin, and inflammatory analyses |
| `06_t_cell_subclustering` | T-cell subclustering, differential expression, cytotoxicity, exhaustion, migration, and apoptosis analyses |
| `07_global_cellchat_analysis` | CellChat inference and comparison of cell-cell communication networks |
| `08_supplementary_analyses` | Supplementary B-cell and CD8 T-cell summaries |
| `09_ifn_source_analysis` | Comparison of interferon expression across cell types |
| `10_chikv_dengue_comparison` | Integration and comparison of CHIKV and public DENV single-cell datasets |

The numbered directories follow the approximate order of the analyses in the
manuscript. Each directory contains one R script.

## Data availability

The CHIKV sequencing data were submitted to the Genome Sequence Archive for
Human under accession **HRA017257**. Human sequencing data may require
controlled-access approval from the repository before download.

The DENV comparison uses the public dataset **GSE220969**.

Participant-level raw or processed data are not stored in this GitHub
repository. Before running the scripts, place authorized inputs under a local
`data/` directory. At minimum, later-stage scripts expect a processed Seurat
object at:

```text
data/combined_seurat.rds
```

The initial-processing script expects Cell Ranger count matrices under
`data/Matrix/`. The DENV comparison script expects downloaded DENV matrices
under `data/dengue/`.

## Software

The analyses were developed in R using established single-cell and
bioinformatics packages:

- Seurat 5
- Matrix
- here
- CellChat
- UCell
- msigdbr
- Harmony
- clusterProfiler, enrichplot, AnnotationDbi, GO.db, and org.Hs.eg.db
- dplyr, tidyr, readr, stringr, ggplot2, ggpubr, ggrepel, pheatmap, and
  RColorBrewer

For strict reproducibility, record the exact R and package versions used for
the final figures with `sessionInfo()` or `renv::snapshot()` before archival.
The repository does not currently contain a lockfile, so package updates may
change numerical or graphical results.

## Running the analysis

1. Clone the repository and open R in the repository root.
2. Create `data/` and obtain the controlled-access CHIKV data through the
   authorized repository process.
3. Install the packages listed above from CRAN, Bioconductor, or the package
   authors' official installation instructions.
4. Run the numbered scripts in order where dependencies exist. Several
   downstream scripts require objects or tables produced by earlier scripts.
5. Write generated figures and intermediate files to a local `results/`
   directory or update the output paths for your environment.

The scripts document the analysis workflow but are not a fully automated
pipeline. Input object metadata must contain the fields referenced by each
script, including sample, condition, and cell-type annotations. Review these
field names before execution.

## Statistical considerations

Biological replication is at the participant or sample level, not the cell
level. Cell-level visualizations should therefore be accompanied by
sample-level summaries or models that account for participant-level
replication. Differential-expression and abundance results should be
interpreted with this design in mind.

## Reviewer access to this private repository

The editorial submission includes a repository-specific, read-only GitHub
deploy key. It permits cloning over SSH but does not provide access to the
GitHub web interface and does not permit changes to the repository.

On macOS or Linux:

```bash
chmod 600 chikv_reviewer_key
GIT_SSH_COMMAND="ssh -i /absolute/path/chikv_reviewer_key -o IdentitiesOnly=yes" \
  git clone git@github.com:Marxxxxxxxxxx/CHIKV-scRNAseq.git
```

On Windows, restrict the private-key file so that it is readable only by the
current user before using the same SSH clone URL. If access fails, contact the
corresponding author through the journal rather than sharing the key through an
unsecured channel.

The deploy key should be removed from the repository settings after peer review
is complete. A permanent public repository or a DOI-backed software archive is
preferable for the final published version.

## Citation

Please cite the associated manuscript when using this code. Full citation
details will be added after publication.

## License

No software license has yet been declared. Until a license is added, copyright
remains with the repository owner and reuse permission is not automatically
granted.
