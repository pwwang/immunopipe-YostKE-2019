# immunopipe-Yost-2019

Reanalysis of the scRNA-seq and scTCR-seq data from [Yost, K.E., et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6689255) using [immunopipe](https://github.com/pwwang/immunopipe).

> [Yost, Kathryn E., et al. "Clonal replacement of tumor-specific T cells following PD-1 blockade." Nature medicine 25.8 (2019): 1251-1259.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6689255)

## Data preparation

The data was downloaded from [GSE123813](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123813). Only basal cell carcinoma patients were included in the analysis, where 11 samples from 6 patients were sequenced with scRNA-seq and scTCR-seq. 21 samples from 11 patients were included with scTCR-seq data and their paired scRNA-seq data.

See `prepare-data.sh` for details.

## Configuration

> [!NOTE]
> This is not a replication of the original paper, primarily due to the irreproducibility of the clustering results. This is a reanalysis of the data using [`immunopipe`](https://github.com/pwwang/immunopipe), showing the potential of the pipeline similar analyses listed in the paper.
>

This examples shows how to load scRNA-seq data from a Seurat object and scTCR-seq data from sample-specific AIRR files.

See details at `Immunopipe.config.toml`.

## Results/Reports

You can find the results in the `Immunopipe-output` directory.

The report can be found at [https://imp-yostke-2020.pwwang.com/REPORTS](https://imp-yostke-2020.pwwang.com/REPORTS).
