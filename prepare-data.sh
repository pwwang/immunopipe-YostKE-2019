#!/usr/bin/env bash

# The URL to the tar file
DATAURL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123813/suppl/"
DATAFILE1="prepared-data/GSE123813_bcc_all_metadata.txt.gz"
DATAFILE2="prepared-data/GSE123813_bcc_scRNA_counts.txt.gz"
DATAFILE3="prepared-data/GSE123813_bcc_tcell_metadata.txt.gz"
DATAFILE4="prepared-data/GSE123813_bcc_tcr.txt.gz"

echo "- Make the directory for prepared data ..."
mkdir -p prepared-data

echo "- Download the data if needed ..."
if [ ! -e $DATAFILE1 ]; then
    wget -O $DATAFILE1 $DATAURL"GSE123813_bcc_all_metadata.txt.gz"
fi
if [ ! -e $DATAFILE2 ]; then
    wget -O $DATAFILE2 $DATAURL"GSE123813_bcc_scRNA_counts.txt.gz"
fi
if [ ! -e $DATAFILE3 ]; then
    wget -O $DATAFILE3 $DATAURL"GSE123813_bcc_tcell_metadata.txt.gz"
fi
if [ ! -e $DATAFILE4 ]; then
    wget -O $DATAFILE4 $DATAURL"GSE123813_bcc_tcr.txt.gz"
fi

echo "- Unzip the data ..."
gunzip -kfc $DATAFILE3 > prepared-data/GSE123813_bcc_tcell_metadata.txt

echo "- Create Seurat object ..."
cat << EOF > /tmp/immunopipe-YostKE-2019-create-seurat-object.R
library(Seurat)
library(rlang)
library(stringr)
library(data.table)
library(biopipen.utils)

to_airr <- function(df) {
  dt <- as.data.table(df)

  # Split chains
  dt[, nt := strsplit(cdr3s_nt, ";", fixed = TRUE)]
  dt[, aa := strsplit(cdr3s_aa, ";", fixed = TRUE)]

  airr <- dt[
    ,
    .(
      nt = unlist(nt, use.names = FALSE),
      aa = unlist(aa, use.names = FALSE)
    ),
    by = .(cell_id = barcode)
  ]

  airr[, c("locus", "junction") := tstrsplit(nt, ":", fixed = TRUE)]
  airr[, c("tmp", "junction_aa") := tstrsplit(aa, ":", fixed = TRUE)]
  airr[, tmp := NULL]

  airr[
    ,
    \`:=\`(
      consensus_count = 1L,
      v_call = NA_character_,
      d_call = NA_character_,
      j_call = NA_character_,
      c_call = NA_character_
    )
  ]

  setcolorder(
    airr,
    c(
      "cell_id",
      "locus",
      "consensus_count",
      "v_call",
      "d_call",
      "j_call",
      "c_call",
      "junction",
      "junction_aa"
    )
  )

  airr[]
}

print("Converting TCR data to AIRR format ...")
tcr_data <- read.table(gzfile("${DATAFILE4}"), sep="\t", header=TRUE, row.names=1)
under_score_samples <- tcr_data\$barcode <- rownames(tcr_data)
under_score_samples <- gsub("\\\\.", "_", under_score_samples)
samples <- sub("_tcell_.*$", "", under_score_samples)
samples <- sub("_[ATGC]+$", "", samples)
samples <- unique(sub("^bcc_", "", samples))
for (sample in samples) {
    print(paste0("- Sample: ", sample))
    td <- tcr_data[grepl(sample, under_score_samples), , drop=FALSE]
    airr_data <- to_airr(td)
    dir.create(file.path("prepared-data", sample), showWarnings=FALSE, recursive=TRUE)
    write.table(airr_data, file=file.path("prepared-data", sample, "airr_rearrangement.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
}

count_file <- "${DATAFILE2}"
metadata_file <- "${DATAFILE1}"

print("Creating Seurat object ...")
counts <- read.table(gzfile(count_file), sep="\t", header=TRUE, row.names=1)
metadata <- read.table(gzfile(metadata_file), sep="\t",header=TRUE, row.names=NULL)
rownames(metadata) <- metadata\$cell.id
obj <- CreateSeuratObject(counts=counts, meta.data=metadata)
obj\$Sample <- paste0(obj\$patient, "_", obj\$treatment)
write_obj(obj, "prepared-data/seurat-object.qs")
EOF

Rscript /tmp/immunopipe-YostKE-2019-create-seurat-object.R

echo "- Done!"
