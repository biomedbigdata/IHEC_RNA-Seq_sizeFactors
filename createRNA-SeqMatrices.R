library(data.table)
library(R.utils)
library(DESeq2)
library(pbmcapply)
options(mc.cores = 10)

rna_seq_dir <- "/nfs/data/IHEC/RNAseq/RNA-Seq"
iso_files <- list.files(rna_seq_dir, "\\.isoforms\\.results$", full.names=TRUE)
gene_files <- list.files(rna_seq_dir, "\\.genes\\.results$", full.names=TRUE)
epirr.uuid_pattern <- "ihec.rna-seq.ihec-grapenf-containerv1.1.0.(IHECRE[[:alnum:]]{8}\\.[[:alnum:]]+\\.[[:alnum:]]{8}-[[:alnum:]]{4}-[[:alnum:]]{4}-[[:alnum:]]{4}-[[:alnum:]]{12})"

write_sample_matrices <- function(file_list, id_col, outfile_prefix, count_col="expected_count", tpm_col="TPM"){
  # read all samples
  sample_list <- pbmcapply::pbmclapply(
    file_list,
    fread,
    stringsAsFactors = TRUE,
    select = c(id_col, count_col, tpm_col)
  )
  # name sample list by epirr id and uuid
  sample_names <- sub(
    "\\.(genes|isoforms)\\.results$",
    "",
    sub(
      "ihec.rna-seq.ihec-grapenf-containerv1.1.0.",
      replacement = "",
      basename(file_list),
      fixed = TRUE
    )
  )
  names(sample_list) <- sample_names
  
  # gather samples into dt
  sample_dt <-
    data.table::rbindlist(sample_list, idcol = "EpiRR.uuid")
  
  # write wide TPM matrix
  data.table::fwrite(data.table::dcast(sample_dt, get(id_col) ~ EpiRR.uuid, value.var = tpm_col), file = paste0(outfile_prefix, "_TPM.csv.gz"))
  
  # gather count data
  counts_wide <- data.table::dcast(sample_dt, get(id_col) ~ EpiRR.uuid, value.var = count_col)
  
  # estimate size factors for count data
  sizeFactors <- DESeq2::estimateSizeFactorsForMatrix(counts_wide[, ..sample_names])
  
  # multiply counts with size factors
  counts_wide[, (sample_names):=Map("*", .SD, sizeFactors), .SDcols=sample_names]
  data.table::fwrite(counts_wide, file = paste0(outfile_prefix, "_DESeq2Counts.csv.gz"))
}

write_sample_matrices(file_list = gene_files, id_col = "gene_id", outfile_prefix = "genes")
write_sample_matrices(file_list = gene_files, id_col = "gene_id", count_col = "posterior_mean_count", tpm_col = "pme_TPM", outfile_prefix = "genes_posterior")

write_sample_matrices(file_list = iso_files, id_col = "transcript_id", outfile_prefix = "isoforms")
write_sample_matrices(file_list = iso_files, id_col = "transcript_id", count_col = "posterior_mean_count", tpm_col = "pme_TPM", outfile_prefix = "isoforms_posterior")
