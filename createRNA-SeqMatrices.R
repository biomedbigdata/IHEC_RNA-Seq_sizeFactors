library(data.table)
library(R.utils)
library(DESeq2)
library(pbmcapply)
library(ggfortify)
library(ggrepel)
options(mc.cores = 10)

rna_seq_dir <- "/nfs/data/IHEC/RNAseq/RNA-Seq"
iso_files <- list.files(rna_seq_dir, "\\.isoforms\\.results$", full.names=TRUE)
gene_files <- list.files(rna_seq_dir, "\\.genes\\.results$", full.names=TRUE)
epirr.uuid_pattern <- "ihec.rna-seq.ihec-grapenf-containerv1.1.0.(IHECRE[[:alnum:]]{8}\\.[[:alnum:]]+\\.[[:alnum:]]{8}-[[:alnum:]]{4}-[[:alnum:]]{4}-[[:alnum:]]{4}-[[:alnum:]]{12})"

metadata <- data.table::fread('IHEC_metadata_harmonization.v1.0.csv')
metadata[, EpiRR_no_version := data.table::tstrsplit(EpiRR, '.', fixed = TRUE)[1]]
outliers <- metadata[harmonized_sample_ontology_intermediate == "germ line cell"]

write_dt_and_md5 <- function(dt, out_csv){
  data.table::fwrite(dt, file = out_csv)
  md5 <- tools::md5sum(out_csv)
  writeLines(text = paste(md5, names(md5), sep="  "), con=paste0(out_csv, ".md5"))
}

pca_plot <- function(dt, outfile_prefix) {
  matrix <- t(dt[, !"id_col"])
  
  # filter out germ line samples
  matrix <- matrix[!Reduce("|" , lapply(outliers[, EpiRR_no_version], function(prefix) startsWith(rownames(matrix), prefix))), ]
  
  # filter out cols without variance
  matrix <- matrix[, colVars(matrix) != 0, drop = FALSE]

  # perform PCA
  pca <- irlba::prcomp_irlba(matrix, n=2, center = TRUE, scale=TRUE)
  
  metadata_tmp <- metadata[sapply(data.table::tstrsplit(rownames(matrix), ".", fixed = TRUE, keep = 1)[[1]], grep, x = EpiRR_no_version, fixed=TRUE)]
  labels <- metadata_tmp[, harmonized_sample_ontology_intermediate]
  ggplot2::autoplot(
    pca,
    data = metadata_tmp,
    colour = "project",
    shape = "harmonized_biomaterial_type"
  ) + 
    scale_shape_manual(values=1:metadata_tmp[, uniqueN(project)]) +
    labs(title = outfile_prefix, color = 'Project', shape = 'Biomaterial Type') + 
    theme_bw() + guides(col=guide_legend(ncol = 1)) + 
    ggrepel::geom_text_repel(aes(label = labels), min.segment.length = 0, max.overlaps = 5, force = 5)
}



write_sample_matrices <- function(file_list, id_col, outfile_prefix, count_cols=c("expected_count", "posterior_mean_count"), tpm_cols=c("TPM", "pme_TPM")){
  # read all samples
  sample_list <- pbmcapply::pbmclapply(
    file_list,
    data.table::fread,
    stringsAsFactors = TRUE,
    select = c(id_col, count_cols, tpm_cols)
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
  # sample_dt[, EpiRR_no_version := data.table::tstrsplit(EpiRR.uuid, '.', fixed = TRUE)[1]]
  
  # write wide TPM matrix
  for (tpm_col in tpm_cols) {
    outfile_col <- paste0(outfile_prefix, "_", tpm_col)
    tpm_wide <- data.table::dcast(sample_dt, get(id_col) ~ EpiRR.uuid, value.var = tpm_col)
    write_dt_and_md5(tpm_wide, paste0(outfile_col, ".csv.gz"))
    
    ggsave(paste0(outfile_col, ".pdf"), pca_plot(tpm_wide, outfile_col), width = 12, height = 10)
  }
  
  # gather count data
  for (count_col in count_cols) {
    outfile_col <- paste0(outfile_prefix, "_", count_col)
    counts_wide <- data.table::dcast(sample_dt, get(id_col) ~ EpiRR.uuid, value.var = count_col)
  
    # estimate size factors for count data
    sizeFactors <- DESeq2::estimateSizeFactorsForMatrix(counts_wide[, ..sample_names])
    
    # multiply counts with size factors
    counts_wide[, (sample_names):=Map("*", .SD, sizeFactors), .SDcols=sample_names]
    
    write_dt_and_md5(counts_wide, paste0(outfile_col, "_DESeq2.csv.gz"))
    
    ggsave(paste0(outfile_col, ".pdf"), pca_plot(counts_wide, outfile_col), width = 12, height = 10)
  }
  
}

write_sample_matrices(file_list = gene_files, id_col = "gene_id", outfile_prefix = "genes")

write_sample_matrices(file_list = iso_files, id_col = "transcript_id", outfile_prefix = "isoforms")
