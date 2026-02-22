# ---- GWAS fine-mapping SNP priors ----

#' Read fine-mapping SNP table
#' @param path tsv/csv path
#' @param chr_col,pos_col,pip_col column names
#' @param sep separator; NULL auto-detects by extension (.csv => ',', else tab)
#' @return data.frame
read_finemap_snps <- function(path, chr_col = "chr", pos_col = "pos", pip_col = "pip", sep = NULL) {
  if (is.null(sep)) {
    sep <- if (grepl("\\.csv$", tolower(path))) "," else "\t"
  }
  df <- utils::read.table(path, header = TRUE, sep = sep, stringsAsFactors = FALSE, check.names = FALSE)
  if (!(chr_col %in% names(df))) stop("chr_col not found")
  if (!(pos_col %in% names(df))) stop("pos_col not found")
  if (!(pip_col %in% names(df))) stop("pip_col not found")
  out <- df[, c(chr_col, pos_col, pip_col), drop=FALSE]
  colnames(out) <- c("chr", "pos", "pip")
  out
}

#' Convert SNP table to GRanges
#' @param snp_df data.frame with chr,pos,pip
#' @param genome genome string stored in GRanges
#' @return GRanges
snps_to_granges <- function(snp_df, genome = NA_character_) {
  .check_cols(snp_df, c("chr", "pos", "pip"))
  chr <- as.character(snp_df$chr)
  # normalize chr prefix: keep as-is if already has 'chr'
  chr <- ifelse(grepl("^chr", chr), chr, paste0("chr", chr))
  pos <- as.integer(snp_df$pos)
  gr <- GenomicRanges::GRanges(seqnames = chr, ranges = IRanges::IRanges(start = pos, end = pos))
  GenomicRanges::mcols(gr)$pip <- as.numeric(snp_df$pip)
  GenomeInfoDb::genome(gr) <- genome
  gr
}

#' Map SNP GRanges to ATAC peaks by overlap
#' @param obj Seurat object
#' @param snps_gr GRanges with mcols$pip
#' @param atac_assay ATAC assay name
#' @return data.frame peak,pip
map_snps_to_peaks <- function(obj, snps_gr, atac_assay = "ATAC") {
  peak_gr <- Signac::granges(obj[[atac_assay]])
  hits <- GenomicRanges::findOverlaps(snps_gr, peak_gr, ignore.strand = TRUE)
  if (length(hits) == 0) {
    return(data.frame(peak = character(), pip = numeric()))
  }
  peak_id <- as.character(GenomicRanges::queryHits(hits))
  peak_nm <- names(peak_gr)[GenomicRanges::subjectHits(hits)]
  pip <- GenomicRanges::mcols(snps_gr)$pip[GenomicRanges::queryHits(hits)]
  data.frame(peak = peak_nm, pip = pip)
}

#' Aggregate SNP PIP into peak-level weights
#' @param snp_peak_df data.frame peak,pip
#' @param method max or sum
#' @return named numeric vector weights per peak
peak_weights_from_snps <- function(snp_peak_df, method = c("max", "sum"), pip_floor = 0) {
  method <- match.arg(method)
  if (nrow(snp_peak_df) == 0) return(setNames(numeric(), character()))
  snp_peak_df <- snp_peak_df[snp_peak_df$pip >= pip_floor, , drop=FALSE]
  if (nrow(snp_peak_df) == 0) return(setNames(numeric(), character()))
  if (method == "max") {
    w <- tapply(snp_peak_df$pip, snp_peak_df$peak, max)
  } else {
    w <- tapply(snp_peak_df$pip, snp_peak_df$peak, sum)
  }
  as.numeric(w) |> setNames(names(w))
}

#' Convert peak weights to penalty factors for glmnet
#' Smaller penalty => less shrinkage => more likely selected.
#'
#' @param peak_weights named numeric (0..1 typical)
#' @param peaks character vector of peaks to align
#' @param alpha strength parameter
#' @param min_penalty lower bound
#' @return numeric penalty.factor aligned to peaks
penalty_from_peak_weights <- function(peak_weights, peaks, alpha = 2, min_penalty = 0.2) {
  w <- rep(0, length(peaks)); names(w) <- peaks
  if (length(peak_weights) > 0) {
    idx <- intersect(names(peak_weights), peaks)
    w[idx] <- peak_weights[idx]
  }
  pen <- exp(-alpha * w)
  pen <- pmax(min_penalty, pen)
  unname(pen)
}
