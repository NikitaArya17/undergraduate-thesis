library(dplyr)

ratio_df <- read.table("NEW_ANALYSIS/NEW_WGS/cnv/recal_BSQR.bam_ratio.txt",
                       header = TRUE,
                       stringsAsFactors = FALSE)
cnv_segments <- read.table("NEW_ANALYSIS/NEW_WGS/cnv/recal_BSQR.bam_ratio.txt.p.value.txt",
                           header = TRUE,
                           sep = "\t",
                           stringsAsFactors = FALSE)

ratio_df$pvalue <- 1.0
ratio_df$gain.loss <- "neutral"

for (i in seq_len(row(cnv_segments))) {
  hits <- which(ratio_df$Chromosome == cnv_segments$chr[i] &
                  ratio_df$Start >= cnv_segments$start[i] &
                  ratio_df$Start <= cnv_segments$end[i])

  if (length(hits) > 0) {
    ratio_df$pvalue[hits] <- cnv_segments$WilcoxonRankSumTestPvalue[i]
    ratio_df$gain.loss[hits] <- cnv_segments$status[i]
  }
}

ratio_df <- ratio_df %>%
  rename(location = Start,
         Ratio_2 = Ratio)

write.table(ratio_df,
            "NEW_ANALYSIS/NEW_WGS/cnv/recal_BSQR_BAM_pvalue_ratio_final.txt",
            row.names = FALSE)