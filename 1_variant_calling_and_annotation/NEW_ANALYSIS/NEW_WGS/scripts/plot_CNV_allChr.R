library(dplyr)
library(ggplot2)

ratio_file  <- "NEW_ANALYSIS/NEW_WGS/cnv/recal_BSQR.bam_ratio.txt"
pvalue_file <- "NEW_ANALYSIS/NEW_WGS/cnv/recal_BSQR.bam_ratio.txt.p.value.txt"
output_pdf  <- "NEW_ANALYSIS/NEW_WGS/cnv/recal_BSQR_CNV_plot.pdf"

ratio_df     <- read.table(ratio_file,
                           header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE)

cnv_segments <- read.table(pvalue_file,
                           header = TRUE,
                           sep = "\t",
                           stringsAsFactors = FALSE,
                           check.names = FALSE)

if(!("Ratio" %in% colnames(ratio_df))) {
  stop("Error: The 'ratio_file' seems wrong. It should contain a 'Ratio' column.")
}
if(!("status" %in% colnames(cnv_segments))) {
  stop("Error: The 'pvalue_file' seems wrong. It should contain a 'status' column.")
}

ratio_df$Start <- as.numeric(ratio_df$Start)
ratio_df$Ratio <- as.numeric(ratio_df$Ratio)
ratio_df$pvalue <- 1.0
ratio_df$gain.loss <- "neutral"

for(i in seq_len(nrow(cnv_segments))) {
  hits <- which(ratio_df$Chromosome == cnv_segments$chr[i] &
                  ratio_df$Start >= cnv_segments$start[i] &
                  ratio_df$Start <= cnv_segments$end[i])

  if(length(hits) > 0) {
    ratio_df$pvalue[hits]    <- cnv_segments$WilcoxonRankSumTestPvalue[i]
    ratio_df$gain.loss[hits] <- cnv_segments$status[i]
  }
}


unique_chrs <- unique(ratio_df$Chromosome)
chr_map <- data.frame(Chromosome = unique_chrs,
                      Chr_ID = seq_len(length(unique_chrs)))

gapminder <- ratio_df %>%
  left_join(chr_map, by = "Chromosome") %>%
  group_by(Chromosome) %>%
  mutate(Max_Pos = max(Start, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    Plot_X = (Chr_ID * 100) + (Start / Max_Pos * 90)
  ) %>%
  rename(Ratio_2 = Ratio)

background_df <- gapminder %>%
  filter(gain.loss == "neutral") %>%
  sample_frac(0.1)

highlight_df <- gapminder %>%
  filter(gain.loss != "neutral", pvalue <= 0.05)

plot_data <- bind_rows(background_df, highlight_df)

pdf(output_pdf, width = 20, height = 6)
theme_set(theme_bw(base_size = 18))

p <- ggplot(plot_data, aes(x = Plot_X, y = Ratio_2)) +
  geom_point(data = background_df, color = "darkgrey", size = 0.8) +
  geom_point(data = highlight_df, aes(color = gain.loss), size = 1.2) +
  scale_x_continuous(name = "Chromosome",
                     breaks = seq(100, length(unique_chrs)*100, 100),
                     labels = seq_len(length(unique_chrs))) +
  scale_y_continuous(name = "Ratio", limits = c(0, 8)) +
  scale_color_manual(values = c("gain" = "red", "loss" = "blue")) +
  geom_vline(xintercept = seq(195, length(unique_chrs) * 100, 100),
             linetype = "dashed", color = "grey80") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

print(p)
dev.off()