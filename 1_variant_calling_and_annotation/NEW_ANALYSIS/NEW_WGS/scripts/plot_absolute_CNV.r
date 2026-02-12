library(dplyr)
library(ggplot2)

ratio_file  <- "NEW_ANALYSIS/NEW_WGS/cnv/recal_BSQR.bam_ratio.txt"
pvalue_file <- "NEW_ANALYSIS/NEW_WGS/cnv/recal_BSQR.bam_ratio.txt.p.value.txt"
output_pdf  <- "NEW_ANALYSIS/NEW_WGS/cnv/recal_BSQR_absolute_CNV.pdf"

ratio_df <- read.table(ratio_file,
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors = FALSE)

cnv_segments <- read.table(pvalue_file,
                           header = TRUE,
                           sep = "\t",
                           stringsAsFactors = FALSE,
                           check.names = FALSE)

ratio_df$Start  <- as.numeric(ratio_df$Start)
ratio_df$CopyNumber <- as.numeric(ratio_df$CopyNumber)
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
  mutate(Max_Pos = max(Start, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(Plot_X = (Chr_ID * 100) + (Start / Max_Pos * 90))

background_df <- gapminder %>%
  filter(gain.loss == "neutral") %>%
  sample_frac(0.1)
highlight_df  <- gapminder %>% filter(gain.loss != "neutral", pvalue <= 0.05)
plot_data     <- bind_rows(background_df, highlight_df)

x_breaks <- seq(100, length(unique_chrs) * 100, 100)
x_labels <- seq_len(length(unique_chrs))
x_separators <- seq(195, length(unique_chrs) * 100, 100)

pdf(output_pdf, width = 20, height = 6)
theme_set(theme_bw(base_size = 18))

p <- ggplot(plot_data, aes(x = Plot_X, y = CopyNumber)) +
  geom_point(data = background_df, color = "darkgrey", size = 0.5) +
  geom_point(data = highlight_df, aes(color = gain.loss), size = 1.0) +
  scale_x_continuous(name = "Chromosome",
                     breaks = x_breaks,
                     labels = x_labels) +
  scale_y_continuous(name = "Absolute Copy Number",
                     breaks = seq(0, 20, 2),
                     limits = c(0, 14)) +

  scale_color_manual(values = c("gain" = "red", "loss" = "blue")) +

  geom_vline(xintercept = x_separators, linetype = "dashed", color = "grey80") +

  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

print(p)
dev.off()
message("Success! Huzzah!")