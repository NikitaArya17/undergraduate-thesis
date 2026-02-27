library(ggplot2)
library(tidyr)
library(dplyr)

input_comparison_file <- "3_ML/output_data/Ensemble_model_comparison.csv"
output_plot_file <- "3_ML/output_data/Ensemble_comparison.png"

df_wide <- read.csv(input_comparison_file, stringsAsFactors = FALSE)

df_long <- df_wide %>%
  pivot_longer(
    cols = -Model,
    names_to = "Metric",
    values_to = "Value"
  )

df_long$Metric <- gsub("_", " ", df_long$Metric)

model_order <- c("SVM", "NB", "ANN", "Ensemble")
df_long <- df_long %>%
  filter(Model %in% model_order) %>%
  mutate(Model = factor(Model, levels = model_order))

p_bar <- ggplot(df_long, aes(x = Metric, y = Value, fill = Model)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.7) +

  scale_fill_manual(values = c(
    "SVM" = "#3498db",
    "NB"  = "#2ecc71",
    "ANN" = "#9b59b6",
    "Ensemble" = "#e74c3c"
  )) +

  scale_y_continuous(limits = c(0, 1.05),
                     breaks = seq(0, 1, 0.2),
                     expand = c(0, 0)) +
  labs(
    title = "E. coli Model Performance Comparison",
    y = "Score Value",
    x = "", 
    fill = "Model Type"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(color = "gray20"),
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +

  geom_text(aes(label = sprintf("%.2f", Value)), 
            position = position_dodge(width = 0.8), 
            vjust = -0.5, size = 3)

print(p_bar)

ggsave(output_plot_file, plot = p_bar, width = 10, height = 6, device = "png")
cat("Saved")