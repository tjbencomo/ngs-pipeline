library(readr)
library(dplyr)
library(ggplot2)

df <- read_csv(snakemake@input[[1]])
df <- df %>%
    mutate(type = factor(type, levels = c("normal", "tumor")))
p <- df %>%
    ggplot(aes(type, mean_depth)) +
    geom_boxplot(aes(fill = type)) +
    labs(x = "", y = "Average Depth") +
    guides(fill = guide_legend(title="")) +
    ggtitle("Sequencing Depths")
ggsave(snakemake@output[[1]], p)
