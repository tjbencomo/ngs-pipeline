Sys.setenv("TZDIR"=paste0(Sys.getenv("CONDA_PREFIX"), "/share/zoneinfo"))

library(readr)
library(dplyr)
library(ggplot2)

df <- read_csv(snakemake@input[[1]])
df <- df %>%
    mutate(sample_type = factor(sample_type, levels = c("normal", "tumor")))
p <- df %>%
    ggplot(aes(sample_type, mean)) +
    geom_boxplot(aes(fill = sample_type)) +
    labs(x = "", y = "Average Depth") +
    guides(fill = guide_legend(title="")) +
    ggtitle("Sequencing Depths")
ggsave(snakemake@output[[1]], p)
