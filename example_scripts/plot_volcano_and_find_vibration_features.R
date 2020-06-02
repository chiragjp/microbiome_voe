
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

args <- commandArgs(trailingOnly = TRUE)

data=readRDS(args[[1]])

plot=ggplot(data, aes(x = estimate, y = -log10(by.p.val), label = feature)) +ggtitle('Example CRC vs. species abundances meta-analytic results')+ geom_point(size=1) +theme_cowplot(12)+ geom_hline(yintercept = -log10(0.05)) + xlim(-2, 2)+ylab('-log10(p.value)')+xlab('Effect Size')
ggsave('./example_pipeline_output/volcano_plot_crc_example.png')

to_vibrate = data %>% filter(by.p.val < 0.05) %>% select(feature) %>% unname %>% unlist
saveRDS(to_vibrate,args[[2]])
