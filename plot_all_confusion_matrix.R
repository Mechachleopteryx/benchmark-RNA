#!/usr/bin/env Rscript
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)< 4) {
  stop("At least 4 arguments must be supplied", call.=FALSE)
}

#file format
#reference correction ratio soft
#exclusion exclusion 1.0 msa_exon
#inclusion inclusion 1.0 msa_exon

ggplot_matrixconfAll <- function()
{
yy <- read.table(args[1], h=T)
df1 <- data.frame(yy[which(yy$coverage == args[3]),])
df <- df1[which(df1$error == args[4]),]

label_names <- c(
                    `10` = "Skipped exon size:10",
                    `50` = "Skipped exon size:50",
                    `100` = "Skipped exon size:100",
                    `LoRMA` = "LoRMA",
                    `LoRDEC` = "LoRDEC",
                    `Proovread` = "Proovread",
                    `MECAT` = "MECAT",
                    `PBDagCon` = "PBDagCon",
                    `daccord` = "daccord",
                    `msa_exon` = "msa_exon",
                    `msa_isoform` = "msa_isoform",
                    `colorMap` = "colorMap"
                    )

ggplot(df, aes(df$reference, df$correction)) + geom_tile(aes(fill=df$ratio)) +  geom_text(aes(label = round(df$ratio, 2))) +  scale_fill_gradient(low = "white", high = "red") + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none") + ggtitle("Ratios of isoform counts in reads before and after correction") + xlab("Reference isoforms") + ylab("Corrected isoforms")  + facet_grid(~df$soft,labeller = as_labeller(label_names) )   + theme(panel.grid.minor = element_blank())

}



outName <- paste(args[2], "/all_confusion_matrix.png", sep="")


ggsave(
  outName,
  ggplot_matrixconfAll(),
  width = 8,
  height = 8,
  dpi = 250
)
