#!/usr/bin/env Rscript
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)< 2) {
  stop("At least 3 arguments must be supplied", call.=FALSE)
}


ggplot_matrixconf <- function()
{
yy <- read.table(args[1], h=T)
df <- data.frame(yy)

ggplot(df, aes(df$V1, df$V2)) + geom_tile(aes(fill=df$V3)) +  geom_text(aes(label = round(df$V3, 2))) +  scale_fill_gradient(low = "white", high = "red") + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none") + ggtitle("HMP gut (N=100, Rs=40M)") + xlab("Truth distances") + ylab("Minhash distances")

ggplot(df, aes(df$reference, df$correction)) + geom_tile(aes(fill=df$ratio)) +  geom_text(aes(label = round(df$ratio, 2))) +  scale_fill_gradient(low = "white", high = "red") + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none") + ggtitle("Ratios of isoform counts in reads before and after correction") + xlab("Reference isoforms") + ylab("Corrected isoforms")
}

outName <- paste(args[3], "/confusion_matrix_", args[2],".png", sep="")
ggsave(
  outName,
  ggplot_matrixconf(),
  width = 8,
  height = 8,
  dpi = 250
)

