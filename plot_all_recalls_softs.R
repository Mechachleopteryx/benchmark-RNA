



#!/usr/bin/env Rscript
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)< 2) {
  stop("At least 2 arguments must be supplied", call.=FALSE)
}

#soft precision coverage error
ggplot_softs_recall <- function()
{


yy <- read.table(args[1],h=T)

  ggplot(data=yy, aes(x=factor(yy$soft), y=yy$value, fill=yy$soft))  + geom_boxplot()  + 
  xlab("Correction Methods") +
  ylab("Precision distribution in reads") +
  ggtitle("Distribution of precision for corrected bases in reads for different correction methods")  + facet_grid(yy$coverage~yy$error )   + theme(panel.grid.minor = element_blank())

}

outName <- paste(args[2], "/all_recall_softs.png", sep="")
ggsave(
  outName,
  ggplot_softs_recall(),
  width = 8,
  height = 8,
  dpi = 250
)




