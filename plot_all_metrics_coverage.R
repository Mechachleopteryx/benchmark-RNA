



#!/usr/bin/env Rscript
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)< 2) {
  stop("At least 2 arguments must be supplied", call.=FALSE)
}


ggplot_all_metrics_cov <- function()
{

 df <- read.table(args[1],h=T)
 ggplot(data=df, aes(x=as.factor(df$coverage), y=df$value, fill=df$metric)) + geom_bar(stat="identity", position=position_dodge()) + coord_cartesian(ylim=c(0.75,1)) + guides(fill=guide_legend(title="Metric"))  + 
  xlab("Coverage for a gene") +
  ylab("Metric mean value") +
  ggtitle("Correction quality over different coverage values per gene")

}

outName <- paste(args[2], "/metrics_function_coverage.png", sep="")
ggsave(
  outName,
  ggplot_all_metrics_cov(),
  width = 8,
  height = 8,
  dpi = 250
)


