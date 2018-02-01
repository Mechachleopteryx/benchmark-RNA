



#!/usr/bin/env Rscript
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)< 3) {
  stop("At least 2 arguments must be supplied", call.=FALSE)
}


ggplot_all_metrics_err <- function()
{

 df <- read.table(args[1],h=T)
 df2 <- df[which(df$coverage == args[2]),]
 ggplot(data=df2, aes(x=as.factor(df2$error), y=df2$value, fill=df2$metric)) + geom_bar(stat="identity", position=position_dodge()) + coord_cartesian(ylim=c(0.75,1))  + guides(fill=guide_legend(title="Metric"))  + 
  xlab("Error rate") +
  ylab("Metric mean value") +
  ggtitle("Correction quality over different error rates in reads")

}

outName <- paste(args[3], "/metrics_function_errorrate.png", sep="")
ggsave(
  outName,
  ggplot_all_metrics_err(),
  width = 8,
  height = 8,
  dpi = 250
)


