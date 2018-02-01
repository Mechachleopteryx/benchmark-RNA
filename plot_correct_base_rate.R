
#!/usr/bin/env Rscript
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("At least 2 arguments must be supplied", call.=FALSE)
}


ggplot_recall <- function()
{
yy <- read.table(args[1])

#~ yy <- yy[order(yy$corrector),]
 ggplot(data=yy, aes(x=factor(yy$V1), y=yy$V2, fill=yy$V1))  + geom_boxplot()  + 
  xlab("Correction Methods") +
  ylab("Correct bases distribution in reads") +
  ggtitle("Distribution correct bases in reads for different correction methods")
}

outName <- paste(args[2], "/correct_base_rate.png", sep="")
ggsave(
  outName,
  ggplot_recall(),
  width = 8,
  height = 8,
  dpi = 250
)
