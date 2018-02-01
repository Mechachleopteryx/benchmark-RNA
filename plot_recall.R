
#!/usr/bin/env Rscript
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) < 2 ) {
  stop("At least 2 argument must be supplied", call.=FALSE)
}


ggplot_recall <- function()
{
yy <- read.table(args[1])

#~ yy <- yy[order(yy$corrector),]
 ggplot(data=yy, aes(x=factor(yy$V1), y=yy$V2, fill=yy$V1))  + geom_boxplot()  + guides(fill=guide_legend(title="Correction Method"))  + 
  xlab("Correction Methods") +
  ylab("Recall distribution in reads") +
  ggtitle("Distribution of recall for erroneous bases in reads for different correction methods")
}

outName <- paste(args[2], "/recall.png", sep="")

ggsave(
  outName,
  ggplot_recall(),
  width = 8,
  height = 8,
  dpi = 250
)
