
#!/usr/bin/env Rscript
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)< 2) {
  stop("At least 2 arguments must be supplied", call.=FALSE)
}


ggplot_size <- function()
{

 x <- read.table(args[1],h=T)
 ggplot(data=x, aes(x=factor(x$soft), y=x$size, fill=x$soft))  + geom_boxplot() + guides(fill=guide_legend(title="Correction Method"))  + 
  xlab("Correction Methods") +
  ylab("Ratio (corrected size/real size)") +
  ggtitle("Distribution of read size ratios for different correction methods")

}

outName <- paste(args[2], "/size.png", sep="")
ggsave(
  outName,
  ggplot_size(),
  width = 8,
  height = 8,
  dpi = 250
)
