#!/usr/bin/env Rscript
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
}



ggplot_hybrid <- function()
{
df1 <- read.table(args[1], h=T)
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
                    `MSA` = "MSA",
                    `colorMap` = "colorMap"
                    )
                    
               
df_hybrid <- df1[which(df1$corrector == "LoRDEC" | df1$corrector == "Proovread" | df1$corrector == "colorMap"),]
ggplot(df_hybrid,aes(x=factor(ratio),y=percent, fill=factor(ratio)), color=factor(ratio))   + geom_bar(stat="identity") + facet_grid(corrector~size_skipped_exon,labeller = as_labeller(label_names) )  + xlab("% of inclusion (major) isoform in expressed RNA")+ ylab("% of reads corrected to inclusion (major) isoform") +scale_fill_brewer(palette = "Paired")    + theme(legend.position='none', text = element_text(size=20)) + scale_y_continuous(breaks = c(0,50,75,90)) + theme(panel.grid.minor = element_blank())

}




ggsave(
  "corrected_to_inclusionmajor_hybrid.png",
  ggplot_hybrid(),
  width = 8,
  height = 8,
  dpi = 250
)
