library(ggplot2)
library(gridExtra)
library(reshape2)
library(grid)

args <- commandArgs(trailingOnly = TRUE)

dir <- getwd()
setwd(dir)
#Bar plot of mean and standard deviatio for each sample
dat1 <- read.delim(args[1], header=T)

sink('depth_summary_report.md', append = T)
library(knitr)
kable((dat1[c(1,3:5)]), type = 'pandoc', digits=3, align='c')
cat('\n')
cat('\n')
sink()


limits <- aes(ymax = MEAN_COVERAGE + SD_COVERAGE, ymin = MEAN_COVERAGE - SD_COVERAGE)
depth_plot <- ggplot(dat1, aes(x=reorder(Sample, MEAN_COVERAGE, min), y = MEAN_COVERAGE)) + geom_bar(stat = "identity", fill="steelblue")
depth_plot <- depth_plot + geom_errorbar(limits, width=0.25) + theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1))
depth_plot <- depth_plot + labs(x = "Sample", y = "Mean depth +/- s.d.")
ggsave("mean_depth.pdf", h=4, w=6)

#Cumulative coverage plot
dat1 <- read.delim("depth_summary.txt", header=T)
precent_cov <- dat1[c(1, 14:26)] 
dat2 <- melt(precent_cov, id.vars="Sample")
dat2$cov <- factor(dat2$variable, labels=c("5x", "10x", "15x", "20x", 
                                           "25x", "30x", "40x", "50x", "60x", 
                                           "70x", "80x", "90x", "100x" ))
depth_percent_cov <- ggplot(dat2, aes(x=cov, y=value)) + geom_bar(stat="identity", fill="steelblue")
depth_percent_cov <- depth_percent_cov + facet_wrap(~ Sample, ncol=2 ) + theme_bw(base_size = 10) + theme(axis.text.x=element_text(angle=45, hjust = 1))
depth_percent_cov <- depth_percent_cov + labs(y="Fraction of genome covered", x="Coverage" ) + geom_text(aes(label=round(value, 2)), vjust=0, size=4)
ggsave("fraction_covered.pdf", h=10, w=8)


# Histogram of coverage
dat3 <- read.delim(args[2], header=T)
depth_hist <- ggplot(dat3, aes(x=factor(coverage), y = count)) + geom_bar(stat="identity", fill="steelblue")
depth_hist <- depth_hist + facet_wrap( ~ Sample, ncol = 2) + coord_cartesian(xlim=c(0, 140)) + labs(x="Coverage")
depth_hist <- depth_hist + theme_bw() + scale_x_discrete(breaks=c("0","20","40","60","80", "100", "120", "140"))
ggsave("depth_histograms.pdf", h=10, w=8)

# write the WGSmetrics results as a table

#pdf('wgsmetric_table.pdf', w=19.5, h=8.75)
#grid.table(dat1[c(1, 3:10)], theme=ttheme_default())
#dev.off()





