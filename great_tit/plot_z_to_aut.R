library(ggplot2)

#setwd('/home/bo1pgc/parus_reseq/scripts')

args<-commandArgs(TRUE)
dat <- read.delim(args[1], header=F)

p <- ggplot(dat, aes(x = reorder(V1, V4, min), y = V4)) 
p <- p + geom_bar(stat = "identity", fill="Steel blue") + geom_hline(yintercept=0.5, linetype=2)
p <- p + theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1)) 
p <- p +  labs(x="Sample", y="Z Coverage / Autosome Coverage ")


ggsave(args[2], h=4, w=6)

