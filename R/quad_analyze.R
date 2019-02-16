#!/usr/bin/Rscript
#  R/quad_analyze.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.30.2019

require(ggplot2)

load("data/sims/quadratic1550158396.01742.RData")

sets <- length(results)
titles <- c("Effort", "NLL", "InSample_Error", "OoS_Error", "Subspace_Error")
for (set in 1:sets) {
    pdf(paste("images/quadsim_", Ps[set], ".pdf", sep = ""), width = 10, height = 10)
    par(mfrow=c(3,2))
    for (r in 1:length(titles)) {
        if (r > 2) {
            boxplot(log(results[[set]][[r]]), main = titles[r])
        } else {
            boxplot(results[[set]][[r]], main = titles[r])
        }
    }
    dev.off()
}
