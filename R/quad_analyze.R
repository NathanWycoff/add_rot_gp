#!/usr/bin/Rscript
#  R/quad_analyze.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.30.2019

require(ggplot2)

load("data/sims/quadratic1548914553.61066.RData")

sets <- length(results)
titles <- c("Effort", "NLL", "InSample_Error", "OoS_Error", "Subspace_Error")
results[[1]][[1]]
for (set in 1:sets) {
    for (r in 1:length(titles)) {
        pdf(paste("images/quadsim_", Ps[set], "_", titles[r], ".pdf", sep = ""))
        if (r > 2) {
            boxplot(log(results[[set]][[r]]), main = titles[r])
        } else {
            boxplot(results[[set]][[r]], main = titles[r])
        }
        dev.off()
    }
}
