# Title     : TODO
# Objective : TODO
# Created by: donggy
# Created on: 2020/5/16

args=commandArgs(T)
parameter1 = args[1]

library(corrplot)
ratio_data = read.csv('ratio.csv', row.names= 1)
q_data = read.csv('q.csv', row.names= 1)
ratio_data = data.matrix(ratio_data)
q_data = data.matrix(q_data)
# Insignificant correlations are leaved blank
pdf(parameter1, width=10, height=10)
corrplot(ratio_data, method='pie', type="lower", order="original", p.mat = q_data, sig.level = 0.05,
         is.corr = FALSE, addgrid.col='lightgrey', insig='blank', tl.col = "black", tl.srt = 45,
         hclust.method = "ward.D2")
dev.off()
