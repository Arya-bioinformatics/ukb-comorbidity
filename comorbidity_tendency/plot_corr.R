# Title     : TODO
# Objective : TODO
# Created by: donggy
# Created on: 2020/5/16

library(corrplot)
ratio_data = read.csv('ratio.csv', row.names= 1)
q_data = read.csv('q.csv', row.names= 1)
ratio_data = data.matrix(ratio_data)
q_data = data.matrix(q_data)
# Insignificant correlations are leaved blank
#png('pathway.png', width=10, height=10,units="in",res=300, pointsize = 15)
pdf(file="a.pdf")
corrplot(ratio_data, method='circle', type="lower", order="original", p.mat = q_data, sig.level = 1,
         is.corr = FALSE, addgrid.col='lightgrey', insig='blank', cl.lim = c(0, 1), tl.col = "black", tl.srt = 45,
         hclust.method = "ward.D2")
dev.off()
