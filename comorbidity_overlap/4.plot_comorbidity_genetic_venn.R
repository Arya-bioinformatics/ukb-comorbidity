library (VennDiagram)
library(grDevices)


snp_overlap = list()
c = file('../overlap/comorbidity_snp.txt',"r")
line=readLines(c,n=1)
line=readLines(c,n=1)
while( length(line) != 0 ) {
  x = strsplit(line, '\t')
  snp_overlap = c(snp_overlap, paste(x[[1]][1], x[[1]][2], sep = '*'))
  line=readLines(c,n=1)
}
close(c)


gene_overlap = list()
c = file('../overlap/comorbidity_gene.txt',"r")
line=readLines(c,n=1)
line=readLines(c,n=1)
while( length(line) != 0 ) {
  x = strsplit(line, '\t')
  gene_overlap = c(gene_overlap, paste(x[[1]][1], x[[1]][2], sep = '*'))
  line=readLines(c,n=1)
}
close(c)


ppi_overlap = list()
c = file('../overlap/comorbidity_ppi.txt',"r")
line=readLines(c,n=1)
line=readLines(c,n=1)
while( length(line) != 0 ) {
  x = strsplit(line, '\t')
  ppi_overlap = c(ppi_overlap, paste(x[[1]][1], x[[1]][2], sep = '*'))
  line=readLines(c,n=1)
}
close(c)



pathway_overlap = list()
c = file('../overlap/comorbidity_pathway.txt',"r")
line=readLines(c,n=1)
line=readLines(c,n=1)
while( length(line) != 0 ) {
  x = strsplit(line, '\t')
  pathway_overlap = c(pathway_overlap, paste(x[[1]][1], x[[1]][2], sep = '*'))
  line=readLines(c,n=1)
}
close(c)

rg_overlap = list()
c = file('../overlap/comorbidity_rg.txt',"r")
line=readLines(c,n=1)
line=readLines(c,n=1)
while( length(line) != 0 ) {
  x = strsplit(line, '\t')
  rg_overlap = c(rg_overlap, paste(x[[1]][1], x[[1]][2], sep = '*'))
  line=readLines(c,n=1)
}
close(c)



# plot venne
pdf(file="venn.pdf")
venn = venn.diagram(x =list(SNP= snp_overlap, Gene = gene_overlap, 'Genetic correlation' = rg_overlap,
                     PPI = ppi_overlap, Pathway = pathway_overlap), filename = NULL,
             height = 1500, width= 1500, resolution =500,imagetype="png", lwd=0.6, 
             fill =c("cornflowerblue","green","yellow", 'chocolate', 'darkorchid1'),cex=1, cat.cex=1,
             cat.just=list(c(0.5,1.1), c(-0.7,-3), c(0.5,0.1), c(1,1), c(1.3,-2.8)),
             cat.fontfamily='Helvetica', fontfamily='Helvetica', alpha=c(0.8,0.8,0.8,0.8,0.8))
grid.draw(venn)
dev.off()