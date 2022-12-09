library(UpSetR)

all = read.delim(snakemake@input[[1]], header=FALSE)
multiply = read.delim(snakemake@input[[2]], header=FALSE)
long = read.delim(snakemake@input[[3]], header=FALSE)

listInput1 = list("Short reads, 1+ exon-exon junctions" = all$V1, "Long reads" = long$V1)
listInput2 = list("Short reads, 2+ exon-exon junctions" = multiply$V1, "Long reads" = long$V1)

tiff(snakemake@output[[1]], units="in", width=5, height=5, res=300)
upset(fromList(listInput1))
dev.off()
tiff(snakemake@output[[2]], units="in", width=5, height=5, res=300)
upset(fromList(listInput2))
dev.off()
