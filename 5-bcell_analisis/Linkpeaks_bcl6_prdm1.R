

bcl6<-df_links_DARs_DEA[df_links_DARs_DEA$gene=="BCL6",]


table(bcl6$gene=="BCL6",bcl6$cluster)

prdm1<-df_links_DARs_DEA[df_links_DARs_DEA$gene=="PRDM1",]

table(prdm1$gene=="PRDM1",prdm1$cluster)


bcl6_prdm1<- rbind(bcl6,prdm1)

tonsil_wnn_bcell_links_DARs_DEA_bcl6_prdm1<-tonsil_wnn_bcell_links

Links(tonsil_wnn_bcell_links_DARs_DEA_bcl6_prdm1)<-makeGRangesFromDataFrame(bcl6_prdm1,
keep.extra.columns=FALSE,
ignore.strand=FALSE,
seqinfo=NULL,
seqnames.field=c("seqnames", "seqname","chromosome",
"chrom","chr","chromosome_name","seqid"),
start.field="start",
end.field=c("end", "stop"),
strand.field="strand",
starts.in.df.are.0based=FALSE)

values(Links(tonsil_wnn_bcell_links_DARs_DEA_bcl6_prdm1))<-DataFrame(score=bcl6_prdm1$score,
gene=bcl6_prdm1$gene,
cluster=bcl6_prdm1$cluster,
peak=bcl6_prdm1$peak,
zscore=bcl6_prdm1$zscore,
pvalue=bcl6_prdm1$pvalue,
p_val.dars=bcl6_prdm1$p_val.dars,
pct.1.dars=bcl6_prdm1$pct.1.dars,
pct.2.dars=bcl6_prdm1$pct.2.dars,
avg_log2FC.dars=bcl6_prdm1$avg_log2FC.dars,
start.peak=bcl6_prdm1$start.peak,
end.peak=bcl6_prdm1$end.peak,
p_val.rna=bcl6_prdm1$p_val.rna,
avg_log2FC.rna=bcl6_prdm1$avg_log2FC.rna,
pct.1.rna=bcl6_prdm1$pct.1.rna,
pct.2.rna=bcl6_prdm1$pct.2.rna,
p_val_adj.rna=bcl6_prdm1$p_val_adj.rna)

table(bcl6_prdm1$gene=="BCL6",bcl6_prdm1$cluster)


coverage_extend <- function(x,y,seuratobject){purrr::map(y, function(y) {
p <- CoveragePlot(object = seuratobject, region = x, features = x, expression.assay = "RNA", idents = idents.plot, extend.upstream = y, extend.downstream = y, region.highlight = ranges.show
#tile = TRUE
)
p & scale_fill_manual(values = cols_cluster)
})}


ranges.show <- StringToGRanges("chr3-187721377-187745725")
cols_cluster <- c("#a6cee3", "#1f78b4","#b2df8a",
"#33a02c","#e31a1c")
coverage_extend("BCL6",c(1e+6,1e+7),tonsil_wnn_bcell_links_DARs_DEA_bcl6_prdm1)

coverage_extend("PRDM1",(5e+6,1e+7),tonsil_wnn_bcell_links_DARs_DEA_bcl6_prdm1)
coverage_extend("PRDM1",c(5e+6,1e+7),tonsil_wnn_bcell_links_DARs_DEA_bcl6_prdm1)
ranges.show <- StringToGRanges("chr6-105993463-106109939")
coverage_extend("PRDM1",c(5e+6,1e+7),tonsil_wnn_bcell_links_DARs_DEA_bcl6_prdm1)
ranges.show <- StringToGRanges("chr3-187721377-187745725")
coverage_extend("BCL6",c(0,2000,1e+4,1e+5,1e+6,1e+7),tonsil_wnn_bcell_links_DARs_DEA_bcl6_prdm1)
ranges.show <- StringToGRanges("chr3-187721377-187745725")




tonsil_wnn_bcell_links_DARs_DEA_bcl6_prdm1 <- LinkPeaks(
object = tonsil_wnn_bcell_links_DARs_DEA_bcl6_prdm1,
peak.assay = "ATAC",
expression.assay = "RNA",
expression.assay = "RNA",
genes.use = c("BCL6","PRDM1"),
distance = 1e+8
)



Idents(tonsil_wnn_bcell_links_DARs_DEA_bcl6_prdm1)
cols_cluster <- c("#a6cee3", "#1f78b4","#b2df8a",
"#33a02c","#e31a1c")
p <- CoveragePlot(
object = tonsil_wnn_bcell_links_DARs_DEA_bcl6_prdm1,
region = "BCL6",
features = "BCL6",
expression.assay = "RNA",
idents = idents.plot,
extend.upstream = 1e+4,
extend.downstream = 0
#region.highlight = ranges.show
#tile = TRUE
)
p <- CoveragePlot(
