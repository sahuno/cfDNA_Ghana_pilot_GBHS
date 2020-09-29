# Code snippet to get chromosomal position and cytoband from gene symbols using
# biomaRt
# Kenneth Daily, 2014

library(biomaRt)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",version="Ensembl Genes 100")

datasets = listDatasets(ensembl)

# Only use standard human chromosomes
normal.chroms <- c(1:22, "X", "Y", "M")

# Filter on HGNC symbol and chromosome, retrieve genomic location and band
my.symbols <- c("RB1", "TP53", "AKT3")

my.regions <- getBM(c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "band"),
                    filters = c("hgnc_symbol", "chromosome_name"),
                    values = list(hgnc_symbol=my.symbols, chromosome_name=normal.chroms),
                    mart = ensembl)


library(biomaRt)
####use hg19 ensemsbl dataset to annotate   ###############################
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
listDatasets(grch37)
my.regions.hg19 <- getBM(c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "band","gene_biotype"),
                    filters = c("hgnc_symbol", "chromosome_name"),
                    values = list(hgnc_symbol=my.symbols, chromosome_name=normal.chroms),
                    mart = grch37)



top.cna.genes2.cyto <- getBM(c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "band","gene_biotype"),
      filters = c("hgnc_symbol", "chromosome_name"),
      values = list(hgnc_symbol=top.cna.genes2, chromosome_name=normal.chroms),
      mart = grch37)

anal_cancer_specific_genes_list <- getBM(c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "band","gene_biotype"),
                                         filters = c("hgnc_symbol", "chromosome_name"),
                                         values = list(hgnc_symbol=anal_cancer_specific_genes, chromosome_name=normal.chroms),
                                         mart = grch37)


###get list of all atributes of all grch37, which i could  filter with
attributes.grch37 = listAttributes(grch37)
View(head(attributes.grch37,100))

##get all esnseml gene list hg19
ensembl.hg19 <- getBM(c("chromosome_name", "start_position", "end_position","strand","hgnc_symbol","band","ensembl_gene_id","gene_biotype"),
      filters = c("chromosome_name"),
      values = list(chromosome_name=normal.chroms),
      mart = grch37)

class(ensembl.hg19) ##data frame
###correct strand information for 
ensembl.hg19.strandChnaged <- ensembl.hg19
ensembl.hg19.strandChnaged$strand <- gsub("-1","-",ensembl.hg19.strandChnaged$strand)
ensembl.hg19.strandChnaged$strand <- gsub("1","+",ensembl.hg19.strandChnaged$strand)

ensembl.hg19.strandChnaged <- makeGRangesFromDataFrame(ensembl.hg19.strandChnaged,start.field= "start_position", end.field= "end_position",keep.extra.columns = T)

str(ensembl.hg19.strandChnaged)

