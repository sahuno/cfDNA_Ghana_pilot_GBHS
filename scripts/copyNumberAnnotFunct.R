###date August 15th 2020
### function to annotate cna seg files with hg19 genes

#### load libraries
library(biomaRt)
library(GenomeInfoDb)
library(GenomicFeatures)
library(data.table)
library(data.table)
library(tidyverse)
library(reshape2)
library(plyr)

####use hg19 ensemsbl dataset to annotate   ###############################
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl") ### use ensemble legacy database for hg19. current is hg38
ensembl.hg19 <- getBM(c("chromosome_name", "start_position", "end_position","strand","hgnc_symbol","band","ensembl_gene_id","gene_biotype"),
                      filters = c("chromosome_name"),
                      values = list(chromosome_name=c(1:22, "X", "Y", "M")),
                      mart = grch37)

###correct strand information for easy conversion Granges
ensembl.hg19.strandChnaged <- ensembl.hg19
ensembl.hg19.strandChnaged$strand <- gsub("-1","-",ensembl.hg19.strandChnaged$strand)
ensembl.hg19.strandChnaged$strand <- gsub("1","+",ensembl.hg19.strandChnaged$strand)

##final Granges Object of ensembl genes hg19
ensembl.hg19.strandChnaged <- makeGRangesFromDataFrame(ensembl.hg19.strandChnaged,start.field= "start_position", end.field= "end_position",keep.extra.columns = T)


##read a single cna.seg file
#file.path = "/Users/samuelahuno/Polak_lab10082019/Ghana_cfDNA/IchorCNA_vr2_RepTime_Results/results_ichor_minSegBin20_AltFractThresh_0.01/ichorCNA/GBR00332_1000kb/GBR00332_1000kb.cna.seg"

########### read a files in a driectory
cna.paths <- "/Users/samuelahuno/Polak_lab10082019/Ghana_cfDNA/IchorCNA_vr2_RepTime_Results/IchorCNA_vr2_RepTime_Results.Backup/results_ichor_minSegBin20_AltFractThresh_0.01/ichorCNA"

cna.paths <- "/Users/samuelahuno/Polak_lab10082019/Ghana_cfDNA/IchorCNA_vr2_RepTime_Results/IchorCNA_vr2_RepTime_Results.Backup/results_ichor_minSegBin20_AltFractThresh_0.01/ichorCNA"
list.cna.files <- list.files(path=cna.paths,full.names=TRUE,pattern = "cna.seg", recursive = TRUE)
str(list.cna.files)

########################################## Begining of function to read .seg file and annotate with ensemble genes ##################
AnnotateEnsembleGenes <- function(geneList=ensembl.hg19.strandChnaged,filePath=""){
  df.seg <- read.delim(filePath, header=TRUE, sep = "\t", strip.white = TRUE)
  dr.seg <- makeGRangesFromDataFrame(df.seg,start.field=  "start", end.field= "end", seqnames.field="chr", keep.extra.columns = T)
  sampleName <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(filePath)))
  ##find overlaps between cna.seg and ensemble genes
  Ov.refGenes.in.CNA.dr.seg <- findOverlaps(query=ensembl.hg19.strandChnaged, subject=dr.seg, type = "within") ##complete overlaps
  unique(queryHits(Ov.refGenes.in.CNA.dr.seg))
  queryLength(Ov.refGenes.in.CNA.dr.seg)  ##sanity checks
  
  ###assigning copy number values to gene hits ##same as previous line
  ensembl.hg19.strandChnaged$copy_number <- NA
  ##Get Subject hits as dataFrame
  sub_hit <- as.data.frame(dr.seg[subjectHits(Ov.refGenes.in.CNA.dr.seg)])
  ensembl.hg19.strandChnaged[queryHits(Ov.refGenes.in.CNA.dr.seg)]$copy_number <- sub_hit[,grepl("Corrected_Copy_Number",names(sub_hit))]
  annotated_Seg_EnsembleGenes <- ensembl.hg19.strandChnaged
  
  #print(tools::file_path_sans_ext(tools::file_path_sans_ext(basename(filePath))))
  #print(annotated_Seg_EnsembleGenes[seqnames(annotated_Seg_EnsembleGenes)=="10", ][85:95])
  #print(annotated_Seg_EnsembleGenes[which(annotated_Seg_EnsembleGenes$Gene.name=="PTEN")])
  #print(annotated_Seg_EnsembleGenes[which(annotated_Seg_EnsembleGenes$Gene.name=="CCND1")])
  write.table(annotated_Seg_EnsembleGenes, file=paste0(sampleName,".annotated.","ensemblGenes.tsv"), quote=FALSE, sep='\t', col.names = NA)
  
  ##### assign returns of function to variables for easy processing
  df.metadata <- as.data.frame(elementMetadata(annotated_Seg_EnsembleGenes))
  df.sample <- as.data.frame(df.metadata[,grepl("copy_number",names(df.metadata))])
  geneNames <- as.data.frame(df.metadata[,grepl("hgnc_symbol",names(df.metadata))])
  cytobands <- as.data.frame(df.metadata[,grepl("band",names(df.metadata))])
  

  df.output <- data.frame(geneNames,df.sample,cytobands)
  names(df.output)[1] <- "Gene_Names" ##rename col. 1  with gene names
  names(df.output)[2] <- sampleName    ###rename col 2. with filenames\sample names 
  names(df.output)[3] <- "cytoband"    #### make availble cytoband infor just in case
  return(df.output)
}
############################################ end of annotation function ####################################

###TEST the function of 
#path_IchorP04_500kb <-"/Users/samuelahuno/Polak_lab10082019/Chimera_local_10092019/GavinHaLab/snakemake_GHBS30x_4MissingGenes/results/ichorCNA/P04/P04.cna.seg"
#AnnotateEnsembleGenes(geneList=gr.ensemblGenes.biomart.hg19,filePath=path_IchorP04_500kb)

### create empty dataFrame to hold copy numbers
data1 <- data.frame(stringsAsFactors=FALSE)

### for loop to run all samples
for (dat in list.cna.files) {
    df.data <- AnnotateEnsembleGenes(geneList=gr.ensemblGenes.biomart.hg19,filePath=dat)
    #data1$v1 <- NA
    #data1$v1 <- df.data[1]
    data1 <- as.data.frame(c(data1,df.data))
 }

##get gene names
getGeneName<- data1[,grepl("Gene",names(data1))][1]
getCytobands <- data1[,grepl("cytoband",names(data1))][1]

removeMatch <- c("Gene","cytoband") ### patterns to remove
##get cna nunbers
df.annotations.genes <- cbind(getGeneName,getCytobands,data1[,!grepl(paste(removeMatch, collapse = "|"),names(data1))])

###insert NA's into empty cells and get ridd of all empty cases get
df.annotations.genes[df.annotations.genes==""]<-NA
df.annotations.genes <- df.annotations.genes[complete.cases(df.annotations.genes), ]

##get unique  cells 
df.annotations.genes <- distinct(df.annotations.genes, Gene_Names, .keep_all = TRUE)
rownames(df.annotations.genes) <- df.annotations.genes[,1] 
df.annotations.genes.final <- df.annotations.genes[-1]


df.annotations.genes %>% filter(Gene_Names %in% top.cna.genes2) #sanity checks

names(df.annotations.genes.final) ## whats is missing? debug
tools::file_path_sans_ext(tools::file_path_sans_ext(basename(list.cna.files)))
list.cna.files


##comut ready data - for plotting
library(reshape2)
#melt.df.annotations.genes <-  melt(head(df.annotations.genes[-2]),id="Gene_Names",id.name="Gene",variable.name = "sample",value.name ="corrected_copy_number") ##test with header
melt.df.annotations.genes <-  melt(df.annotations.genes[-2],id="Gene_Names",id.name="Gene",variable.name = "sample",value.name ="corrected_copy_number")







top.cna.genes2 <- c("ZNF703","MYC","CCNE1","MDM2","RAD21","ERBB2","NDRG1","BRCA1","BRCA2","TP53","MDM1","PIK3CA","CDH1","GATA3","TERT")
sub.melt.df.annotations.cna.genes = melt.df.annotations.genes[which(melt.df.annotations.genes$Gene_Names %in% top.cna.genes2),]


###sumamry stats on 
sub.melt.df.annotations.cna.genes.summary <- sub.melt.df.annotations.cna.genes %>% group_by(Gene_Names) %>% dplyr::summarise(count=n(),minCN=min(corrected_copy_number),
                                                                                medianCN=median(corrected_copy_number),maxCN=max(corrected_copy_number),
                                                                                CNgain=sum(corrected_copy_number>2 & corrected_copy_number<4),
                                                                                CNgain_Percent=((sum(corrected_copy_number>2 & corrected_copy_number<4))/n())*100,
                                                                                CNamp_Percent=((sum(corrected_copy_number==4))/n())*100,
                                                                                CNhamp_Percent=((sum(corrected_copy_number>4))/n())*100,
                                                                                CNhamp_numb=(sum(corrected_copy_number>4)),
                                                                                CNneutr_Percent=((sum(corrected_copy_number==2))/n())*100,
                                                                                CNdel_Perc=((sum(corrected_copy_number<2))/n())*100,
                                                                                CN_anyCN_above1_Percent=((sum(corrected_copy_number>1))/n())*100,
                                                                                CN_atLeast_gain=(sum(corrected_copy_number>2)),
                                                                                CN_atLeast_gain_Percent=((sum(corrected_copy_number>2))/n())*100)

###plot percentage of individuals with any cn higher tan neutral in the cohort
ggplot(sub.melt.df.annotations.cna.genes.summary, aes(x=as.factor(Gene_Names),y=CN_anyCN_above1_Percent)) + geom_col(color="blue")





#### create CNA event dictionary to define coppy numbers
cna_dictionary <- function(x, na.rm = FALSE){
  ifelse(x==2,"Neutral",
         ifelse(x == 1, "Homodeletion",
                ifelse(x == 4, "Amp",
                       ifelse(x > 4, "High amplification",
                              ifelse(x == 3, "Gain", cnaEvent)))))
}
cna_dictionary(20) 

melt.df.annotations.genes.cna.diction <- melt.df.annotations.genes %>% filter (Gene_Names %in% top.cna.genes2)### subset seppecific genes if you want
melt.df.annotations.genes.cna.diction$cnaEvent <- NA
test.cna.event <- melt.df.annotations.genes.cna.diction %>% mutate(cnaEvent = ifelse(corrected_copy_number==2,"Neutral",ifelse(corrected_copy_number == 1, "Homodeletion",ifelse(corrected_copy_number == 4, "Amplification",ifelse(corrected_copy_number > 4, "High amplification",
                                                                                                                                                                                                                                    ifelse(corrected_copy_number == 3, "Gain", cnaEvent))))))
### errm when i add one more doesn't work anymore
test.cna.event %>% filter(cnaEvent =="High amplification")
test.cna.event.cmut <- test.cna.event[,c(2,1,4)]

write.table(test.cna.event.cmut, file="coMutPlot_ready_data_allSamples.tsv",quote=FALSE, sep='\t', col.names = NA)


####anonymous, 
melt.df.annotations.genes.cna.diction$sample <- gsub("_1000kb","",melt.df.annotations.genes.cna.diction$sample)
melt.df.annotations.genes.cna.diction$sample <- mapvalues(melt.df.annotations.genes.cna.diction$sample, from=c("GBR05352", "GBR06151", "GBR06256", "GBR06127", "GBR00332", "GBR05161", "GBR06243", "GBR05809", 
                                             "GBR05104", "GBR03120", "GBR05085", "GBR05016", "GBR03686", "GBR05110", "GBR05132"),
                         to=c("Participant_01", "Participant_02", "Participant_03", "Participant_04", "Participant_05", "Participant_06", "Participant_07", 
                              "Participant_08", "Participant_09", "Participant_10", "Participant_11", "Participant_12", "Participant_13", "Participant_14", "Participant_15"))

##if you already have ensemble genes
# ##load ensemble genes (gene stable id) from biomart
# ensemblGenes.biomart.hg19 <- read.delim(file = "/Users/samuelahuno/Polak_lab10082019/Ghana_cfDNA/data/mart_export.txt", sep="\t",strip.white = TRUE)
# ##make GRanges object from `ensemblGenes.biomart.hg19` data frame
# gr.ensemblGenes.biomart.hg19 <- makeGRangesFromDataFrame(ensemblGenes.biomart.hg19, start.field="Gene.start..bp.", seqnames.field="Chromosome.scaffold.name",end.field="Gene.end..bp.",keep.extra.columns=TRUE, ignore.strand = TRUE)
# #remove patches
# gr.ensemblGenes.biomart.hg19 <- keepSeqlevels(gr.ensemblGenes.biomart.hg19, value =c(1:22,"X","Y"), pruning.mode = "coarse") ##keep standard chromosomes and sex (X,Y) chrom only
# ##sort seqNmaes
# gr.ensemblGenes.biomart.hg19 <- sortSeqlevels(gr.ensemblGenes.biomart.hg19, X.is.sexchrom=TRUE)  ##reorder chrom to start from chr1 to chrY
# seqlevels(gr.ensemblGenes.biomart.hg19) <- sort(seqlevels(gr.ensemblGenes.biomart.hg19))
# gr.ensemblGenes.biomart.hg19 <- sort(gr.ensemblGenes.biomart.hg19)
# ###esemble genes 
# gr.ensemblGenes.biomart.hg19
# length(gr.ensemblGenes.biomart.hg19) ##[1] 57812



#################################find duplcated cells, only use unque genes ##############################
##find duplicated genes 
#df.annotations.genes[duplicated(df.annotations.genes$Gene_Names),]
##### find number of duplicates
#dim(df.annotations.genes[duplicated(df.annotations.genes$Gene_Names),])[1]
###find all duplicated cells
#df.annotations.genes.duplicated.cells <- df.annotations.genes[duplicated(df.annotations.genes$Gene_Names) | duplicated(df.annotations.genes$Gene_Names, fromLast=TRUE),]

######## try function to get idea of hw values are thrown around################
# addfun <- function(num1,num2){
#   w=num1+num2
#   z=num1*num2
#   df.add <- data.frame(w,z)
#   return(df.add)
# }
# ffff <- addfun(2,5)
# ffff[,1]
########################################## end of function #####################