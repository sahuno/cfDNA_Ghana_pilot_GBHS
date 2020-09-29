##merge data set cna data and patients informations
library(reshape2)

df.top.cna.ghana <- read.delim(file="/Users/samuelahuno/Polak_lab10082019/Ghana_cfDNA/data/df.top.cna.genes.30xGhana.4pandas.tsv",header=TRUE, sep = "\t", strip.white = TRUE)
pid.infor <- read.delim(file="/Users/samuelahuno/Polak_lab10082019/Ghana_cfDNA/data/final_data/data_Ghana_cfDNA_study_anonymous_07182020.tsv",header=TRUE, sep = "\t", strip.white = TRUE)
pid.infor$PID_anonymous <- gsub("P","Participant_",pid.infor$PID_anonymous)

##reshape for merging
dcast.df.top.cna.ghana <- reshape2::dcast(df.top.cna.ghana, sample ~ category, value.var = "value")





##merge
merged.cna.pid.infor <- cbind(pid.infor,dcast.df.top.cna.ghana)

write.table(merged.cna.pid.infor, file="/Users/samuelahuno/Polak_lab10082019/Ghana_cfDNA/data/final_data/merged.cna.and.pid.infor.ghana.cfDNA.tsv",quote=FALSE, sep='\t', col.names = NA)



melt.df.annotations.genes.cna.inter <- reshape2::dcast(melt.df.annotations.genes.cna.diction, sample ~ Gene_Names, value.var = "corrected_copy_number")
merged.cna.pid.infor.integerCNA <- cbind(pid.infor,dcast.df.top.cna.ghana,melt.df.annotations.genes.cna.inter)
write.table(merged.cna.pid.infor.integerCNA, file="/Users/samuelahuno/Polak_lab10082019/Ghana_cfDNA/data/final_data/merged.Integer_copyNumber.and.pid.infor.ghana.cfDNA.tsv",quote=FALSE, sep='\t', col.names = NA)
