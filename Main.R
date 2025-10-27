install.packages("BiocManager")
BiocManager::install(version = "3.21")
BiocManager::install("ChAMP")
BiocManager::install("ChAMPdata")
install.packages("/Volumes/Elements/pd-methyl/R-packages/kpmt_0.1.0.tar.gz", repos = NULL, type="source")

library(ChAMP)
library(dplyr)

ppmi_meta <- read.csv("/Volumes/Elements/pd-methyl/ppmi/metaDataIR3.csv", stringsAsFactors = F)

ppmi_meth_120_txt <- read.delim("/Volumes/Elements/pd-methyl/ppmi/Project120_IDATS_n524final_toLONI_030718/PPMI_Meth_n524_for_LONI030718.txt", header=T, sep="\t")

ppmi_meta_baseline <- ppmi_meta[ppmi_meta$CLINICAL_EVENT == "BL",]

clinical_meta <- ppmi_meta_baseline %>% select(PATNO, DIAGNOSIS, GENDER)

sample_sheet <- merge(clinical_meta, ppmi_meth_120_txt, by = "PATNO")

sample_sheet$Basename <- paste0(sample_sheet$Sentrix.ID, "_", sample_sheet$Sentrix.Position)

sample_sheet$Sample_Group <- sample_sheet$DIAGNOSIS

# sample_sheet <- sample_sheet %>% select(Sample_Name = PATNO, Sample_Group, Basename, Gender = GENDER)

sample_sheet <- sample_sheet %>% rename(Sentrix_Position = Sentrix.Position, Sentrix_ID = Sentrix.ID)

write.csv(sample_sheet, "/Volumes/Elements/pd-methyl/ppmi/Project120_IDATS_n524final_toLONI_030718/sample_sheet.csv", row.names = F)

methyl_load <- champ.load(directory = "/Volumes/Elements/pd-methyl/ppmi/Project120_IDATS_n524final_toLONI_030718", arraytype="EPIC")
 

