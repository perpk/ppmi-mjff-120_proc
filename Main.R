BiocManager::install("minfi")
BiocManager::install("limma")
BiocManager::install("DMRcate")
BiocManager::install("missMethyl")

install.packages("openxlsx", dependencies = TRUE)

library(tidyverse)
library(limma)
library(GEOquery)
library(illuminaio)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(DMRcate)
library(missMethyl)
library(openxlsx)

basedir <- "./GSE111629_RAW"
targets <- read.metharray.sheet(basedir, pattern = "GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz")
View(targets)

xlsx <- read.xlsx('./GSE111629_RAW/GPL13534_450K_Manifest_header_Descriptions.xlsx')

idat_files <- list.files("./GSE111629_RAW/", pattern = "GSM.*.idat.gz$")
gsm_ids <- unique(gsub("_(Red|Grn).idat.gz", "", idat_files))
gsm_ids <- unique(gsub("_.*", "", gsm_ids))
print(paste("Found", length(gsm_ids), "unique GSM IDs:", paste(gsm_ids, collapse = ", ")))

gse <- getGEO('GSE111629', GSEMatrix = TRUE, getGPL = FALSE)
pdata <- pData(gse[[1]])

basedir <- './GSE111629_RAW/'
all_idat_files <- list.files(basedir, pattern = "\\.idat\\.gz$", full.names = FALSE)
basenames <- unique(gsub("_(Grn|Red).idat.gz", "", all_idat_files))

targets <- data.frame(
  Sample_Name = gsub("_.*", "", basenames),
  Basename = file.path(basedir, basenames),
  stringsAsFactors = FALSE
)

targets <- merge(targets, pdata, by.x = "Sample_Name", by.y = "row.names", all.x = TRUE)
View(targets)
colnames(targets)
targets$Sample_Group <- factor(targets$`disease state:ch1`)
targets <- targets %>% mutate(Sample_Group = case_when(Sample_Group == "Parkinson's disease (PD)" ~ "PD", Sample_Group == "PD-free control" ~ "Control"))
targets$Sample_Group <- as.factor(targets$Sample_Group)

files_exist <- file.exists(paste0(targets$Basename, "_Grn.idat.gz")) & 
  file.exists(paste0(targets$Basename, "_Red.idat.gz"))
print(paste("All files found:", all(files_exist)))
unique(targets$`disease state:ch1`)

rg_set <- read.metharray.exp(targets = targets, extended = TRUE)
saveRDS(rg_set, "rg_set.rds")

qc <- getQC(preprocessRaw(rg_set))

plotQC(qc)

densityPlot(rg_set, sampGroups = targets$Sample_Group, main = "Beta Value Distribution")

rm(list = setdiff(ls(), "rg_set"))
gc(full = TRUE)

methyl_set <- preprocessNoob(rg_set)
saveRDS(methyl_set, "methyl_set.rds")
rm(rg_set)
gc(full = TRUE)
# Option 2: Functional normalization (good for large studies)
# methyl_set <- preprocessFunnorm(rg_set)

detP <- detectionP(readRDS("rg_set.rds"))
keep <- rowSums(detP < 0.01) == ncol(methyl_set)
methyl_set <- methyl_set[keep, ]
print(paste("Probes after detection p-value filter:", nrow(methyl_set)))
saveRDS(methyl_set, "methyl_set_filtered.rds")
rm(detP, keep)
gc(full = TRUE)

# 2. Remove cross-reactive probes
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
cross_reactive <- ann_450k$Name[ann_450k$Probe_Type %in% c("I") & ann_450k$MASK_general == TRUE]
methyl_set <- readRDS("methyl_set_filtered.rds")
methyl_set <- methyl_set[!rownames(methyl_set) %in% cross_reactive, ]
print(paste("Probes after removing cross-reactive probes:", nrow(methyl_set)))
rm(methyl_set)
gc(full = TRUE)

# 3. Remove probes on sex chromosomes (optional)
methyl_set <- readRDS("methyl_set_filtered.rds")
ann_450k_sub <- ann_450k[rownames(methyl_set), ]
keep_auto <- !(ann_450k_sub$chr %in% c("chrX", "chrY"))
methyl_set <- methyl_set[keep_auto, ]
print(paste("Probes after removing sex chromosomes:", nrow(methyl_set)))
saveRDS(methyl_set, "methyl_set_filtered_chrom.rds")
rm(methyl_set)
gc(full = TRUE)

methyl_set <- readRDS("methyl_set_filtered_chrom.rds")
beta_matrix <- getBeta(methyl_set)
m_values <- getM(methyl_set)

saveRDS(beta_matrix, "beta_matrix.rds")
saveRDS(m_values, "m_values.rds")
rm(methyl_set)
gc(full = TRUE)

# Ensure your group variable is a factor
targets$Sample_Group <- factor(targets$`disease state:ch1`)
targets <- targets %>% mutate(Sample_Group = case_when(Sample_Group == "Parkinson's disease (PD)" ~ "PD", Sample_Group == "PD-free control" ~ "Control"))
targets$Sample_Group <- as.factor(targets$Sample_Group)

View(targets)

# Create design matrix
design <- model.matrix(~0 + targets$Sample_Group)
colnames(design) <- levels(targets$Sample_Group)
# levels(targets$Sample_Group)
# targets$Sample_Group

# View the design matrix to ensure it's correct
design

# Fit linear model
fit <- lmFit(m_values, design)

# Define your contrast (ADJUST THIS BASED ON YOUR GROUP NAMES)
# Example: If your groups are "Control" and "Alzheimer"
cont.matrix <- makeContrasts(Parkinsons_vs_Control = PD - Control, levels = design)

# Apply contrasts
fit2 <- contrasts.fit(fit, cont.matrix)

# Apply empirical Bayes moderation
fit2 <- eBayes(fit2)

# Extract results
results <- topTable(fit2, number = Inf, coef = "Parkinsons_vs_Control")
head(results)

# Check how many probes are significant at FDR < 0.05
sum(results$adj.P.Val < 0.05)
dim(results)

# Add annotation to your results
results$ProbeID <- rownames(results)
ann_sub <- ann_450k[rownames(results), c("chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_Island")]

annotated_results <- merge(results, ann_sub, by.x = "ProbeID", by.y = "row.names")

# Sort by significance
annotated_results <- annotated_results[order(annotated_results$adj.P.Val), ]

# View the top results
head(annotated_results, 20)

# Save the results
write.csv(annotated_results, "differential_methylation_results.csv", row.names = FALSE)
View(annotated_results)

# Volcano plot
library(ggplot2)
volcano_plot <- ggplot(annotated_results, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.6, aes(color = adj.P.Val < 0.05)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  ggtitle("Volcano Plot of Differential Methylation") +
  theme_minimal()
print(volcano_plot)

# Heatmap of top 50 significant probes
library(pheatmap)
top_probes <- annotated_results$ProbeID[1:10]
top_beta <- beta_matrix[rownames(beta_matrix) %in% top_probes, ]

annotation_col <- data.frame(Group = targets$Sample_Group)
rownames(annotation_col) <- colnames(top_beta)

pheatmap(top_beta,
         annotation_col = annotation_col,
         show_rownames = T,
         show_colnames = F,
         main = "Top 10 Differentially Methylated Probes",
         scale = "row")
