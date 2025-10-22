# Set working directory
setwd("/Users/jasminealqassar/Library/CloudStorage/GoogleDrive-j.alqassar@gwmail.gwu.edu/.shortcut-targets-by-id/1Hi7WIp_ha7vnyQUclvt-AJRl6AF8y1hj/Martin Lab/@LabData/2025_Skipper_WntA/SSS_Differential_Expression/R_working_dir")

# Install necessary packages
#install.packages("tidyverse"); install.packages("stringr"); install.packages("dplyr"); install.packages("DESeq2")

# Load necessary packages
library(tidyverse); library(stringr); library(dplyr); library(DESeq2); library(ggplot2); library(readxl)

#########################
#36h pupal wings FW vs HW


# read the two counts files in 

seqdata_se <- read_tsv("skipper_rnaseq_se.featurecounts.txt", comment="#")
seqdata_pe <- read_tsv("skipper_rnaseq_pe.featurecounts.txt", comment="#")

# remove columns other than Geneid and sample counts
seqdata_se = seqdata_se %>%
  select(-Chr, -Start, -End, -Strand, -Length)

seqdata_pe = seqdata_pe %>%
  select(-Chr, -Start, -End, -Strand, -Length)

seqdata <- seqdata_se %>%
  left_join(seqdata_pe, by = c("Geneid" = "Geneid")) 

# Now rename the "Geneid" column to "Symbol" because that is the appropriate name. NOTE this is the original NCBI symbol so we need to replace it with Geneid before we go further so everything has a unique accession number 
seqdata <- seqdata %>%
  rename("Geneid" = "old_symbol")

#read in the original NCBI gene annotation table 
Symbol_to_Geneid <- read_excel("Ecla_gene_annotation_table_orig.xlsx")

Symbol_to_Geneid <- Symbol_to_Geneid %>%
  select(Symbol, Geneid)

seqdata <- seqdata %>%
  left_join(Symbol_to_Geneid, by = c("old_symbol" = "Symbol"))  %>%
  drop_na() %>%
  relocate("Geneid", .after = "old_symbol") %>%
  select(-"old_symbol")

#rename by sample 
seqdata <- seqdata %>%
  rename_with(~ ifelse(
    . %in% c("Geneid"),
    .,
    str_extract(., "[^/]+(?=_pass2_mappedAligned\\.sortedByCoord\\.out\\.bam)")
  ))

#transform raw data into a matrix of counts
countdata <- seqdata  %>%
  group_by(Geneid) %>%
  column_to_rownames("Geneid") %>%
  as.matrix()

#Remove the genes that are not/lowly expressed:

keep <- rowSums(countdata) > 0 
head(keep)
table(keep)

countdata <-countdata[keep, ]
dim(countdata)
head(countdata)

# QC of counts 
summary(countdata)
boxplot(countdata, las=2)

# Define %out% function

`%out%` <- function(x, y) !(x %in% y)

### before doing analysis by wing type we need to collapse the different compartments of the wing for each individual 

sampleinfo <- read_tsv("sample_info.txt")
sampleinfo

sampleinfo <- sampleinfo %>%
  mutate(collapse_group = paste(individual, wing_type, sep = "_"))

countdata <- as.data.frame(countdata)
countdata$gene_id <- rownames(countdata)
countdata <- as_tibble(countdata)

group_map <- sampleinfo$collapse_group
names(group_map) <- sampleinfo$sample

# now I want to do FW vs HW analysis so I will drop the head sample and collapse compartments 
collapsed_counts <- countdata %>%
  select(-"5d_female_ind1_head") %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "count") %>%
  mutate(group = group_map[sample]) %>%
  group_by(gene_id, group) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = count)

collapsed_counts <- collapsed_counts  %>%
  group_by(gene_id) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()


sampleinfo_hw_fw <- read_tsv("fw_hw_samples.txt")

collapsed_counts_36h_only <- as.data.frame(collapsed_counts) %>%
  rownames_to_column("Geneid") 

collapsed_counts_36h_only <- collapsed_counts_36h_only %>%
  select(Geneid, starts_with("36"))

sampleinfo_hw_fw_36h <- sampleinfo_hw_fw %>%
  filter(grepl("^36", .[[1]]))

collapsed_counts_36h_only <- collapsed_counts_36h_only  %>%
  group_by(Geneid) %>%
  column_to_rownames("Geneid") %>%
  as.matrix()

# order the sample info and count data the same 
sampleinfo_hw_fw_36h <- sampleinfo_hw_fw_36h[match(colnames(collapsed_counts_36h_only),
                                                   sampleinfo_hw_fw_36h[[1]]), ]
# check the sample info and count data are in the same order
all(colnames(collapsed_counts_36h_only) == sampleinfo_hw_fw_36h[[1]])


dds <- DESeqDataSetFromMatrix(
  countData = collapsed_counts_36h_only,
  colData = sampleinfo_hw_fw_36h,
  design = ~ wing_type
)

ddsObj.raw <- DESeq(dds)
ddsObj <- estimateSizeFactors(ddsObj.raw)
colData(ddsObj.raw)
colData(ddsObj)
ddsObj <- DESeq(ddsObj.raw)
res <- results(ddsObj, alpha=0.05) #adding p-value cutoff
res

resultsNames(ddsObj) 



#now the top genes diff exp between FW and HW
wing_type_HW_vs_FW_36h <- results(ddsObj, alpha = 0.05, contrast=c("wing_type","FW","HW"))

sum(wing_type_HW_vs_FW_36h$padj < 0.05, na.rm = TRUE)
allGenesHW_vs_FW_36h <- as.data.frame(wing_type_HW_vs_FW_36h) %>%
  rownames_to_column("Geneid") %>% 
  arrange(padj) 

manual_annotations <- read_excel("working_Ecla_annotation_table_with_flybase_names.xlsx")


allGenesHW_vs_FW_36h <- allGenesHW_vs_FW_36h %>%
  left_join(manual_annotations, by = c("Geneid" = "Geneid"))

# Now make a vol plot for FW vs HW

allGenesHW_vs_FW_36h <- allGenesHW_vs_FW_36h %>%
  mutate(direction = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "up",
    padj < 0.05 & log2FoldChange < 0 ~ "down",
    TRUE ~ "ns"
  ))

threshold = 5
vol_plot <- allGenesHW_vs_FW_36h %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = direction)) +  
  geom_point() +  
  scale_color_manual(values = c("down" = "#26b3ff", "ns" = "grey", "up" = "#bb0c00"),  
                     labels = c("down" = "HW", "ns" = "Not significant", "up" = "FW")) +
  theme(legend.title=element_blank())+
  ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) +
  #geom_text(aes(label =Geneid), color = "black", size = 3)
  
  geom_text(aes(label = ifelse(-log10(padj) > threshold, Symbol, "")),
            color = "black", size = 3, nudge_y = 1)  
vol_plot  


library(ggrepel)

vol_plot <- allGenesHW_vs_FW_36h %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = direction)) +  
  scale_y_continuous( limits=c(-2, 70), expand=c(0,0)) +
  scale_x_continuous( limits=c(-15, 10), expand=c(0,0)) +
  geom_point(data = filter(allGenesHW_vs_FW_36h, direction == "ns"),
             aes(x = log2FoldChange, y = -log10(padj), color = direction),
             show.legend = TRUE) +
  geom_point(data = filter(allGenesHW_vs_FW_36h, direction != "ns"),
             aes(x = log2FoldChange, y = -log10(padj), color = direction)) +
  
  scale_color_manual(
    values = c("down" = "#26b3ff", "up" = "#bb0c00", "ns" = "grey"),
    labels = c("down" = "HW", "up" = "FW", "ns" = "Non-Significant")
  ) +
  theme_classic() +
  theme(legend.title=element_blank())+
  ggtitle("FW vs HW in 36hr pupae") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_label_repel(
    data = filter(allGenesHW_vs_FW_36h,  -log10(padj) > 20 | padj < 0.05 & log2FoldChange > 2),
    aes(label = gene_label), # here edit where you are pulling the label from 
    color = "black",
    size = 3,
    nudge_y = 1,                  # Nudges label vertically
    max.overlaps = 8, #can increase to Inf
    segment.color = "black",
    segment.size = 0.3,
    fill = alpha("white", 0),
    label.r = unit(0, "lines"),     # No rounded corners
    label.size = 0,                 # No label border
    point.padding = 0,          # Ensures label does not overlap point
    box.padding = 0             # Controls spacing around labels
    
  ) 

vol_plot

#rotated volcano plot 

vol_plot <- allGenesHW_vs_FW_36h %>%
  ggplot(aes(x = -log10(padj), y = log2FoldChange, color = direction)) +  
  #scale_y_continuous( limits=c(-2, 70), expand=c(0,0)) +
  scale_y_continuous( limits=c(-5, 15), expand=c(0,0)) +
  geom_point(data = filter(allGenesHW_vs_FW_36h, direction == "ns"),
             aes(x = -log10(padj), y = log2FoldChange, color = direction),
             show.legend = TRUE) +
  geom_point(data = filter(allGenesHW_vs_FW_36h, direction != "ns"),
             aes(x = -log10(padj), y = log2FoldChange, color = direction)) +
  
  scale_color_manual(
    values = c("down" = "#26b3ff", "up" = "#bb0c00", "ns" = "grey"),
    labels = c("down" = "HW", "up" = "FW", "ns" = "Non-Significant")
  ) +
  theme_classic() +
  theme(legend.title=element_blank())+
  ggtitle("FW vs HW in 36hr pupae") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_label_repel(
    data = filter(allGenesHW_vs_FW_36h,  -log10(padj) > 20 | padj < 0.05 & log2FoldChange < -2 | padj < 0.05 & log2FoldChange > 10),
    aes(label = gene_label), # here edit where you are pulling the label from 
    color = "black",
    size = 3,
    nudge_y = 1,                  # Nudges label vertically
    max.overlaps = 8, #can increase to Inf
    segment.color = "black",
    segment.size = 0.3,
    fill = alpha("white", 0),
    label.r = unit(0, "lines"),     # No rounded corners
    label.size = 0,                 # No label border
    point.padding = 0,          # Ensures label does not overlap point
    box.padding = 0             # Controls spacing around labels
    
  ) 

vol_plot

# re-arrange stuff before writing final table 

allGenesHW_vs_FW_36h <- allGenesHW_vs_FW_36h %>%
  relocate("Symbol", .after = "Geneid") %>%
  relocate("Name", .after = "Symbol") %>%
  relocate("FB_gene_symbol", .after = "Name") %>%
  relocate("gene_fullname", .after = "FB_gene_symbol")


write.table(allGenesHW_vs_FW_36h, file="DESeq2_Results_36hr_HW_vs_FW.tsv", quote=F, sep="\t", row.names=FALSE, na="")

#get norm counts from deseq
counts_ddsObj_FW_HW <- counts(ddsObj, normalized=TRUE) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Geneid")

manual_annotations <- read_excel("working_Ecla_annotation_table_with_flybase_names.xlsx")
counts_ddsObj_FW_HW  <- counts_ddsObj_FW_HW  %>%
  left_join(manual_annotations, by = c("Geneid" = "Geneid")) %>%
  relocate("Symbol", .after = "Geneid") %>%
  relocate("Name", .after = "Symbol") %>%
  relocate("FB_gene_symbol", .after = "Name") %>%
  relocate("gene_fullname", .after = "FB_gene_symbol")

write.table(counts_ddsObj_FW_HW , file="DESeq2_normcounts_36hr_FW_vs_HW.tsv", quote=F, sep="\t", row.names=FALSE, na="")

# MA plot 
label_data <- allGenesHW_vs_FW_36h %>%
  filter(padj < 0.05 & !is.na(gene_label))


MA_plot <- allGenesHW_vs_FW_36h %>%
  ggplot(aes(x = log2(baseMean), y = log2FoldChange, color = direction)) +  
  scale_y_continuous( limits=c(-5, 15), expand=c(0,0)) +
  scale_x_continuous( limits=c(-4, 25), expand=c(0,0)) +
  geom_point(data = filter(allGenesHW_vs_FW_36h, direction == "ns"),
             aes(x = log2(baseMean), y = log2FoldChange, color = direction),
             show.legend = TRUE) +
  geom_point(data = filter(allGenesHW_vs_FW_36h, direction != "ns"),
             aes(x = log2(baseMean), y = log2FoldChange, color = direction)) +
  
  scale_color_manual(
    values = c("down" = "#26b3ff",  "up" = "#bb0c00", "ns" = "grey"),
    labels = c("down" = "HW", "up" = "FW", "ns" = "Non-Significant")
  ) +
  theme_classic() +
  theme(legend.title=element_blank())+
  ggtitle("MA Plot FW vs HW in 36hr pupae") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_label_repel(
    data = filter(allGenesHW_vs_FW_36h, padj < 0.05 & log2(baseMean) > 15 | padj < 0.05 & log2FoldChange < -2 | padj < 0.05 & log2FoldChange > 10),
    aes(label = gene_label),
    color = "black",
    size = 3,
    nudge_y = 1,                  # Nudges label vertically
    max.overlaps = 8,
    segment.color = "black",
    segment.size = 0.3,
    fill = alpha("white", 0),
    label.r = unit(0, "lines"),     # No rounded corners
    label.size = 0,                 # No label border
    point.padding = 0,          # Ensures label does not overlap point
    box.padding = 0             # Controls spacing around labels
    
  ) 
MA_plot




###### 36hr compartment analysis 

# read the two counts files in 

seqdata_se <- read_tsv("skipper_rnaseq_se.featurecounts.txt", comment="#")
seqdata_pe <- read_tsv("skipper_rnaseq_pe.featurecounts.txt", comment="#")

# remove columns other than Geneid and sample counts
seqdata_se = seqdata_se %>%
  select(-Chr, -Start, -End, -Strand, -Length)

seqdata_pe = seqdata_pe %>%
  select(-Chr, -Start, -End, -Strand, -Length)

seqdata <- seqdata_se %>%
  left_join(seqdata_pe, by = c("Geneid" = "Geneid")) 

# Now rename the "Geneid" column to "Symbol" because that is the appropriate name. NOTE this is the original NCBI symbol so we need to replace it with Geneid before we go further so everything has a unique accession number 
seqdata <- seqdata %>%
  rename("Geneid" = "old_symbol")

#read in the original NCBI gene annotation table 
Symbol_to_Geneid <- read_excel("Ecla_gene_annotation_table_orig.xlsx")

Symbol_to_Geneid <- Symbol_to_Geneid %>%
  select(Symbol, Geneid)

seqdata <- seqdata %>%
  left_join(Symbol_to_Geneid, by = c("old_symbol" = "Symbol"))  %>%
  drop_na() %>%
  relocate("Geneid", .after = "old_symbol") %>%
  select(-"old_symbol")

#rename by sample 
seqdata <- seqdata %>%
  rename_with(~ ifelse(
    . %in% c("Geneid"),
    .,
    str_extract(., "[^/]+(?=_pass2_mappedAligned\\.sortedByCoord\\.out\\.bam)")
  ))

#transform raw data into a matrix of counts
countdata <- seqdata  %>%
  group_by(Geneid) %>%
  column_to_rownames("Geneid") %>%
  as.matrix()

#Remove the genes that are not/lowly expressed:

keep <- rowSums(countdata) > 0 
head(keep)
table(keep)

countdata <-countdata[keep, ]
dim(countdata)
head(countdata)

# QC of counts 
summary(countdata)
boxplot(countdata, las=2)

# Define %out% function

`%out%` <- function(x, y) !(x %in% y)

countdata <- as.data.frame(countdata)
countdata$gene_id <- rownames(countdata)
countdata <- as_tibble(countdata)

sampleinfo <- read_tsv("sample_info.txt")

countdata_silver <- countdata %>%
  rename("gene_id" = "Geneid") 


countdata_silver_36h_only <- countdata_silver  %>%
  select(Geneid, starts_with("36")) 

sampleinfo_silver_36h <- sampleinfo %>%
  filter(grepl("^36", .[[1]])) 

countdata_silver_36h_only <- countdata_silver_36h_only  %>%
  group_by(Geneid) %>%
  column_to_rownames("Geneid") %>%
  as.matrix()

# order the sample info and count data the same 
sampleinfo_silver_36h <- sampleinfo_silver_36h[match(colnames(countdata_silver_36h_only),
                                                     sampleinfo_silver_36h[[1]]), ]
# check the sample info and count data are in the same order
all(colnames(countdata_silver_36h_only) == sampleinfo_silver_36h[[1]])


dds <- DESeqDataSetFromMatrix(
  countData = countdata_silver_36h_only,
  colData = sampleinfo_silver_36h,
  design = ~ compartment
)

ddsObj.raw <- DESeq(dds)
ddsObj <- estimateSizeFactors(ddsObj.raw)
colData(ddsObj.raw)
colData(ddsObj)
ddsObj <- DESeq(ddsObj.raw)
res <- results(ddsObj, alpha=0.05) #adding p-value cutoff
res

resultsNames(ddsObj) 


#relevel for HWM to be the control 
ddsObj$compartment <- relevel(ddsObj$compartment, ref = "HWM")

ddsObj <- DESeq(ddsObj)
resultsNames(ddsObj) 

plotCounts(ddsObj, "LOC140755181", "compartment", normalized=TRUE)
plotCounts(ddsObj, "LOC140744733", "compartment", normalized=TRUE)

#get norm counts from deseq
counts_ddsObj_HWM_control <- counts(ddsObj, normalized=TRUE) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Geneid")

manual_annotations <- read_excel("working_Ecla_annotation_table_with_flybase_names.xlsx")
counts_ddsObj_HWM_control <- counts_ddsObj_HWM_control %>%
  left_join(manual_annotations, by = c("Geneid" = "Geneid")) %>%
  relocate("Symbol", .after = "Geneid") %>%
  relocate("Name", .after = "Symbol") %>%
  relocate("FB_gene_symbol", .after = "Name") %>%
  relocate("gene_fullname", .after = "FB_gene_symbol")

write.table(counts_ddsObj_HWM_control, file="DESeq2_normcounts_36hr_HWM_vs_compartments.tsv", quote=F, sep="\t", row.names=FALSE, na="")


# Load results
res_hwd_hwm <- results(ddsObj, alpha = 0.05, contrast = c("compartment", "HWD", "HWM"))
res_hwp_hwm <- results(ddsObj, alpha = 0.05, contrast = c("compartment", "HWP", "HWM"))

manual_annotations <- read_excel("working_Ecla_annotation_table_with_flybase_names.xlsx")

sum(res_hwp_hwm$padj < 0.05, na.rm = TRUE)

# Make result tables for all pairwise comparisons 
# HWP vs HWM 
allGenes_HWP_VS_HWM_36h <- as.data.frame(res_hwp_hwm) %>%
  rownames_to_column("Geneid") %>% 
  arrange(padj) %>%
  left_join(manual_annotations, by = c("Geneid" = "Geneid")) %>%
  mutate(direction = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "up",
    padj < 0.05 & log2FoldChange < 0 ~ "down",
    TRUE ~ "ns"
  )) %>%
  relocate("Symbol", .after = "Geneid") %>%
  relocate("Name", .after = "Symbol") %>%
  relocate("FB_gene_symbol", .after = "Name") %>%
  relocate("gene_fullname", .after = "FB_gene_symbol")

write.table(allGenes_HWP_VS_HWM_36h, file="DESeq2_Results_36hr_HWM_vs_HWP.tsv", quote=F, sep="\t", row.names=FALSE, na="")
#HWD
allGenes_HWD_VS_HWM_36h <- as.data.frame(res_hwd_hwm) %>%
  rownames_to_column("Geneid") %>% 
  arrange(padj) %>%
  left_join(manual_annotations, by = c("Geneid" = "Geneid")) %>%
  mutate(direction = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "up",
    padj < 0.05 & log2FoldChange < 0 ~ "down",
    TRUE ~ "ns"
  )) %>%
  relocate("Symbol", .after = "Geneid") %>%
  relocate("Name", .after = "Symbol") %>%
  relocate("FB_gene_symbol", .after = "Name") %>%
  relocate("gene_fullname", .after = "FB_gene_symbol")

write.table(allGenes_HWD_VS_HWM_36h, file="DESeq2_Results_36hr_HWM_vs_HWD.tsv", quote=F, sep="\t", row.names=FALSE, na="")


#relevel for HWD to be the control to get HWP vs HWD 
ddsObj$compartment <- relevel(ddsObj$compartment, ref = "HWD")

ddsObj <- DESeq(ddsObj)
resultsNames(ddsObj) 


#get norm counts from deseq
counts_ddsObj_HWD_control <- counts(ddsObj, normalized=TRUE) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Geneid")

manual_annotations <- read_excel("working_Ecla_annotation_table_with_flybase_names.xlsx")
counts_ddsObj_HWD_control <- counts_ddsObj_HWD_control %>%
  left_join(manual_annotations, by = c("Geneid" = "Geneid")) %>%
  relocate("Symbol", .after = "Geneid") %>%
  relocate("Name", .after = "Symbol") %>%
  relocate("FB_gene_symbol", .after = "Name") %>%
  relocate("gene_fullname", .after = "FB_gene_symbol")

write.table(counts_ddsObj_HWD_control, file="DESeq2_normcounts_36hr_HWD_vs_compartments.tsv", quote=F, sep="\t", row.names=FALSE, na="")


# Load results
res_hwp_hwd <- results(ddsObj, alpha = 0.05, contrast = c("compartment", "HWP", "HWD"))

manual_annotations <- read_excel("working_Ecla_annotation_table_with_flybase_names.xlsx")

sum(res_hwp_hwd$padj < 0.05, na.rm = TRUE)

# Make result tables for all pairwise comparisons 
# HWP vs HWM 
allGenes_HWP_VS_HWD_36h <- as.data.frame(res_hwp_hwd) %>%
  rownames_to_column("Geneid") %>% 
  arrange(padj) %>%
  left_join(manual_annotations, by = c("Geneid" = "Geneid")) %>%
  mutate(direction = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "up",
    padj < 0.05 & log2FoldChange < 0 ~ "down",
    TRUE ~ "ns"
  )) %>%
  relocate("Symbol", .after = "Geneid") %>%
  relocate("Name", .after = "Symbol") %>%
  relocate("FB_gene_symbol", .after = "Name") %>%
  relocate("gene_fullname", .after = "FB_gene_symbol")

write.table(allGenes_HWP_VS_HWD_36h , file="DESeq2_Results_36hr_HWD_vs_HWP.tsv", quote=F, sep="\t", row.names=FALSE, na="")

# Now do the three FW comparisons by releveling
# first FWD vs FWM and FWP vs FWM
ddsObj$compartment <- relevel(ddsObj$compartment, ref = "FWM")

ddsObj <- DESeq(ddsObj)
resultsNames(ddsObj) 


#get norm counts from deseq
counts_ddsObj_FWM_control <- counts(ddsObj, normalized=TRUE) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Geneid")

manual_annotations <- read_excel("working_Ecla_annotation_table_with_flybase_names.xlsx")
counts_ddsObj_FWM_control <- counts_ddsObj_FWM_control %>%
  left_join(manual_annotations, by = c("Geneid" = "Geneid")) %>%
  relocate("Symbol", .after = "Geneid") %>%
  relocate("Name", .after = "Symbol") %>%
  relocate("FB_gene_symbol", .after = "Name") %>%
  relocate("gene_fullname", .after = "FB_gene_symbol")

write.table(counts_ddsObj_FWM_control, file="DESeq2_normcounts_36hr_FWM_vs_compartments.tsv", quote=F, sep="\t", row.names=FALSE, na="")


# Load results
res_fwp_fwm <- results(ddsObj, alpha = 0.05, contrast = c("compartment", "FWP", "FWM"))
res_fwd_fwm <- results(ddsObj, alpha = 0.05, contrast = c("compartment", "FWD", "FWM"))


manual_annotations <- read_excel("working_Ecla_annotation_table_with_flybase_names.xlsx")

sum(res_fwp_fwm$padj < 0.05, na.rm = TRUE)

# Make result tables for all pairwise comparisons 
# FWP vs FWM 
allGenes_FWP_VS_FWM_36h <- as.data.frame(res_fwp_fwm) %>%
  rownames_to_column("Geneid") %>% 
  arrange(padj) %>%
  left_join(manual_annotations, by = c("Geneid" = "Geneid")) %>%
  mutate(direction = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "up",
    padj < 0.05 & log2FoldChange < 0 ~ "down",
    TRUE ~ "ns"
  )) %>%
  relocate("Symbol", .after = "Geneid") %>%
  relocate("Name", .after = "Symbol") %>%
  relocate("FB_gene_symbol", .after = "Name") %>%
  relocate("gene_fullname", .after = "FB_gene_symbol")

write.table(allGenes_FWP_VS_FWM_36h , file="DESeq2_Results_36hr_FWP_vs_FWM.tsv", quote=F, sep="\t", row.names=FALSE, na="")


# Make result tables for all pairwise comparisons 
# FWD vs FWM 
allGenes_FWD_VS_FWM_36h <- as.data.frame(res_fwd_fwm) %>%
  rownames_to_column("Geneid") %>% 
  arrange(padj) %>%
  left_join(manual_annotations, by = c("Geneid" = "Geneid")) %>%
  mutate(direction = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "up",
    padj < 0.05 & log2FoldChange < 0 ~ "down",
    TRUE ~ "ns"
  )) %>%
  relocate("Symbol", .after = "Geneid") %>%
  relocate("Name", .after = "Symbol") %>%
  relocate("FB_gene_symbol", .after = "Name") %>%
  relocate("gene_fullname", .after = "FB_gene_symbol")

write.table(allGenes_FWD_VS_FWM_36h , file="DESeq2_Results_36hr_FWD_vs_FWM.tsv", quote=F, sep="\t", row.names=FALSE, na="")


# revlevel to get FWP vs FWD
ddsObj$compartment <- relevel(ddsObj$compartment, ref = "FWD")

ddsObj <- DESeq(ddsObj)
resultsNames(ddsObj) 


#get norm counts from deseq
counts_ddsObj_FWD_control <- counts(ddsObj, normalized=TRUE) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Geneid")

manual_annotations <- read_excel("working_Ecla_annotation_table_with_flybase_names.xlsx")
counts_ddsObj_FWD_control <- counts_ddsObj_FWD_control %>%
  left_join(manual_annotations, by = c("Geneid" = "Geneid")) %>%
  relocate("Symbol", .after = "Geneid") %>%
  relocate("Name", .after = "Symbol") %>%
  relocate("FB_gene_symbol", .after = "Name") %>%
  relocate("gene_fullname", .after = "FB_gene_symbol")

write.table(counts_ddsObj_FWD_control, file="DESeq2_normcounts_36hr_FWD_vs_compartments.tsv", quote=F, sep="\t", row.names=FALSE, na="")


# Load results
res_fwp_fwd <- results(ddsObj, alpha = 0.05, contrast = c("compartment", "FWP", "FWD"))


manual_annotations <- read_excel("working_Ecla_annotation_table_with_flybase_names.xlsx")

sum(res_fwp_fwd$padj < 0.05, na.rm = TRUE)

# Make result tables for all pairwise comparisons 
# FWP vs FWM 
allGenes_FWP_VS_FWD_36h <- as.data.frame(res_fwp_fwd) %>%
  rownames_to_column("Geneid") %>% 
  arrange(padj) %>%
  left_join(manual_annotations, by = c("Geneid" = "Geneid")) %>%
  mutate(direction = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "up",
    padj < 0.05 & log2FoldChange < 0 ~ "down",
    TRUE ~ "ns"
  )) %>%
  relocate("Symbol", .after = "Geneid") %>%
  relocate("Name", .after = "Symbol") %>%
  relocate("FB_gene_symbol", .after = "Name") %>%
  relocate("gene_fullname", .after = "FB_gene_symbol")

write.table(allGenes_FWP_VS_FWD_36h , file="DESeq2_Results_36hr_FWP_vs_FWD.tsv", quote=F, sep="\t", row.names=FALSE, na="")



# Now make a heatmap highlighting genes DE in at least one the 3 HW comparisons 

#first make a list of the LOC gene ids of differentially expressed genes 

sigGenes_HWD_vs_HWM <- allGenes_HWD_VS_HWM_36h[!is.na(allGenes_HWD_VS_HWM_36h$padj) & allGenes_HWD_VS_HWM_36h$padj < 0.05, ]
sigGenes_HWP_vs_HWM <- allGenes_HWP_VS_HWM_36h[!is.na(allGenes_HWP_VS_HWM_36h$padj) & allGenes_HWP_VS_HWM_36h$padj < 0.05, ]
sigGenes_HWP_vs_HWD <- allGenes_HWP_VS_HWD_36h[!is.na(allGenes_HWP_VS_HWD_36h$padj) & allGenes_HWP_VS_HWD_36h$padj < 0.05, ]



sigGenes_HW_df <- bind_rows(sigGenes_HWD_vs_HWM, sigGenes_HWP_vs_HWM, sigGenes_HWP_vs_HWD) #keep duplicates for now

# Now do the same for the FW Analyses 

sigGenes_FWD_vs_FWM <- allGenes_FWD_VS_FWM_36h[!is.na(allGenes_FWD_VS_FWM_36h$padj) & allGenes_FWD_VS_FWM_36h$padj < 0.05, ]
sigGenes_FWP_vs_FWM <- allGenes_FWP_VS_FWM_36h[!is.na(allGenes_FWP_VS_FWM_36h$padj) & allGenes_FWP_VS_FWM_36h$padj < 0.05, ]
sigGenes_FWP_vs_FWD <- allGenes_FWP_VS_FWD_36h[!is.na(allGenes_FWP_VS_FWD_36h$padj) & allGenes_FWP_VS_FWD_36h$padj < 0.05, ]

sigGenes_FW_df <- bind_rows(sigGenes_FWD_vs_FWM, sigGenes_FWP_vs_FWM, sigGenes_FWP_vs_FWD)  #keep duplicates for now


# combine to one dataset HW and FW df 
sigGenes_HW_FW_df <- bind_rows(sigGenes_HW_df, sigGenes_FW_df) 

sigGenes_HW_FW_df <- sigGenes_HW_FW_df[!is.na(sigGenes_HW_FW_df$padj) & sigGenes_HW_FW_df$padj < 0.05, ]


sigGenes_HW_FW_df  <- sigGenes_HW_FW_df %>%
  select(Geneid) %>%
  distinct(Geneid, .keep_all = TRUE)


sigGenes_FW_df <- sigGenes_FW_df %>%
  select(Geneid) %>%
distinct(Geneid, .keep_all = TRUE)

sigGenes_HW_df <- sigGenes_HW_df %>%
  select(Geneid) %>%
  distinct(Geneid, .keep_all = TRUE)

# now re-run the DEseq object just to get the counts for the heatmap

sampleinfo_all <- read_tsv("sample_info_new.txt")
sampleinfo_all <- sampleinfo_all %>% 
  filter(!grepl("^5d", .[[1]])) %>% 
  filter(!grepl("^48h_male_ind2", .[[1]])) %>% 
  filter(!grepl("^48h_male_ind3_FWP", .[[1]])) %>% 
  filter(!grepl("^48h_male_ind3_FWM", .[[1]])) %>% 
  filter(!grepl("^48h_male_ind4_FWD", .[[1]])) %>% 
  filter(!grepl("^48h_male_ind5_FWD", .[[1]])) %>% 
  filter(!grepl("^48h_male_ind1", .[[1]]))

countdata_all <- countdata %>%
  rename("gene_id" = "Geneid") %>%
  select(-starts_with("5d")) %>%
  select(-starts_with("48h_male_ind2"))%>%
  select(-starts_with("48h_male_ind3_FWP"))%>%
  select(-starts_with("48h_male_ind3_FWM"))%>%
  select(-starts_with("48h_male_ind4_FWD")) %>%
  select(-starts_with("48h_male_ind5_FWD")) %>%
  select(-starts_with("48h_male_ind1")) 
  
  
  
countdata_all <- countdata_all %>%
  group_by(Geneid) %>%
  column_to_rownames("Geneid") %>%
  as.matrix()

# order the sample info and count data the same 
sampleinfo_all <- sampleinfo_all[match(colnames(countdata_all),
                                             sampleinfo_all[[1]]), ]
# check the sample info and count data are in the same order
all(colnames(countdata_all) == sampleinfo_all[[1]])


dds <- DESeqDataSetFromMatrix(
  countData = countdata_all,
  colData = sampleinfo_all,
  design = ~ compartment
)

ddsObj.raw <- DESeq(dds)
ddsObj <- estimateSizeFactors(ddsObj.raw)
colData(ddsObj.raw)
colData(ddsObj)
ddsObj <- DESeq(ddsObj.raw)
res <- results(ddsObj, alpha=0.05) #adding p-value cutoff
res

resultsNames(ddsObj) 

ddsObj$compartment <- relevel(ddsObj$compartment, ref = "36h_HWM")

ddsObj <- DESeq(ddsObj)
resultsNames(ddsObj) 



#Make a PCA

vsd <- vst(ddsObj, blind=FALSE)

pcaData<- plotPCA(vsd, intgroup=c("dev_stage", "sex", "wing_type"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=dev_stage, shape=sex, alpha =wing_type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_color_manual(values = c("L5" = "#1babe2", "36h" = "#f89c30", "48h" = "#b0207c")) +
  scale_alpha_manual(values = c("HW" = 0.3, "FW" = 1))

#get norm counts from deseq
counts_ddsObj_for_heatmap <- counts(ddsObj, normalized=TRUE) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Geneid")

manual_annotations <- read_excel("working_Ecla_annotation_table_with_flybase_names.xlsx")
counts_ddsObj_for_heatmap <- counts_ddsObj_for_heatmap %>%
  left_join(manual_annotations, by = c("Geneid" = "Geneid")) %>%
  relocate("Symbol", .after = "Geneid") %>%
  relocate("Name", .after = "Symbol") %>%
  relocate("FB_gene_symbol", .after = "Name") %>%
  relocate("gene_fullname", .after = "FB_gene_symbol")

write.table(counts_ddsObj_for_heatmap, file="DESeq2_normcounts_for_heatmap.tsv", quote=F, sep="\t", row.names=FALSE, na="")


# Now make individual count plots for specific genes (for example the gene WntA)
plotCounts(ddsObj, "LOC140755181", normalized=TRUE, intgroup = "compartment") #default DESeq plot

# extract counts for WntA
WntA_normalized_counts <- counts(ddsObj, normalized=TRUE)["LOC140755181", ]
normcounts <- counts(ddsObj, normalized=TRUE)
WntA_counts_plot_data <- data.frame(
  gene_counts = WntA_normalized_counts,
  condition = colData(dds)$compartment )


WntA_counts_plot_data$condition <- factor(WntA_counts_plot_data$condition, levels = c("L5_FW", "L5_HW", "36h_FWP", "36h_FWM", "36h_FWD", "36h_HWP", "36h_HWM", "36h_HWD", "48h_FWP", "48h_FWM", "48h_FWD", "48h_HWP", "48h_HWM", "48h_HWD"))
WntA_counts_plot_data <- WntA_counts_plot_data %>%
  mutate(group = case_when(
    grepl("^36h_FW", condition) ~ "36h_FW",
    grepl("^36h_HW", condition) ~ "36h_HW",
    grepl("^48h_FW", condition) ~ "48h_FW",
    grepl("^48h_HW", condition) ~ "48h_HW",
    grepl("^L5_HW", condition) ~ "L5_HW",
    grepl("^L5_FW", condition) ~ "L5_FW",
  ))

WntA_counts_plot_data <- WntA_counts_plot_data %>%
  mutate(color_col = case_when(
    grepl("^36h_FW", condition) ~ "#4477AA",
    grepl("^36h_HW", condition) ~ "#66CCEE",
    grepl("^48h_FW", condition) ~ "#CCBB44",
    grepl("^48h_HW", condition) ~ "#EE6677",
  ))

ggplot(WntA_counts_plot_data, aes(x = condition, y = gene_counts)) +
  #scale_y_continuous( limits=c(0, 3000), expand=c(0,0)) +
  geom_point() +
  #theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  xlab("Compartment") +
  ylab("Counts") +
  ggtitle("WntA") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_summary(fun =mean,geom="line",lwd=1, aes(group=group, color= color_col) ) +
  scale_color_identity(guide = "none") +   # use the hex codes as-is 
  theme(legend.position = "none")  # remove the legend


######## A for loop to generate many gene count plots at once 
normcounts <- counts(ddsObj, normalized=TRUE)
# Now try to automate this in a for loop to make many plots at once 
gene_list_for_count_plot <- read_tsv("gene_list.txt")

gene_list_for_count_plot <- gene_list_for_count_plot %>%
  select(Geneid, Gene_label)

 for (i in 1:nrow(gene_list_for_count_plot)) {
   gene_id <- gene_list_for_count_plot$Geneid[i] #create list
   gene_name <- gene_list_for_count_plot$Gene_label[i] #create list
   
   gene_counts <- normcounts[gene_id, ] #get the counts for only that gene 
  
    counts_plot_data <- data.frame(
     gene_counts = gene_counts,
     condition   = colData(ddsObj)$compartment
   )
   
    counts_plot_data$condition <- factor(
      counts_plot_data$condition,
      levels = c("L5_FW", "L5_HW", 
                 "36h_FWP", "36h_FWM", "36h_FWD", 
                 "36h_HWP", "36h_HWM", "36h_HWD",
                 "48h_FWP", "48h_FWM", "48h_FWD", 
                 "48h_HWP", "48h_HWM", "48h_HWD"))
    
    counts_plot_data <- counts_plot_data %>%
      mutate(group = case_when(
        grepl("^36h_FW", condition) ~ "36h_FW",
        grepl("^36h_HW", condition) ~ "36h_HW",
        grepl("^48h_FW", condition) ~ "48h_FW",
        grepl("^48h_HW", condition) ~ "48h_HW",
        grepl("^L5_HW", condition)  ~ "L5_HW",
        grepl("^L5_FW", condition)  ~ "L5_FW"
      )) %>%
      mutate(color_col = case_when(
        grepl("^36h_FW", condition) ~ "#4477AA",
        grepl("^36h_HW", condition) ~ "#66CCEE",
        grepl("^48h_FW", condition) ~ "#CCBB44",
        grepl("^48h_HW", condition) ~ "#EE6677"
      ))
    plot <- ggplot(counts_plot_data, aes(x = condition, y = gene_counts)) +
      geom_point() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
      xlab("Compartment") +
      ylab("Counts") +
      ggtitle(gene_name) +
      theme(plot.title = element_text(hjust = 0.5)) +
      stat_summary(fun = mean, geom = "line", lwd = 1,
                   aes(group = group, color = color_col)) +
      scale_color_identity(guide = "none") +
      theme(legend.position = "none") +
      scale_y_continuous(
        breaks = function(x) pretty(x),         # generate standard breaks
        labels = function(b) format(b, trim=TRUE)  # force all breaks (incl. first) to show to add more labels 
      )
    
    # save each counts plot per gene
    ggsave(filename = paste0(gene_name, "_counts_plot.eps"), plot = plot, width = 1492, height = 1292, units="px") #size specified in pixels 
 }

# What I ended up doing for the final figure 

normcounts <- counts(ddsObj, normalized=TRUE)

# Now try to make a large matrix with the plots aligned with 4 plots per row
gene_list_for_count_plot <- read_tsv("new_gene_list.txt")

gene_list_for_count_plot <- gene_list_for_count_plot %>%
  select(Geneid, Gene_label)

counts_plot_data <- data.frame(normcounts, check.names = FALSE)
counts_plot_data <- counts_plot_data  %>%
  rownames_to_column("Geneid")

counts_plot_data <- counts_plot_data %>%
  filter(Geneid %in% gene_list_for_count_plot$Geneid)  


# Make a longer count matrix with a row for the count data for each sample for the gene included in the count plots 
counts_long <- counts_plot_data %>%
  pivot_longer(
    cols = -Geneid,
    names_to = "sample_id",
    values_to = "gene_counts"
  )

#Add the gene labels you want for the plot
counts_long <- counts_long %>%
  left_join(gene_list_for_count_plot, by = "Geneid")

#Get the sample info from DESeq2 object
sample_info <- as.data.frame(colData(ddsObj)) %>%
  rownames_to_column("sample_id")

counts_long <- counts_long %>%
  left_join(sample_info, by = "sample_id")

# Now add groupings that you want and color codes for the regression lines 
counts_long <- counts_long %>%
  mutate(condition = factor(compartment,
                            levels = c("L5_FW", "L5_HW", 
                                       "36h_FWP", "36h_FWM", "36h_FWD", 
                                       "36h_HWP", "36h_HWM", "36h_HWD",
                                       "48h_FWP", "48h_FWM", "48h_FWD", 
                                       "48h_HWP", "48h_HWM", "48h_HWD")))%>%
  mutate(group = case_when(
    grepl("^36h_FW", condition) ~ "36h_FW",
    grepl("^36h_HW", condition) ~ "36h_HW",
    grepl("^48h_FW", condition) ~ "48h_FW",
    grepl("^48h_HW", condition) ~ "48h_HW",
    grepl("^L5_HW", condition)  ~ "L5_HW",
    grepl("^L5_FW", condition)  ~ "L5_FW"
  )) %>%
  mutate(color_col = case_when(
    grepl("^36h_FW", condition) ~ "#4477AA",
    grepl("^36h_HW", condition) ~ "#66CCEE",
    grepl("^48h_FW", condition) ~ "#CCBB44",
    grepl("^48h_HW", condition) ~ "#EE6677"
  ))
#This line then forces it to be plotted in the order you specified before 
counts_long <- counts_long %>%
  mutate(Gene_label = factor(Gene_label, 
                             levels = gene_list_for_count_plot$Gene_label))
plot <- ggplot(counts_long, aes(x = condition, y = gene_counts)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  xlab("Compartment") +
  ylab("Counts") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_summary(fun = mean, geom = "line", lwd = 1.5,
               aes(group = group, color = color_col)) +
  scale_color_identity(guide = "none") +
  facet_wrap(~ Gene_label, scales='free', drop = TRUE, ncol = 4) +
  theme(legend.position = "none") +
  scale_y_continuous(
    breaks = function(x) pretty(x),        
    labels = function(b) format(b, trim=TRUE)  # force all breaks to show more labels 
  )
plot  


# Now generate the heatmap for the all genes differentially expressed in the compartment analysis 
#load necessary packages
library(readxl);library(DESeq2);library(reshape2);library(ggdendro);library(khroma);library(viridis);library(gridExtra);library(cowplot);library(pals)

# perform variance stabilizing transformation 
vsd <- vst(ddsObj)
vst_mat <- assay(vsd)

vst_df <- as.data.frame(vst_mat) %>%
  rownames_to_column("Geneid")

#Set rownames and remove Geneid column before scaling

mat <- vst_df
rownames(mat) <- mat$Geneid
mat$Geneid <- NULL

# Z-score across compartments (i.e., scale by row)
Z <- t(scale(t(as.matrix(mat))))

# Convert to dataframe and keep gene IDs
Z <- as.data.frame(Z) %>%
  rownames_to_column("Geneid")

#filter by significant genes
Z <- Z[Z$Geneid %in% sigGenes_HW_FW_df$Geneid, ]

# add manual annotations 

manual_annotations <- manual_annotations %>%
  select(`Geneid`, Symbol)

Z <- Z %>%
  left_join(manual_annotations, by = c("Geneid" = "Geneid")) %>%
  select(-Geneid)

# Melt to long format for ggplot2

library(reshape2); library(ggdendro)
Z_df <- melt(Z)
Z_df <- na.omit(Z_df)
colnames(Z_df) <- c("Gene", "Sample", "Expression")


# Convert to wide matrix format for clustering

Z_df_matrix <- dcast(Z_df, Gene ~ Sample, value.var = "Expression")
rownames(Z_df_matrix) <- Z_df_matrix$Gene
Z_df_matrix$Gene <- NULL

# Compute distances and gene and sample clusters
distanceGene <- dist(Z_df_matrix)
clusterGene <- hclust(distanceGene, method = "average")
gene_order <- clusterGene$labels[clusterGene$order] 
Z_df$Gene <- factor(Z_df$Gene, levels = gene_order)

# Make dendogram for genes
geneModel <- as.dendrogram(clusterGene)
geneDendrogramData <- segment(dendro_data(geneModel, type = "rectangle"))
geneDendrogram <- ggplot(geneDendrogramData) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_flip() +
  scale_y_reverse(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_dendro()


#Re-factor samples for ggplot2, ordering by compartment and time 
Z_df$Sample <- factor(Z_df$Sample, levels = c("L5_ind1_FW", "L5_ind2_FW", "L5_ind3_FW", "L5_ind1_HW", "L5_ind2_HW", "L5_ind3_HW", "36h_female_ind1_FWP", "36h_female_ind2_FWP", "36h_male_ind3_FWP",  "36h_female_ind2_HWP", "36h_male_ind3_HWP", "36h_female_ind1_FWM", "36h_female_ind2_FWM", "36h_male_ind3_FWM", "36h_female_ind1_HWM", "36h_female_ind2_HWM", "36h_male_ind3_HWM", "36h_female_ind1_FWD", "36h_female_ind2_FWD", "36h_male_ind3_FWD", "36h_female_ind1_HWD", "36h_female_ind2_HWD", "36h_male_ind3_HWD", "48h_male_ind1_FWP", "48h_male_ind5_FWP", "48h_male_ind4_HWP", "48h_male_ind5_HWP", "48h_male_ind1_FWM",  "48h_male_ind4_FWM", "48h_male_ind5_FWM", "48h_male_ind4_HWM", "48h_male_ind5_HWM", "48h_male_ind3_FWD",  "48h_male_ind5_FWD", "48h_male_ind4_HWD", "48h_male_ind5_HWD"))

# Define your color palettes you might use and load libraries
library(pals)
ocean.balance <- ocean.balance
coolwarm <- coolwarm
library(gridExtra); library(patchwork); library(cowplot)

#Construct the heatmap
heatmap <- ggplot(Z_df, aes(x=Sample, y=Gene, fill=Expression)) + geom_raster() + scale_fill_gradientn(colours =coolwarm(256)) + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y = element_blank(), axis.ticks.y=element_blank())
heatmap
grid.arrange(geneDendrogram, heatmap, ncol=1, heights=c(1,5))


# Create heatmap
heatmap <- ggplot(Z_df, aes(x = Sample, y = Gene, fill = Expression)) +
  geom_raster() +
  scale_fill_gradientn(colours = coolwarm(256)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1),
    axis.text.y = element_text(size = 6),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 10, face = "bold"
    ))

heatmap
# Remove duplicated legend
heatmap_clean <- heatmap + theme(legend.position = "none")

# Arrange plots

mid_row <- plot_grid(geneDendrogram, heatmap_clean, ncol = 2, rel_widths = c(0.2, 1))
legend <- get_legend(heatmap)

# Final layout
final_plot <- plot_grid(mid_row, ncol = 1, rel_heights = c(0.2, 1))
plot_grid(final_plot, legend, rel_widths = c(1, 0.12))



manual_annotations_heatmap <- read_excel("working_Ecla_annotation_table_with_flybase_names.xlsx")

manual_annotations_heatmap  <- manual_annotations_heatmap  %>%
  select(Symbol, gene_label)

Z_df <- Z_df %>%
  left_join(manual_annotations_heatmap, by = c("Gene" = "Symbol"))

# Re-enforce clustering order
Z_df$Gene <- factor(Z_df$Gene, levels = gene_order)


# Generate label mapping safely (including NAs)
gene_label_map <- setNames(Z_df$gene_label, as.character(Z_df$Gene)) 

# Plot
heatmap <- ggplot(Z_df, aes(x = Sample, y = Gene, fill = Expression)) +
  geom_raster() +
  scale_fill_gradientn(colours = coolwarm(256)) +
  scale_y_discrete(labels = gene_label_map) +  # override tick labels
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1),
    axis.text.y = element_text(size = 6),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 10, face = "bold")
  )

heatmap

heatmap <- ggplot(Z_df, aes(x=Sample, y=Gene, fill=Expression)) + geom_raster() + scale_fill_gradientn(colours =coolwarm(256)) + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y = element_blank(), axis.ticks.y=element_blank())
heatmap
grid.arrange(geneDendrogram, heatmap, ncol=1, heights=c(1,5))


# Create heatmap
heatmap <- ggplot(Z_df, aes(x = Sample, y = Gene, fill = Expression)) +
  geom_raster() +
  scale_fill_gradientn(colours = coolwarm(256)) +
  scale_y_discrete(labels = gene_label_map) +  # override tick labels
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1),
    axis.text.y = element_text(size = 6),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 10, face = "bold"
    ))

heatmap
# Remove duplicated legend
heatmap_clean <- heatmap + theme(legend.position = "none")

# Arrange plots

mid_row <- plot_grid(geneDendrogram, heatmap_clean, ncol = 2, rel_widths = c(0.2, 1))
legend <- get_legend(heatmap)

# Final layout
final_plot <- plot_grid(mid_row, ncol = 1, rel_heights = c(0.2, 1))
plot_grid(final_plot, legend, rel_widths = c(1, 0.12))



manual_annotations <- read_excel("working_Ecla_annotation_table_with_flybase_names.xlsx")

Z_df<- as.data.frame(Z_df_matrix) %>%
  rownames_to_column("Symbol")%>%
  left_join(manual_annotations, by = c("Symbol" = "Symbol")) %>%
  relocate("Symbol", .after = "Geneid") %>%
  relocate("Name", .after = "Symbol") %>%
  relocate("FB_gene_symbol", .after = "Name") %>%
  relocate("gene_fullname", .after = "FB_gene_symbol")

write.table(Z_df, file="heatmap_HW_FW_DEgenes_36hr_all_without48hcontamind.tsv", quote=F, sep="\t",row.names=FALSE, na="")

#sample PCA 
# make a PCA
vsd <- vst(ddsObj, blind=FALSE)
pcaData<- plotPCA(vsd, intgroup=c("dev_stage", "compartment"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=compartment, shape=dev_stage)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() 


pcaData<- plotPCA(vsd, intgroup=c("sample", "compartment"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))


ggplot(pcaData, aes(PC1, PC2, color=sample, shape=compartment)) +
  scale_color_brewer(palette = "Paired") +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


###############################################
# GO Enrichment Analysis for Compartments 
# inspired by tutorial here https://www.nature.com/articles/s41596-024-01020-z#Sec43

library(clusterProfiler)
library(BaseSet)

gaf_data <- getGAF("GCF_041222505.1-RS_2025_04_slim_GO.mapped.gaf")
gaf_df <- as.data.frame(gaf_data)


gaf_df <- gaf_df %>%
  relocate("DB_Object_ID", .after = "elements") %>%
  relocate("elements", .after = "sets") %>%
  relocate("sets", .after = "DB_Object_ID")
  

DE_genes <- sigGenes_HW_FW_df %>%
  select(Geneid)

manual_annotations <- read_excel("working_Ecla_annotation_table_with_flybase_names.xlsx")

Geneid_to_ID <- manual_annotations %>%
  select("Geneid", `old Gene ID`)

DE_genes <- DE_genes %>%
  left_join(Geneid_to_ID, by = c("Geneid" = "Geneid")) %>%
  select(-Geneid)

enrich_result <- compareCluster(DE_genes,
                                          fun = "enricher", TERM2GENE = gaf_df[, c(2, 1)],
)




#install.packages("ontologyIndex")
library(ontologyIndex)
download.file("http://purl.obolibrary.org/obo/go.obo", destfile = "go.obo", mode = "wb")
go <- get_ontology("go.obo", extract_tags = "everything")

go$id[1:5]         # GO IDs
go$name[1:5]       # GO descriptions

go_df <- data.frame(
  ID   = go$id,
  Term = unname(go$name[go$id]),   # use go$id to order the terms to match IDs
  stringsAsFactors = FALSE
)



GO_plot <- dotplot(enrich_result, by = "Count",
                      showCategory = 7,
                      label_format = 40) +
  theme_minimal() +
  theme(axis.text.x = element_text(vjust = 1, hjust = 1,
                                   angle = 30, size = 10)) +
  xlab(NULL) + ggtitle(NULL)

GO_plot

# add term 
go_map <- setNames(go_df$Term, go_df$ID)

df <- as.data.frame(enrich_result)

custom_order <- c("GO:0003700", "GO:0038023", "GO:0008092", "GO:0030154", "GO:0032502")


df <- df %>%
  # Make Description a factor with levels ordered by custom_order of IDs
  mutate(Description = factor(Description, levels = Description[match(custom_order, ID)])) %>%
  arrange(factor(ID, levels = custom_order))
         
# Now plot
ggplot(df[1:5, ], aes(x = FoldEnrichment, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_continuous(low = "red", high = "blue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
    xlab("Fold Enrichment") +
  scale_y_discrete( labels = function(ids) {
    labs <- go_map[ids]
    labs[is.na(labs)] <- ids[is.na(labs)] # fallback to ID if missing
    labs
  }) +
  ylab(NULL)


# Now do the GO enrichment analysis with the whole HW vs FW results at 36hr 

sigGenes_whole_HW_vs_FW<- allGenesHW_vs_FW_36h[!is.na(allGenesHW_vs_FW_36h$padj) & allGenesHW_vs_FW_36h$padj < 0.05, ]


#######################

# GO Analysis 


library(clusterProfiler)

#install.packages("BaseSet")
library(BaseSet)

gaf_data <- getGAF("GCF_041222505.1-RS_2025_04_slim_GO.mapped.gaf")
gaf_df <- as.data.frame(gaf_data)


gaf_df <- gaf_df %>%
  relocate("DB_Object_ID", .after = "elements") %>%
  relocate("elements", .after = "sets") %>%
  relocate("sets", .after = "DB_Object_ID")


DE_genes_HW_FW <- sigGenes_whole_HW_vs_FW %>%
  select(Geneid)

manual_annotations <- read_excel("working_Ecla_annotation_table_with_flybase_names.xlsx")

Geneid_to_ID <- manual_annotations %>%
  select("Geneid", `old Gene ID`)

DE_genes_HW_FW <- DE_genes_HW_FW %>%
  left_join(Geneid_to_ID, by = c("Geneid" = "Geneid")) %>%
  select(-Geneid)

enrich_result_whole_HW_FW <- compareCluster(DE_genes_HW_FW,
                                fun = "enricher", TERM2GENE = gaf_df[, c(2, 1)],
)


#install.packages("ontologyIndex")
library(ontologyIndex)
download.file("http://purl.obolibrary.org/obo/go.obo", destfile = "go.obo", mode = "wb")
go <- get_ontology("go.obo", extract_tags = "everything")

go$id[1:5]         # GO IDs
go$name[1:5]       # GO descriptions

go_df <- data.frame(
  ID   = go$id,
  Term = unname(go$name[go$id]),   # use go$id to order the terms to match IDs
  stringsAsFactors = FALSE
)



GO_plot <- dotplot(enrich_result_whole_HW_FW, by = "Count",
                   showCategory = 6,
                   label_format = 40) +
  theme_minimal() +
  theme(axis.text.x = element_text(vjust = 1, hjust = 1,
                                   angle = 30, size = 10)) +
  xlab(NULL) + ggtitle(NULL)

GO_plot

# add term 
go_map <- setNames(go_df$Term, go_df$ID)

custom_order <- c("GO:0045202", "GO:0005773", "GO:0008289", "GO:0006629","GO:0005783","GO:0032502")

df <- df %>%
  # Make Description a factor with levels ordered by custom_order of IDs
 mutate(Description = factor(Description, levels = Description[match(custom_order, ID)])) %>%
arrange(factor(ID, levels = custom_order))

# Now plot
ggplot(df[1:6, ], aes(x = FoldEnrichment, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_continuous(low = "red", high = "blue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  xlab("Fold Enrichment") +
  scale_y_discrete( labels = function(ids) {
    labs <- go_map[ids]
    labs[is.na(labs)] <- ids[is.na(labs)] # fallback to ID if missing
    labs
  }) +
  ylab(NULL)


# Now make heatmaps with only certain GO terms from the DEGs from compartment analysis

# first transcription factors 

sigGenes_HW_FW_GO_TF <-  sigGenes_HW_FW_df %>%
  left_join(Geneid_to_ID, by = c("Geneid" = "Geneid")) %>%
  left_join(gaf_df, by = c(`old Gene ID` = "DB_Object_ID"))

sigGenes_HW_FW_GO_TF <- sigGenes_HW_FW_GO_TF %>%
  filter(sets %in% c("GO:0003700", "GO:0008134")) %>%
  select(Geneid)
  
dds <- DESeqDataSetFromMatrix(
  countData = countdata_all,
  colData = sampleinfo_all,
  design = ~ compartment
)

ddsObj.raw <- DESeq(dds)
ddsObj <- estimateSizeFactors(ddsObj.raw)
colData(ddsObj.raw)
colData(ddsObj)
ddsObj <- DESeq(ddsObj.raw)
res <- results(ddsObj, alpha=0.05) #adding p-value cutoff
res

resultsNames(ddsObj) 

ddsObj$compartment <- relevel(ddsObj$compartment, ref = "36h_HWM")

ddsObj <- DESeq(ddsObj)
resultsNames(ddsObj) 

# perform variance stabilizing transformation
vsd <- vst(ddsObj)
vst_mat <- assay(vsd)

vst_df <- as.data.frame(vst_mat) %>%
  rownames_to_column("Geneid")

#Set rownames and remove Geneid column before scaling

mat <- vst_df
rownames(mat) <- mat$Geneid
mat$Geneid <- NULL

# Z-score across compartments (i.e., scale by row)
Z <- t(scale(t(as.matrix(mat))))

# Convert to dataframe and keep gene IDs
Z <- as.data.frame(Z) %>%
  rownames_to_column("Geneid")

#filter by significant genes
Z <- Z[Z$Geneid %in% sigGenes_HW_FW_GO_TF$Geneid, ]

# add manual annotations 

manual_annotations <- manual_annotations %>%
  select(`Geneid`, Symbol)

Z <- Z %>%
  left_join(manual_annotations, by = c("Geneid" = "Geneid")) %>%
  select(-Geneid)

# Melt to long format for ggplot2

library(reshape2); library(ggdendro)
Z_df <- melt(Z)
Z_df <- na.omit(Z_df)
colnames(Z_df) <- c("Gene", "Sample", "Expression")


# Convert to wide matrix format for clustering

Z_df_matrix <- dcast(Z_df, Gene ~ Sample, value.var = "Expression")
rownames(Z_df_matrix) <- Z_df_matrix$Gene
Z_df_matrix$Gene <- NULL

# Compute distances and gene and sample clusters
distanceGene <- dist(Z_df_matrix)
clusterGene <- hclust(distanceGene, method = "average")
gene_order <- clusterGene$labels[clusterGene$order] 
Z_df$Gene <- factor(Z_df$Gene, levels = gene_order)

# Make dendogram for genes
geneModel <- as.dendrogram(clusterGene)
geneDendrogramData <- segment(dendro_data(geneModel, type = "rectangle"))
geneDendrogram <- ggplot(geneDendrogramData) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_flip() +
  scale_y_reverse(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_dendro()


#Re-factor samples for ggplot2, ordering by compartment and time 
Z_df$Sample <- factor(Z_df$Sample, levels = c("L5_ind1_FW", "L5_ind2_FW", "L5_ind3_FW", "L5_ind1_HW", "L5_ind2_HW", "L5_ind3_HW", "36h_female_ind1_FWP", "36h_female_ind2_FWP", "36h_male_ind3_FWP",  "36h_female_ind2_HWP", "36h_male_ind3_HWP", "36h_female_ind1_FWM", "36h_female_ind2_FWM", "36h_male_ind3_FWM", "36h_female_ind1_HWM", "36h_female_ind2_HWM", "36h_male_ind3_HWM", "36h_female_ind1_FWD", "36h_female_ind2_FWD", "36h_male_ind3_FWD", "36h_female_ind1_HWD", "36h_female_ind2_HWD", "36h_male_ind3_HWD", "48h_male_ind1_FWP", "48h_male_ind5_FWP", "48h_male_ind4_HWP", "48h_male_ind5_HWP", "48h_male_ind1_FWM",  "48h_male_ind4_FWM", "48h_male_ind5_FWM", "48h_male_ind4_HWM", "48h_male_ind5_HWM", "48h_male_ind3_FWD",  "48h_male_ind5_FWD", "48h_male_ind4_HWD", "48h_male_ind5_HWD"))

# Define your color palettes you might use and load libraries
library(pals)
ocean.balance <- ocean.balance
coolwarm <- coolwarm
library(gridExtra); library(patchwork); library(cowplot)

#Construct the heatmap
heatmap <- ggplot(Z_df, aes(x=Sample, y=Gene, fill=Expression)) + geom_raster() + scale_fill_gradientn(colours =coolwarm(256)) + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y = element_blank(), axis.ticks.y=element_blank())
heatmap
grid.arrange(geneDendrogram, heatmap, ncol=1, heights=c(1,5))


# Create heatmap
heatmap <- ggplot(Z_df, aes(x = Sample, y = Gene, fill = Expression)) +
  geom_raster() +
  scale_fill_gradientn(colours = coolwarm(256)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1),
    axis.text.y = element_text(size = 6),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 10, face = "bold"
    ))

heatmap
# Remove duplicated legend
heatmap_clean <- heatmap + theme(legend.position = "none")

# Arrange plots

mid_row <- plot_grid(geneDendrogram, heatmap_clean, ncol = 2, rel_widths = c(0.2, 1))
legend <- get_legend(heatmap)

# Final layout
final_plot <- plot_grid(mid_row, ncol = 1, rel_heights = c(0.2, 1))
plot_grid(final_plot, legend, rel_widths = c(1, 0.12))




manual_annotations <- read_excel("working_Ecla_annotation_table_with_flybase_names.xlsx")

Z_df<- as.data.frame(Z_df_matrix) %>%
  rownames_to_column("Symbol")%>%
  left_join(manual_annotations, by = c("Symbol" = "Symbol")) %>%
  relocate("Symbol", .after = "Geneid") %>%
  relocate("Name", .after = "Symbol") %>%
  relocate("FB_gene_symbol", .after = "Name") %>%
  relocate("gene_fullname", .after = "FB_gene_symbol")

# now a heatmap for signaling related genes 

sigGenes_HW_FW_GO_signaling <-  sigGenes_HW_FW_df %>%
  left_join(Geneid_to_ID, by = c("Geneid" = "Geneid")) %>%
  left_join(gaf_df, by = c(`old Gene ID` = "DB_Object_ID"))

sigGenes_HW_FW_GO_signaling <- sigGenes_HW_FW_GO_signaling %>%
  filter(sets %in% c("GO:0023052", "GO:0038023", "GO:0005102")) %>%
  select(Geneid)

mat <- vst_df
rownames(mat) <- mat$Geneid
mat$Geneid <- NULL

# Z-score across compartments (i.e., scale by row)
Z <- t(scale(t(as.matrix(mat))))

# Convert to dataframe and keep gene IDs
Z <- as.data.frame(Z) %>%
  rownames_to_column("Geneid")

#filter by significant genes
Z <- Z[Z$Geneid %in% sigGenes_HW_FW_GO_signaling$Geneid, ]

# add manual annotations 

manual_annotations <- manual_annotations %>%
  select(`Geneid`, Symbol)

Z <- Z %>%
  left_join(manual_annotations, by = c("Geneid" = "Geneid")) %>%
  select(-Geneid)

# Melt to long format for ggplot2

library(reshape2); library(ggdendro)
Z_df <- melt(Z)
Z_df <- na.omit(Z_df)
colnames(Z_df) <- c("Gene", "Sample", "Expression")


# Convert to wide matrix format for clustering

Z_df_matrix <- dcast(Z_df, Gene ~ Sample, value.var = "Expression")
rownames(Z_df_matrix) <- Z_df_matrix$Gene
Z_df_matrix$Gene <- NULL

# Compute distances and gene and sample clusters
distanceGene <- dist(Z_df_matrix)
clusterGene <- hclust(distanceGene, method = "average")
gene_order <- clusterGene$labels[clusterGene$order] 
Z_df$Gene <- factor(Z_df$Gene, levels = gene_order)

# Make dendogram for genes
geneModel <- as.dendrogram(clusterGene)
geneDendrogramData <- segment(dendro_data(geneModel, type = "rectangle"))
geneDendrogram <- ggplot(geneDendrogramData) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_flip() +
  scale_y_reverse(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_dendro()


#Re-factor samples for ggplot2, ordering by compartment and time 
Z_df$Sample <- factor(Z_df$Sample, levels = c("L5_ind1_FW", "L5_ind2_FW", "L5_ind3_FW", "L5_ind1_HW", "L5_ind2_HW", "L5_ind3_HW", "36h_female_ind1_FWP", "36h_female_ind2_FWP", "36h_male_ind3_FWP",  "36h_female_ind2_HWP", "36h_male_ind3_HWP", "36h_female_ind1_FWM", "36h_female_ind2_FWM", "36h_male_ind3_FWM", "36h_female_ind1_HWM", "36h_female_ind2_HWM", "36h_male_ind3_HWM", "36h_female_ind1_FWD", "36h_female_ind2_FWD", "36h_male_ind3_FWD", "36h_female_ind1_HWD", "36h_female_ind2_HWD", "36h_male_ind3_HWD", "48h_male_ind1_FWP", "48h_male_ind5_FWP", "48h_male_ind4_HWP", "48h_male_ind5_HWP", "48h_male_ind1_FWM",  "48h_male_ind4_FWM", "48h_male_ind5_FWM", "48h_male_ind4_HWM", "48h_male_ind5_HWM", "48h_male_ind3_FWD",  "48h_male_ind5_FWD", "48h_male_ind4_HWD", "48h_male_ind5_HWD"))

# Define your color palettes you might use and load libraries
library(pals)
ocean.balance <- ocean.balance
coolwarm <- coolwarm
library(gridExtra); library(patchwork); library(cowplot)

#Construct the heatmap
heatmap <- ggplot(Z_df, aes(x=Sample, y=Gene, fill=Expression)) + geom_raster() + scale_fill_gradientn(colours =coolwarm(256)) + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y = element_blank(), axis.ticks.y=element_blank())
heatmap
grid.arrange(geneDendrogram, heatmap, ncol=1, heights=c(1,5))


# Create heatmap
heatmap <- ggplot(Z_df, aes(x = Sample, y = Gene, fill = Expression)) +
  geom_raster() +
  scale_fill_gradientn(colours = coolwarm(256)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1),
    axis.text.y = element_text(size = 6),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 10, face = "bold"
    ))

heatmap
# Remove duplicated legend
heatmap_clean <- heatmap + theme(legend.position = "none")

# Arrange plots

mid_row <- plot_grid(geneDendrogram, heatmap_clean, ncol = 2, rel_widths = c(0.2, 1))
legend <- get_legend(heatmap)

# Final layout
final_plot <- plot_grid(mid_row, ncol = 1, rel_heights = c(0.2, 1))
plot_grid(final_plot, legend, rel_widths = c(1, 0.12))


manual_annotations_heatmap <- read_excel("working_Ecla_annotation_table_with_flybase_names.xlsx")

manual_annotations_heatmap  <- manual_annotations_heatmap  %>%
  select(Symbol, gene_label)

Z_df <- Z_df %>%
  left_join(manual_annotations_heatmap, by = c("Gene" = "Symbol"))

# Re-enforce clustering order
Z_df$Gene <- factor(Z_df$Gene, levels = gene_order)


# Generate label mapping safely (including NAs)
gene_label_map <- setNames(Z_df$gene_label, as.character(Z_df$Gene)) 

# Plot
heatmap <- ggplot(Z_df, aes(x = Sample, y = Gene, fill = Expression)) +
  geom_raster() +
  scale_fill_gradientn(colours = coolwarm(256)) +
  scale_y_discrete(labels = gene_label_map) +  # override tick labels
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1),
    axis.text.y = element_text(size = 6),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 10, face = "bold")
  )

heatmap

heatmap <- ggplot(Z_df, aes(x=Sample, y=Gene, fill=Expression)) + geom_raster() + scale_fill_gradientn(colours =coolwarm(256)) + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y = element_blank(), axis.ticks.y=element_blank())
heatmap
grid.arrange(geneDendrogram, heatmap, ncol=1, heights=c(1,5))


# Create heatmap
heatmap <- ggplot(Z_df, aes(x = Sample, y = Gene, fill = Expression)) +
  geom_raster() +
  scale_fill_gradientn(colours = coolwarm(256)) +
  scale_y_discrete(labels = gene_label_map) +  # override tick labels
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1),
    axis.text.y = element_text(size = 6),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 10, face = "bold"
    ))

heatmap
# Remove duplicated legend
heatmap_clean <- heatmap + theme(legend.position = "none")

# Arrange plots

mid_row <- plot_grid(geneDendrogram, heatmap_clean, ncol = 2, rel_widths = c(0.2, 1))
legend <- get_legend(heatmap)

# Final layout
final_plot <- plot_grid(mid_row, ncol = 1, rel_heights = c(0.2, 1))
plot_grid(final_plot, legend, rel_widths = c(1, 0.12))



manual_annotations <- read_excel("working_Ecla_annotation_table_with_flybase_names.xlsx")

Z_df<- as.data.frame(Z_df_matrix) %>%
  rownames_to_column("Symbol")%>%
  left_join(manual_annotations, by = c("Symbol" = "Symbol")) %>%
  relocate("Symbol", .after = "Geneid") %>%
  relocate("Name", .after = "Symbol") %>%
  relocate("FB_gene_symbol", .after = "Name") %>%
  relocate("gene_fullname", .after = "FB_gene_symbol")

################################
# Now run the enrichment analysis for the whole HW and whole FW DEGs seperately with Drosophila slim terms


# GO Enrichment Analysis for Compartments 
# inspired by tutorial here https://www.nature.com/articles/s41596-024-01020-z#Sec43

library(clusterProfiler)

#install.packages("BaseSet")
library(BaseSet)

dmel_gaf_data <- getGAF("GCF_041222505.1-RS_2025_04_drosophila_slim_GO.mapped.gaf")
dmel_gaf_df <- as.data.frame(dmel_gaf_data)


dmel_gaf_df <- dmel_gaf_df %>%
  relocate("DB_Object_ID", .after = "elements") %>%
  relocate("elements", .after = "sets") %>%
  relocate("sets", .after = "DB_Object_ID")


DE_genes_whole_HW <- sigGenes_whole_HW_vs_FW %>%
  filter(direction %in% c("down")) %>%
  select(Geneid)

DE_genes_whole_FW <- sigGenes_whole_HW_vs_FW %>%
  filter(direction %in% c("up")) %>%
  select(Geneid)

manual_annotations <- read_excel("working_Ecla_annotation_table_with_flybase_names.xlsx")

Geneid_to_ID <- manual_annotations %>%
  select("Geneid", `old Gene ID`)

DE_genes_whole_FW <- DE_genes_whole_FW %>%
  left_join(Geneid_to_ID, by = c("Geneid" = "Geneid")) %>%
  select(-Geneid)

DE_genes_whole_HW <- DE_genes_whole_HW %>%
left_join(Geneid_to_ID, by = c("Geneid" = "Geneid")) %>%
  select(-Geneid)

enrich_result_whole_FW <- compareCluster(DE_genes_whole_FW,
                                fun = "enricher", TERM2GENE = dmel_gaf_df[, c(2, 1)],)

                                            
enrich_result_whole_HW <- compareCluster(DE_genes_whole_HW,
                                         fun = "enricher", TERM2GENE = dmel_gaf_df[, c(2, 1)],)


#install.packages("ontologyIndex")
library(ontologyIndex)
download.file("http://purl.obolibrary.org/obo/go.obo", destfile = "go.obo", mode = "wb")
go <- get_ontology("go.obo", extract_tags = "everything")

go$id[1:5]         # GO IDs
go$name[1:5]       # GO descriptions

go_df <- data.frame(
  ID   = go$id,
  Term = unname(go$name[go$id]),   # use go$id to order the terms to match IDs
  stringsAsFactors = FALSE
)

# Now visualize FW whole genes

GO_plot <- dotplot(enrich_result_whole_FW, by = "Count",
                   showCategory = 9,
                   label_format = 40) +
  theme_minimal() +
  theme(axis.text.x = element_text(vjust = 1, hjust = 1,
                                   angle = 30, size = 10)) +
  xlab(NULL) + ggtitle(NULL)

GO_plot

# add term 
go_map <- setNames(go_df$Term, go_df$ID)

GO_plot <- dotplot(enrich_result_whole_FW, by = "Count",
                   showCategory = 9,
                   label_format = 40) +
  theme_minimal() +
  theme(axis.text.x = element_text(vjust = 1, hjust = 1,
                                   angle = 30, size = 10)) +
  scale_y_discrete(labels = function(ids) {
    labs <- go_map[ids]
    labs[is.na(labs)] <- ids[is.na(labs)] 
    labs
  }) +
  ggtitle(NULL) + xlab("GeneRatio") 

GO_plot

#Gene Ratio x axis 

df <- as.data.frame(enrich_result_whole_FW)


custom_order <- c("GO:0042302", "GO:0008061", "GO:0016746", "GO:0005777", "GO:0031012", "GO:0005615", "GO:0006629", "GO:0030198", "GO:0030705", "GO:0044782")


df <- df %>%
  # Make Description a factor with levels ordered by custom_order of IDs
  mutate(Description = factor(Description, levels = Description[match(custom_order, ID)])) %>%
  arrange(factor(ID, levels = custom_order))

# Now plot
ggplot(df[1:9, ], aes(x = FoldEnrichment, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_continuous(low = "red", high = "blue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  xlab("Fold Enrichment") +
  scale_y_discrete( labels = function(ids) {
    labs <- go_map[ids]
    labs[is.na(labs)] <- ids[is.na(labs)] 
    labs
  }) +
  ylab(NULL)



# Now visualize HW whole genes

GO_plot <- dotplot(enrich_result_whole_HW, by = "Count",
                   showCategory = 17,
                   label_format = 40) +
  theme_minimal() +
  theme(axis.text.x = element_text(vjust = 1, hjust = 1,
                                   angle = 30, size = 10)) +
  xlab(NULL) + ggtitle(NULL)

GO_plot

# add term 
go_map <- setNames(go_df$Term, go_df$ID)

GO_plot <- dotplot(enrich_result_whole_HW, by = "Count",
                   showCategory = 17,
                   label_format = 40) +
  theme_minimal() +
  theme(axis.text.x = element_text(vjust = 1, hjust = 1,
                                   angle = 30, size = 10)) +
  scale_y_discrete(labels = function(ids) {
    labs <- go_map[ids]
    labs[is.na(labs)] <- ids[is.na(labs)] 
    labs
  }) +
  ggtitle(NULL) + xlab("GeneRatio") 

GO_plot

#Gene Ratio x axis 

df <- as.data.frame(enrich_result_whole_HW)


custom_order <- c("GO:0042302", "GO:0008061", "GO:0016746", "GO:0005777", "GO:0005615", "GO:0006629", "GO:0030198", "GO:0030705", "GO:0044782")


df <- df %>%
# Make Description a factor with levels ordered by custom_order of IDs
  mutate(Description = factor(Description, levels = Description[match(custom_order, ID)])) %>%
  arrange(factor(ID, levels = custom_order))

# Now plot
ggplot(df[1:17, ], aes(x = FoldEnrichment, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_continuous(low = "red", high = "blue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  xlab("Fold Enrichment") +
  scale_y_discrete( labels = function(ids) {
    labs <- go_map[ids]
    labs[is.na(labs)] <- ids[is.na(labs)] 
    labs
  }) +
  ylab(NULL)

