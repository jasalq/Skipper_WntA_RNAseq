# Skipper_WntA_RNAseq
<body>  
This repository contains all code associated with the differential expression analyses performed in the manuscript Alqassar et al. 2025 <a href=""> WntA expression and wing transcriptomics illuminate the evolution of stripe patterns in skipper butterflies </a>. The code found in this repository processes, maps, and performs differential expression analyses on RNA sequencing data generated during the study where we aimed to understand differences in gene expression between wing tissues in skipper butterflies.		
<body/>	
<br><br>	
<strong>Contact info</strong>: Jasmine Alqassar (j.alqassar@gwu.edu) and Arnaud Martin (arnaud@gwu.edu)
<br><br>
<strong>Data Availability</strong>: All sequencing data generated during this study has been deposited in the NCBI Sequence Read Archive under BioProject PRJNA660444. 

## Pipeline Outline 
1. [RNAseq Quality Assessment](#assessment-of-rna-sequencing-quality-using-fastqc) 
2.[Adding Drosophila Homology Evidence to the Genome Annotation](#homology-searches-using-reciprocal-blastp-to-flybase)
3.[Read Mapping to the Reference Genome](#star-mapping-of-rnaseq-data-to-the-e-clarus-reference-genome-gcf_0412225051)
4.[Read Counting](#read-counts-with-featurecounts)
5.[Differential Expression Analysis DESeq2](#differential-expression-analysis-with-deseq2-in-rstudio)
6.[Heatmap and Count Plots Data Visualization](#differential-expression-analysis-with-deseq2-in-rstudio)
7.[GO Enrichment Analysis with GO Subsets](#go-enrichment-analysis-with-go-subsets)


### Tools Used 
* Blast+ v.2.16.0+
* Fastp v.0.21.0	
* FastQC v.0.11.8
* JDK v.21.0.1
* NCBI Datasets
* OWLtools
* RStudio v2024.04.2+764
* R/ BaseSet
* R/ clusterProfiler
* R/ cowplot
* R/ DESeq2
* R/ dplyr
* R/ ggdendro
* R/ ggplot2
* R/ ggrepel
* R/ gridExtra
* R/ ontologyIndex
* R/ pals
* R/ patchwork
* R/ readxl
* R/ reshape2
* R/ stringr
* R/ tidyverse
* STAR v.2.7.11b
* Subread v.2.0.8

## Assessment of RNA Sequencing Quality using FastQC
Quality of RNA sequencing was accessed using FastQC
```
#!/bin/sh
#SBATCH -J SSS_wing_RNAseq_fastqc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=j.alqassar@gwu.edu
#SBATCH -o SSS_wing_RNAseq_fastqc.out #redirecting stdout
#SBATCH -p defq #queue 
#SBATCH -N 1 #amount of nodes 
#SBATCH -n 16 #amount of cores 
#SBATCH -t 24:00:00


echo "=========================================================="
echo "Running on node : $HOSTNAME"
echo "Current directory : $PWD"
echo "Current job ID : $SLURM_JOB_ID"
echo "=========================================================="

cd /scratch/martinlab/jasmine/Skipper_data

module load jdk/21.0.1
module load fastQC/0.11.8

mkdir /scratch/martinlab/jasmine/Skipper_data/fastqc_output

for i in *fastq.gz; do
        fastqc -f fastq -t 24 -o /scratch/martinlab/jasmine/Skipper_data/fastqc_output $i;
        done
```

## Homology Searches Using Reciprocal BLASTp to FlyBase
**use the protein.faa file from the NCBI download to blast to flybase**
**find this file -> dmel-all-translation.fasta from the latest FlyBase release and unzip**

```
cd db
module load blast+/2.16.0+
makeblastdb -in dmelAA -input_type fasta -dbtype prot -title dmelAA
```
**Now run BLASTp on the HPC** 
```
#!/bin/sh
#SBATCH -J flybase_blastp_Ecla_genome
#SBATCH --mail-type=END
#SBATCH --mail-user=j.alqassar@gwu.edu
#SBATCH -o flybase_blastp_Ecla_genome.out #redirecting stdout
#SBATCH -p nano #queue 
#SBATCH -n 40 #amount of cores 
#SBATCH -t 00:30:00

echo "=========================================================="
echo "Running on node : $HOSTNAME"
echo "Current directory : $PWD"
echo "Current job ID : $SLURM_JOB_ID"
echo "Job Started:"
date
echo "=========================================================="

module load blast+/2.16.0+

DB=/gpfs/automountdir/gpfs/scratch/martinlab/jasmine/skipper_diff_exp/get_flybase_names/db/dmelAA
QUERY=/scratch/martinlab/jasmine/skipper_diff_exp/get_flybase_names/Ecla_GCF_041222505.1_protein.faa
PREFIX=all_Ecla_proteins_for_flybase

blastp -db ${DB} \
-query ${QUERY} \
-out ${PREFIX}.out \
-outfmt 6 \
-evalue 1e-5 \
-num_threads 40

echo "=========================================================="
echo "Job Finished  $SLURM_JOB_ID:"
date
echo "=========================================================="
```
**Now back to RStudio to add the results to your data matrix**	 

Set working directory
```
setwd("/Users/jasminealqassar/Library/CloudStorage/GoogleDrive-j.alqassar@gwmail.gwu.edu/.shortcut-targets-by-id/1Hi7WIp_ha7vnyQUclvt-AJRl6AF8y1hj/Martin Lab/@LabData/2025_Skipper_WntA/SSS_Differential_Expression/R_working_dir")
```
Install and load necessary packages
```
install.packages("tidyverse"); install.packages("stringr"); install.packages("dplyr"); install.packages("DESeq2")
library(tidyverse); library(stringr); library(dplyr); library(DESeq2); library(ggplot2); library(readxl)
```
Read in the results
```
flybase_results <- read_tsv("all_Ecla_proteins_for_flybase.out", col_names=FALSE)
```
Download and read in the following tables from the latest [FlyBase genome release GUI page ](#https://flybase.org/downloads/bulkdata) and clean up column names before reading in
```
  flybase_prot_to_Symbol <- read_tsv("dmel_unique_protein_isoforms_fb_2025_02.tsv", col_names=TRUE, comment="#")
  flybase_gn_to_bpp <- read_tsv("fbgn_fbtr_fbpp_expanded_fb_2025_02.tsv", col_names=TRUE, comment="#") 
  flybase_gn_summary <- read_tsv("automated_gene_summaries_fb_2025_02.tsv", col_names=TRUE, comment="#") 
```
Now organize and add the annotations to your gene annotation table
```
  flybase_prot_to_Symbol<- flybase_prot_to_Symbol %>%
    select(FBgn, FB_gene_symbol) %>%
    distinct(FBgn, .keep_all = TRUE) %>%
    left_join(flybase_gn_summary, by = c("FBgn" = "FBgn"))
  
  flybase_gn_to_bpp <- flybase_gn_to_bpp %>%
    select(gene_ID, gene_fullname, polypeptide_ID) %>%
    distinct(polypeptide_ID, .keep_all = TRUE) %>%
    rename("gene_ID" = "FlyBase_FBgn") %>%
    drop_na()
  
  flybase_results <- flybase_results %>%
    select(X1,X2,X11)
  
  colnames(flybase_results)  <- c("Protein", "Flybase_prot", "Eval")
  
Geneid_XP <- read_excel("Ecla_gene_annotation_table_orig.xlsx") #gene annotation list 
  Geneid_XP <- Geneid_XP %>%
    select('Geneid', 'Protein accession') %>%
    drop_na()
  
  flybase_results <- flybase_results %>%
    left_join(Geneid_XP, by = c("Protein" = "Protein accession")) %>%
    relocate("Geneid", .after = "Protein") %>%
    group_by(Geneid) %>%
    slice_min(order_by = .data$Eval, with_ties = FALSE) %>%
    ungroup() %>%
    left_join(flybase_gn_to_bpp, by = c("Flybase_prot" = "polypeptide_ID")) %>%
    left_join(flybase_prot_to_Symbol, by = c("FlyBase_FBgn" = "FBgn")) 
  
  flybase_results_clean <- flybase_results %>%
    select(-Protein)
  
Ecla_annotation_table <- read_excel("Ecla_gene_annotation_table_orig.xlsx", col_names=TRUE) 

Ecla_annotation_table  <- Ecla_annotation_table %>%
  left_join(flybase_results_clean, by = c("Geneid" = "Geneid"))

write.table(Ecla_annotation_table, file="Ecla_annotation_table_with_flybase_names.tsv", quote=F, sep="\t",row.names=FALSE, na="")
```

## STAR Mapping of RNAseq data to the <em>E. clarus </em> reference genome (GCF_041222505.1)

**Download <em>E. clarus</em> RefSeq genome and annotation from NCBI using NCBI Datasets tool** 
```
mamba activate ncbi_datsets
datasets download genome accession GCF_041222505.1 --include gtf,rna,cds,protein,genome,seq-report
unzip ncbi_dataset
mv /scratch/martinlab/jasmine/skipper_diff_exp/ncbi_skipper_genome/ncbi_dataset/data/GCF_041222505.1 ../../
```
**In addition download the list of all annotated genes by clicking "view annotated genes" on the genome release page selecting all and hit download table while clicking "one sequence per gene" option**
I manually edited the GeneID column to have a LOC prefix	
**Download <em>E. clarus</em> RefSeq genome from NCBI using NCBI Datasets tool** 
```
mamba activate NCBI_datsets
datasets download genome accession GCF_041222505.1
unzip ncbi_dataset
mv ncbi_dataset/data/GCF_041222505.1/GCF_041222505.1_WU_Ecla_fem_2.2_genomic.fna ../../../
```
**Create a STAR genome index for the <em>E. clarus</em> genome**  
Prior to running, to determine the correct value for --genomeSAindexNbases I used the formula provided where the genome is 451.8 Mb: min(14, log2(451,800,000) / 2 - 1).

Also activate mamba enviroment with STAR installed before submitting 
```
mamba activate star_2.7.11b_JDA
```

```
#!/bin/sh
#SBATCH -J skipper_genome_index
#SBATCH --mail-type=ALL
#SBATCH --mail-user=j.alqassar@gwu.edu
#SBATCH -o skipper_genome_index.out #redirecting stdout
#SBATCH -p nano #queue 
#SBATCH -N 1 #amount of nodes 
#SBATCH -n 40 #amount of cores 
#SBATCH -t 0:30:00

echo "=========================================================="
echo "Running on node : $HOSTNAME"
echo "Current directory : $PWD"
echo "Current job ID : $SLURM_JOB_ID"
echo "=========================================================="


cd /scratch/martinlab/jasmine/skipper_diff_exp/star_runs

CORES=40
GENOME_DIR=/scratch/martinlab/jasmine/skipper_diff_exp/star_runs/skipper_genome_index
GENOME=/scratch/martinlab/jasmine/skipper_diff_exp/ncbi_skipper_genome/GCF_041222505.1/GCF_041222505.1_WU_Ecla_fem_2.2_genomic.fna
ANNOTATION=/scratch/martinlab/jasmine/skipper_diff_exp/GCF_041222505.1_WU_Ecla_fem_2.2_genomic_cortexedited.gtf

mkdir $GENOME_DIR
STAR --runMode genomeGenerate --runThreadN $CORES --genomeDir $GENOME_DIR --genomeFastaFiles $GENOME --sjdbGTFfile $ANNOTATION --sjdbOverhang 149 --genomeSAindexNbases 14
```


Before mapping I renamed files using this
```
while read -r old new; do   for f in *"$old"*; do     [ -e "$f" ] || continue;      newname="${f//$old/$new}";      mv -- "$f" "$newname";     echo "Renamed: $f -> $newname";   done; done < file_names.txt
```
**Run the first pass of STAR mapping to the <em>E. clarus</em> genome** 
```
#!/bin/sh
#SBATCH -J skipper_rnaseq_star_pass_1_se
#SBATCH --mail-type=END
#SBATCH --mail-user=j.alqassar@gwu.edu
#SBATCH -o star_pass_1_%A_%a.out #redirecting stdout
#SBATCH -p nano #queue 
#SBATCH -n 40 #amount of cores 
#SBATCH -t 0:30:00
#SBATCH --array=0-31%10

echo "=========================================================="
echo "Running on node : $HOSTNAME"
echo "Current directory : $PWD"
echo "Current job ID : $SLURM_JOB_ID"
echo "Job Started:"
date
echo "=========================================================="

ulimit -n 10000

PREFIX=/scratch/martinlab/jasmine/skipper_diff_exp/star_runs
names=($(cat ${PREFIX}/samples))
echo ${names[${SLURM_ARRAY_TASK_ID}]} 

CORES=40
GENOME_DIR=/scratch/martinlab/jasmine/skipper_diff_exp/star_runs/skipper_genome_index
GENOME=/scratch/martinlab/jasmine/skipper_diff_exp/ncbi_skipper_genome/GCF_041222505.1/GCF_041222505.1_WU_Ecla_fem_2.2_genomic.fna
ANNOTATION=/scratch/martinlab/jasmine/skipper_diff_exp/GCF_041222505.1_WU_Ecla_fem_2.2_genomic_cortexedited.gtf
RNAseq_FILES_PATH=/scratch/martinlab/jasmine/skipper_diff_exp/skipper_data_renamed/single_end
OUT_DIR=/scratch/martinlab/jasmine/skipper_diff_exp/star_runs/star_pass1

cd ${OUT_DIR}

STAR --runMode alignReads --runThreadN $CORES --genomeDir $GENOME_DIR --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${OUT_DIR}/${names[${SLURM_ARRAY_TASK_ID}]} --readFilesCommand zcat \
        --sjdbGTFfile $ANNOTATION --limitBAMsortRAM 60000000000 \
        --sjdbOverhang 79 --outSAMstrandField intronMotif --alignSoftClipAtReferenceEnds No --sjdbGTFtagExonParentTranscript gene_id \
        --readFilesIn ${RNAseq_FILES_PATH}/${names[${SLURM_ARRAY_TASK_ID}]}_se_R1.fastq.gz
      
echo "=========================================================="
echo "Job Finished  $SLURM_JOB_ID:"
date
echo "=========================================================="
```
After job is finished rename each log file by sample name
```
for i in star_pass_1*.out; do
	line8=$(sed -n '8p' ${i})
	mv ${i} "star_pass_1_${line8}.out";
	done
```
**Run the second pass of STAR mapping to the <em>E. clarus</em> genome**

```
#!/bin/sh
#SBATCH -J skipper_rnaseq_star_pass_2_pe
#SBATCH --mail-type=END
#SBATCH --mail-user=j.alqassar@gwu.edu
#SBATCH -o star_pass_2_%A_%a.out #redirecting stdout
#SBATCH -p nano #queue 
#SBATCH -n 40 #amount of cores 
#SBATCH -t 0:30:00
#SBATCH --array=0-17%5

echo "=========================================================="
echo "Running on node : $HOSTNAME"
echo "Current directory : $PWD"
echo "Current job ID : $SLURM_JOB_ID"
echo "Job Started:"
date
echo "=========================================================="

ulimit -n 10000

PREFIX=/scratch/martinlab/jasmine/skipper_diff_exp/star_runs
names=($(cat ${PREFIX}/pe_samples))
echo ${names[${SLURM_ARRAY_TASK_ID}]} 

CORES=40
GENOME_DIR=/scratch/martinlab/jasmine/skipper_diff_exp/star_runs/skipper_genome_index_for_paired_end
GENOME=/scratch/martinlab/jasmine/skipper_diff_exp/ncbi_skipper_genome/GCF_041222505.1/GCF_041222505.1_WU_Ecla_fem_2.2_genomic.fna
ANNOTATION=/scratch/martinlab/jasmine/skipper_diff_exp/GCF_041222505.1_WU_Ecla_fem_2.2_genomic_cortexedited.gtf
RNAseq_FILES_PATH=/scratch/martinlab/jasmine/skipper_diff_exp/skipper_data_renamed/paired_end
OUT_DIR=/scratch/martinlab/jasmine/skipper_diff_exp/star_runs/star_pass2
PASS1_DIR=/scratch/martinlab/jasmine/skipper_diff_exp/star_runs/star_pass1

cd ${OUT_DIR}

STAR --runMode alignReads --runThreadN $CORES --genomeDir $GENOME_DIR --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${OUT_DIR}/${names[${SLURM_ARRAY_TASK_ID}]}_pass2_mapped --readFilesCommand zcat \
        --sjdbGTFfile $ANNOTATION --limitBAMsortRAM 60000000000 \
        --sjdbOverhang 149 --outSAMstrandField intronMotif --alignSoftClipAtReferenceEnds No --sjdbGTFtagExonParentTranscript gene_id \
        --outReadsUnmapped Fastx \
        --quantMode GeneCounts \
        --sjdbFileChrStartEnd ${PASS1_DIR}/${names[${SLURM_ARRAY_TASK_ID}]}_STARgenome/sjdbList.out.tab \
        --readFilesIn ${RNAseq_FILES_PATH}/${names[${SLURM_ARRAY_TASK_ID}]}_pe_R1.fastq.gz ${RNAseq_FILES_PATH}/${names[${SLURM_ARRAY_TASK_ID}]}_pe_R2.fastq.gz      

echo "=========================================================="
echo "Job Finished  $SLURM_JOB_ID:"
date
echo "=========================================================="
```
Rename the log files 
```
for i in star_pass_2*.out; do
	line8=$(sed -n '8p' ${i})
	mv ${i} "star_pass_2_${line8}.out";
	done
```
## Read Counts with FeatureCounts
Install subread software which contains FeatureCounts
```
conda create -n subread_2.0.8 -c conda-forge mamba
conda activate subread_2.0.8
conda install -c conda-forge -c bioconda subread
```
Activate the subread enviroment and seperate single and paired-end files so that you can run seperate counting jobs for them 
```
mamba activate subread_2.0.8_JDA
mkdir se_bams 

while read -r prefix; do
  mv "${prefix}"*.bam se_bams/
done < /scratch/martinlab/jasmine/skipper_diff_exp/star_runs/se_samples

mkdir pe_bams 

while read -r prefix; do
  mv "${prefix}"*.bam pe_bams/
done < /scratch/martinlab/jasmine/skipper_diff_exp/star_runs/pe_samples
```
Now run the counts job for the paired-end files using your GTF file with manual annotations 
```
#!/bin/sh
#SBATCH -J skipper_rnaseq_counts_pe
#SBATCH --mail-type=END
#SBATCH --mail-user=j.alqassar@gwu.edu
#SBATCH -o counts_pe.out #redirecting stdout
#SBATCH -p nano #queue 
#SBATCH -n 40 #amount of cores 
#SBATCH -t 0:30:00

echo "=========================================================="
echo "Running on node : $HOSTNAME"
echo "Current directory : $PWD"
echo "Current job ID : $SLURM_JOB_ID"
echo "Job Started:"
date
echo "=========================================================="

ANNOTATION=/scratch/martinlab/jasmine/skipper_diff_exp/GCF_041222505.1_WU_Ecla_fem_2.2_genomic_cortexedited.gtf
OUT=skipper_rnaseq_pe.featurecounts.txt
BAM_FILES=/scratch/martinlab/jasmine/skipper_diff_exp/star_runs/star_pass2/pe_bams
CORES=40

featureCounts -a $ANNOTATION -o $OUT  -T $CORES -p --primary -t exon -g gene_id ${BAM_FILES}/*.bam

echo "=========================================================="
echo "Job Finished  $SLURM_JOB_ID:"
date
echo "=========================================================="
```
Now run the counts job for the single-end files using your GTF file with manual annotations
```
#!/bin/sh
#SBATCH -J skipper_rnaseq_counts_se
#SBATCH --mail-type=END
#SBATCH --mail-user=j.alqassar@gwu.edu
#SBATCH -o counts_se.out #redirecting stdout
#SBATCH -p nano #queue 
#SBATCH -n 40 #amount of cores 
#SBATCH -t 0:30:00

echo "=========================================================="
echo "Running on node : $HOSTNAME"
echo "Current directory : $PWD"
echo "Current job ID : $SLURM_JOB_ID"
echo "Job Started:"
date
echo "=========================================================="

ANNOTATION=/scratch/martinlab/jasmine/skipper_diff_exp/GCF_041222505.1_WU_Ecla_fem_2.2_genomic_cortexedited.gtf
OUT=skipper_rnaseq_se.featurecounts.txt
BAM_FILES=/scratch/martinlab/jasmine/skipper_diff_exp/star_runs/star_pass2/se_bams
CORES=40

featureCounts -a $ANNOTATION -o $OUT  -T $CORES --primary -t exon -g gene_id ${BAM_FILES}/*.bam

echo "=========================================================="
echo "Job Finished  $SLURM_JOB_ID:"
date
echo "=========================================================="
```
## Differential Expression Analysis with DESeq2 in RStudio
An R script with all of the code detailed below is [here](https://github.com/jasalq/)

Set working directory
```
setwd("/Users/jasminealqassar/Library/CloudStorage/GoogleDrive-j.alqassar@gwmail.gwu.edu/.shortcut-targets-by-id/1Hi7WIp_ha7vnyQUclvt-AJRl6AF8y1hj/Martin Lab/@LabData/2025_Skipper_WntA/SSS_Differential_Expression/R_working_dir")
```
Install and load all necessary packages
```
install.packages("tidyverse"); install.packages("stringr"); install.packages("dplyr"); install.packages("DESeq2")
library(tidyverse); library(stringr); library(dplyr); library(DESeq2); library(ggplot2); library(readxl)
```
### Differential expression analysis of wing type (Forewing vs Hindwing) at 36hrs (12% pupal development)
```
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
```


### Differential expression analysis of wing compartments at 36hrs (12% pupal development)
```
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
```
## Heatmap and Count Plots Data Visualization 
### First I need to re-run DEseq to get the counts for all samples in the heatmap except the samples we removed from the entire analysis because of poor quality
```
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
```
### Generating Count Plots

Now make an individual count plot for a specific gene 
```
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
```
Generate one large matrix with many count plots
```
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
```
### Generating A Heatmap with the differentially expressed genes found in the compartment analysis 
```
#load necessary packages
library(readxl);library(DESeq2);library(reshape2);library(ggdendro);library(gridExtra);library(cowplot);library(pals)

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
```
## GO Enrichment Analysis with GO Subsets
Inspired by tutorial [here](https://www.nature.com/articles/s41596-024-01020-z#Sec43)
### Analysis of DEGs using the GO subset slimGO_agr 
Download the slimGO_agr dataset [here](https://geneontology.org/docs/go-subset-guide/)	

**First with Owltools map the slim terms to the existing NCBI GO annotation for the Ecla genome (located On the NCBI FTP Server)**
Install Owltools
```
download release file from here https://github.com/owlcollab/owltools
chmod +x owltools
module load jdk
export PATH=$PATH:/CCAS/groups/martinlab/jasmine/software
owltools -h
```
Now run the mapping of the slim GO terms to the NCBI GO annotation 
```
#!/bin/sh
#SBATCH -J SSS_GO_slim_mapping
#SBATCH --mail-type=ALL
#SBATCH --mail-user=j.alqassar@gwu.edu
#SBATCH -o SSS_GO_slim_mapping.out #redirecting stdout
#SBATCH -p nano #queue 
#SBATCH -N 1 #amount of nodes 
#SBATCH -n 40 #amount of cores 
#SBATCH -t 00:30:00


echo "=========================================================="
echo "Running on node : $HOSTNAME"
echo "Current directory : $PWD"
echo "Current job ID : $SLURM_JOB_ID"
echo "=========================================================="

PATH=$PATH:/CCAS/groups/martinlab/jasmine/software 
module load jdk

cd /scratch/martinlab/jasmine/skipper_diff_exp/GO_term_analyses

owltools go.obo --gaf GCF_041222505.1-RS_2025_04_gene_ontology.gaf --map2slim --subset goslim_agr --write-gaf GCF_041222505.1-RS_2025_04_slim_GO.mapped.gaf
```
**Now back in RStudio perform the GO enrichment analysis with <em>clusterProfiler</em>**

**GO enrichment analysis with the whole HW vs FW results at 36hr**
```
sigGenes_whole_HW_vs_FW<- allGenesHW_vs_FW_36h[!is.na(allGenesHW_vs_FW_36h$padj) & allGenesHW_vs_FW_36h$padj < 0.05, ]

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

# Now plot as a bubble plot
ggplot(df[1:6, ], aes(x = FoldEnrichment, y = Description)) +
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
```
**Perform the analysis on compartment DEGs**
```
library(clusterProfiler)
library(BaseSet)
library(ontologyIndex)

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

custom_order <- c("GO:0003700", "GO:0038023", "GO:0008092", "GO:0030154", "GO:0032502")

df <- df %>%
  # Make Description a factor with levels ordered by custom_order of IDs
  mutate(Description = factor(Description, levels = Description[match(custom_order, ID)])) %>%
  arrange(factor(ID, levels = custom_order))
         
# Now plot a bubble plot
ggplot(df[1:5, ], aes(x = FoldEnrichment, y = Description)) +
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
```
**Now make heatmaps with only certain GO terms from the DEGs from compartment analysis**
First for transcription factors 
```
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
```
Now for signaling related genes 
```
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
```
### Analysis of FW vs HW DEGs using the GO subset slimGO_Drosophila
Download the slimGO_agr dataset [here](https://geneontology.org/docs/go-subset-guide/)	

**First with Owltools map the slim terms to the existing NCBI GO annotation for the <em>E. clarus</em> genome (located On the NCBI FTP Server)**
```
#!/bin/sh
#SBATCH -J SSS_GO_fly_slim_mapping
#SBATCH --mail-type=ALL
#SBATCH --mail-user=j.alqassar@gwu.edu
#SBATCH -o SSS_GO_fly_slim_mapping.out #redirecting stdout
#SBATCH -p nano #queue 
#SBATCH -N 1 #amount of nodes 
#SBATCH -n 40 #amount of cores 
#SBATCH -t 00:30:00


echo "=========================================================="
echo "Running on node : $HOSTNAME"
echo "Current directory : $PWD"
echo "Current job ID : $SLURM_JOB_ID"
echo "=========================================================="

module load jdk
PATH=$PATH:/CCAS/groups/martinlab/jasmine/software/

cd /scratch/martinlab/jasmine/skipper_diff_exp/GO_term_analyses

owltools go.obo --gaf GCF_041222505.1-RS_2025_04_gene_ontology.gaf --map2slim --subset goslim_drosophila --write-gaf GCF_041222505.1-RS_2025_04_drosophila_slim_GO.mapped.gaf
```
**Now back in RStudio perform the GO enrichment analysis with <em>clusterProfiler</em>**
```
library(clusterProfiler)
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
```
