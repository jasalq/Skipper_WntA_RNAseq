# Skipper_WntA
## Pipeline Outline 
1. [RNAseq Quality Assessment](#assessment-of-rna-sequencing-quality-using-fastqc) 
2. [Adapter Trimming](#adapter-trimming-using-fastp)
3. [Read Mapping to the Reference Genome]

### Tools Used 

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
## Adapter Trimming using Fastp
Trim the adapters and polyG tails with Fastp

## Read Mapping

**Download <em>E. clarus</em> RefSeq genome and annotation from NCBI using NCBI Datasets tool** 
```
mamba activate ncbi_datsets
datasets download genome accession GCF_041222505.1 --include gtf,rna,cds,protein,genome,seq-report
unzip ncbi_dataset
mv /scratch/martinlab/jasmine/skipper_diff_exp/ncbi_skipper_genome/ncbi_dataset/data/GCF_041222505.1 ../../
```
**In addition download the list of all annotated genes by clicking "view annotated genes" on the genome release page selecting all and hit download table while clicking "one sequence per gene" option**
I manually edited the GeneID column to have a LOC prefix

## Homology Searches Using Reciprocal BLASTp to FlyBase
**use the protein.faa file from the NCBI download to blast to flybase**
**find this file -> dmel-all-translation.fasta from the latest FlyBase release and unzip**

``
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

# Set working directory
setwd("/Users/jasminealqassar/Library/CloudStorage/GoogleDrive-j.alqassar@gwmail.gwu.edu/.shortcut-targets-by-id/1Hi7WIp_ha7vnyQUclvt-AJRl6AF8y1hj/Martin Lab/@LabData/2025_Skipper_WntA/SSS_Differential_Expression/R_working_dir")

# Install necessary packages
install.packages("tidyverse"); install.packages("stringr"); install.packages("dplyr"); install.packages("DESeq2")

# Load necessary packages
library(tidyverse); library(stringr); library(dplyr); library(DESeq2); library(ggplot2); library(readxl)

# Read in the results
flybase_results <- read_tsv("all_Ecla_proteins_for_flybase.out", col_names=FALSE)

#Download and read in the following tables from the latest [FlyBase genome release GUI page ](#https://flybase.org/downloads/bulkdata) and clean up column names before reading in
  flybase_prot_to_Symbol <- read_tsv("dmel_unique_protein_isoforms_fb_2025_02.tsv", col_names=TRUE, comment="#")
  flybase_gn_to_bpp <- read_tsv("fbgn_fbtr_fbpp_expanded_fb_2025_02.tsv", col_names=TRUE, comment="#") 
  flybase_gn_summary <- read_tsv("automated_gene_summaries_fb_2025_02.tsv", col_names=TRUE, comment="#") 
  
  
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

##STAR MAPPING
mamba activate star_2.7.11b_JDA 

```
## STAR Mapping of RNAseq data to the <em>Plodia</em> reference genome (GCF_027563975.2)

**Download <em>Plodia</em> RefSeq genome from NCBI using NCBI Datasets tool** 
```
mamba activate NCBI_datsets
datasets download genome accession GCF_027563975.2 
unzip ncbi_dataset
mv ncbi_dataset/data/GCF_027563975.2/GCF_027563975.2_ilPloInte3.2_genomic.fna ../../../
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

# when you do star mapping try counting by gene_id field (really symbol in the annotation table) in the GTF so that it collapses different XM splice isoforms by loci using the option --sjdbGTFtagExonParentTranscript gene_id
