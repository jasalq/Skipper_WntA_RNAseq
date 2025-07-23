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
