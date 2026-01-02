# Pipeline for painting Merian elements on a genome
Notes by Jasmine Alqassar 2025 using the scripts built by Charlotte Wright: https://github.com/charlottewright/lep_busco_painter.git
## Installation
### Use the shared Martin Lab conda enviroment
```
mamba activate buscopaint_JDA
```
### Or create your own enviroment 
```
mamba create -n buscopaint python=3.9 
mamba activate buscopaint 
mamba install samtools 
mamba install -c conda-forge r-base
mamba install -c r r-tidyverse
mamba install -c bioconda r-optparse
```
### Move the necessary scripts to your working directory 
```
wget https://github.com/charlottewright/lep_busco_painter/blob/904ec735c9e78b3747dfe5dbed2c680c593f3d2a/buscopainter.py
wget https://github.com/charlottewright/lep_busco_painter/blob/904ec735c9e78b3747dfe5dbed2c680c593f3d2a/reference_data/Merian_elements_full_table.tsv
cp /CCAS/groups/martinlab/jasmine/scripts/plot_buscopainter_editedJA.R . 
```
**Note: You can also download the R script which I have edited from <a href="https://github.com/DNAcrobatics/Martin-Lab-Bioinformatics/blob/8d869d052bd958bb1bc5fb072fb6a859e782d887/Genome_Assembly/plot_buscopainter_editedJA.R">here</a> if it is easier**
## Run BUSCO via Compleasm 
**Use the shared Martin Lab conda enviroment**
```
mamba activate compleasm_0.2.6_JDA
```
**Run Compleasm**
```
#!/bin/sh
#SBATCH -J sss_compleasm
#SBATCH --mail-type=ALL
#SBATCH --mail-user=j.alqassar@gwu.edu
#SBATCH -o sss.out #redirecting stdout
#SBATCH -p defq #queue 
#SBATCH -N 1 #amount of nodes 
#SBATCH -n 32 #amount of cores 
#SBATCH -t 2:00:00

echo "=========================================================="
echo "Running on node : $HOSTNAME"
echo "Current directory : $PWD"
echo "Current job ID : $SLURM_JOB_ID"
echo "Job started at: $(date)"
echo "=========================================================="

#make sure to activate your compleasm conda enviroment prior to job submission: conda activate compleasm_0.2.6 

#compleasm run -a [assembly] -o [output_dir] -t [threads] \
        #-l lepidoptera -L [lineage_path]
cd /scratch/martinlab/jasmine/sss_genome_compleasm
compleasm run -a /scratch/martinlab/jasmine/sss_genome_compleasm/e_clarus_purged.fa \
        -o sss_output_dir -t 32 \
        -l lepidoptera -L /scratch/martinlab/jasmine/sss_genome_compleasm/mb_downloads


echo "=========================================================="
echo "Job finished at: $(date)"
echo "=========================================================="
```
## Use the BUSCO output (full_table.tsv) to map Merian elements 
```
mamba activate buscopaint_JDA
./buscopainter.py -r Merian_elements_full_table.tsv -q full_table.tsv
./plot_buscopainter_editedJA.R -f buscopainter_complete_location.tsv -p "Ecla" -i e_clarus_purged.fa.fai -m TRUE
```
* If you want only the rearranged Merian elements painted add **-d TRUE**  
* If you don't have a fasta index file before this use Samtools to make one: 
```
samtools faidx [fasta file]
```

