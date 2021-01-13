#!/bin/sh
set -e

conda install -y fastqc
conda install -y fastp
conda install -y trimmomatic
conda install -y cutadapt
conda install -y megahit
conda install -y spades
conda install -y velvet
conda install -y quast
conda install -y blast
conda install -y blat
conda install -y prokka
conda install -y prodigal
conda install -y snap-aligner
conda install -y bowtie2
conda install -y samtools
conda install -y minimap2
conda install -y biopython bcbio-gff xlrd pyyaml
pip install datasketch

exit 0
