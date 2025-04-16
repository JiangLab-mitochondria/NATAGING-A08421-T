#!/bin/bash

##-------------------------
# @Author: Chai guoshi
# @Date: 2024/09/07
##-------------------------


#slurm options
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH -J zhang-le-ping-scRNA-seq
#SBATCH -c 16
#SBATCH -o %j_cellranger-WT742_20240907.log

## user's own commands below
echo -e "=================\nStart.\n================="

module load cellranger/7.1.0

mm10="/storage/publicdata/ref/cellranger/refdata-gex-mm10-2020-A"
fq_file="/storage/jiangminLab/chaiguoshi/projects/zhang-le-ping-projects/raw-data/scRNA/"
out_dir="/storage/jiangminLab/chaiguoshi/projects/zhang-le-ping-projects/results/scRNA/cellranger-results/"

cd ${out_dir}

##WT742
## 75W-TE-G14102A-742-0

echo "WT742 75W-TE-G14102A-742-0 sample start"
cellranger count --id kidney-75W-TE-G14102A-742-0 \
                 --description kidney-75W-TE-G14102A-742-0 \
                 --transcriptome ${mm10} \
                 --fastqs ${fq_file}X101SC24086742-Z01-J001/Rawdata/WT742 \
                 --sample WT742-1,WT742-2,WT742-3,WT742-4 \
                 --localcores 16
echo "WT742 75W-TE-G14102A-742-0 sample end"  





