#!/bin/bash

##-------------------------
# @Author: Chai guoshi
# @Date: 2024/9/22
##-------------------------


#slurm options
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH -J zhang-le-ping-mtscATAC-3W-TE-G14102A-4475-80
#SBATCH -c 16
#SBATCH -o %j_cellranger-atac_20240922.log

## user's own commands below
echo -e "=================\nStart.\n================="

module load cellranger-atac/2.1.0


# 3W-TE-G14102A-4475-80
# mut4475

mm10_masked="/storage/jiangminLab/chaiguoshi/references/cellranger-atac/refdata-cellranger-arc-mm10-mitochondria-masked"
#mm10="/storage/publicdata/ref/cellranger-atac/refdata-cellranger-arc-mm10-2020-A-2.0.0"

fq_file="/storage/jiangminLab/chaiguoshi/projects/zhang-le-ping-projects/raw-data/mtscATAC/X101SC24076191-Z01-J004/Rawdata/mut4475/"
out_dir="/storage/jiangminLab/chaiguoshi/projects/zhang-le-ping-projects/results/mtscATAC/cellranger-atac-results/"

cd ${out_dir}


echo "############################"
echo "cellranger-atac count start"
echo "3W-TE-G14102A-4475-80 mut4475 sample start"
cellranger-atac count --id=kidney-3W-TE-G14102A-4475-80-masked \
                      --reference=${mm10_masked} \
                      --fastqs=${fq_file} \
                      --sample=mut4475-1 \
                      --localcores=16
                        

echo "3W-TE-G14102A-4475-80 mut4475 sample end"  
echo "cellranger-atac count end"
echo "############################"
echo ""




