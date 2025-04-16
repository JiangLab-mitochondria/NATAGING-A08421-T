## mgatk
set -e
echo ""
echo "kidney-3W-TE-G14102A-4475-80 mgatk start"
echo ""
outdir="/mnt/transposon1/zhangyanxiaoLab/chaiguoshi/projects/jiangmin-lab-project/work/zhang-le-ping-project/resutls/mtscATAC/mgatk-results/"
bamfile="/mnt/transposon1/zhangyanxiaoLab/chaiguoshi/projects/jiangmin-lab-project/work/zhang-le-ping-project/resutls/mtscATAC/cellranger-atac-results/kidney-3W-TE-G14102A-4475-80-masked/outs/possorted_bam_rmdup_DNA.bam"
barcode="/mnt/transposon1/zhangyanxiaoLab/chaiguoshi/projects/jiangmin-lab-project/work/zhang-le-ping-project/resutls/mtscATAC/archr-results/clustering/"
mgatk tenx -g mm10 \
           -bt CB \
           --ncores 6 \
           -o ${outdir}kidney-3W-TE-G14102A-4475-80 \
           -n kidney-3W-TE-G14102A-4475-80 \
           -i ${bamfile} \
           -b ${barcode}kidney-3W-TE-G14102A-4475-80-cell-barcode.txt
echo "kidney-3W-TE-G14102A-4475-80 mgatk end"

