#!/bin/bash

cd /home/data/t070331/WES/LINXIN/Results/MSIH/Coverage/
# mkdir /home/data/t070331/WES/LINXIN/PureCN/OUTPUT_DIR/denovo_seg
cat /home/data/t070331/WES/LINXIN/Results/MSIH/Coverage/tasks.tsv | while read task_file
do
    arr=(${task_file})
    normal_file=${arr[1]}
    tumor_file=${arr[0]}
    sex=${arr[2]}

    if [ ! -f /home/data/t070331/WES/LINXIN/PureCN/OUTPUT_DIR/msih/${tumor_file}.rds ] && [ -f /home/data/t070331/WES/LINXIN/Results/MSIH/Coverage/${tumor_file}_coverage_loess_qc.txt ]
    then
        singularity exec /home/data/t070331/WES/LINXIN/PureCN/purecn_latest.sif Rscript /usr/local/lib/R/site-library/PureCN/extdata/PureCN.R --out /home/data/t070331/WES/LINXIN/PureCN/OUTPUT_DIR/msih \
        --tumor /home/data/t070331/WES/LINXIN/Results/MSIH/Coverage/${tumor_file}_coverage_loess.txt.gz \
        --sampleid ${tumor_file} \
        --vcf /home/data/t070331/WES/LINXIN/Results/MSIH/PureCN_VCFS/${tumor_file}_mutect2.vcf \
        --fun-segmentation PSCBS \
        --normaldb /home/data/t070331/WES/LINXIN/PureCN/reference/normalDB_agilent_v7_hg19.rds \
        --mapping-bias-file /home/data/t070331/WES/LINXIN/PureCN/reference/mapping_bias_agilent_v7_hg19.rds \
        --intervals /home/data/t070331/WES/LINXIN/PureCN/reference/baits_hg19_intervals.txt \
        --snp-blacklist /home/data/t070331/WES/LINXIN/PureCN/reference/hg19.simpleRepeat.bed \
        --genome hg19 \
        --model betabin \
        --sex ${sex} \
        --force --post-optimize --seed 123 --cores 35
    fi
done