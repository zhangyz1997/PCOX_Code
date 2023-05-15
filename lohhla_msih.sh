#!/bin/bash
## 计算杀千刀的LoH, 在lohhla环境下跑

cd /home/data/t070331/WES/LINXIN/LOHHLA_Result
cat /home/data/t070331/WES/LINXIN/Results/MSIH/Coverage/tasks.tsv | while read task_file
do
    arr=(${task_file})
    normal_file=${arr[1]}
    tumor_file=${arr[0]}
    if [ ! -e /home/data/t070331/WES/LINXIN/LOHHLA_Result/lohhla_output/${tumor_file}/${tumor_file}.*.DNA.HLAlossPrediction_CI.*.tsv ]
    then
        if [ ! -f /home/data/t070331/WES/LINXIN/LOHHLA_Result/$tumor_file.bam ]
        then
            genounzip /home/data/t070331/WES/LINXIN/Results/Pending_BAMS/${tumor_file}.bam.genozip --output /home/data/t070331/WES/LINXIN/LOHHLA_Result/${tumor_file}.bam -@ 35 
            samtools index -@ 35 ${tumor_file}.bam
        fi
        if [ ! -f /home/data/t070331/WES/LINXIN/LOHHLA_Result/$normal_file.bam ]
        then
            genounzip /home/data/t070331/WES/LINXIN/Results/Pending_BAMS/${normal_file}.bam.genozip --output /home/data/t070331/WES/LINXIN/LOHHLA_Result/${normal_file}.bam -@ 35 
            samtools index -@ 35 ${normal_file}.bam
        fi
        if [ -f /home/data/t070331/WES/LINXIN/LOHHLA_Result/${normal_file}.bam ] && [ -f /home/data/t070331/WES/LINXIN/LOHHLA_Result/${tumor_file}.bam ] && [ -f /home/data/t070331/WES/LINXIN/LOHHLA_Result/${normal_file}.bam.bai ] && [ -f /home/data/t070331/WES/LINXIN/LOHHLA_Result/${tumor_file}.bam.bai ]
        then
        # mv /home/data/t070331/WES/LINXIN/LOHHLA_Result/lohhla_output/flagstat/*.flagstat /home/data/t070331/WES/LINXIN/LOHHLA_Result/lohhla_output/flagstat/archived_flagstat/
        mkdir /home/data/t070331/WES/LINXIN/LOHHLA_Result/lohhla_output/${tumor_file}
        Rscript --verbose /home/data/t070331/Software/LOHHLA/LOHHLAscript.R \
            --patientId $tumor_file \
            --outputDir /home/data/t070331/WES/LINXIN/LOHHLA_Result/lohhla_output/${tumor_file}/ \
            --LOHHLA_loc /home/data/t070331/Software/LOHHLA \
            --normalBAMfile /home/data/t070331/WES/LINXIN/LOHHLA_Result/${normal_file}.bam \
            --tumorBAMfile /home/data/t070331/WES/LINXIN/LOHHLA_Result/${tumor_file}.bam \
            --hlaPath /home/data/t070331/WES/LINXIN/NeoAntigen/HLA_Subtypes/Tumor_Based_HLA_MSIH/${tumor_file}_mutect2_somatic_filtered/winners.hla.txt \
            --HLAfastaLoc /home/data/t070331/Software/hla-polysolver/data/abc_complete.fasta \
            --CopyNumLoc /home/data/t070331/WES/LINXIN/LOHHLA_Result/solutions_msih.txt \
            --mappingStep TRUE --minCoverageFilter 10 --fishingStep TRUE --cleanUp TRUE \
            --gatkDir /home/data/t070331/Software/picard-2.27.5-0/ \
            --novoDir /home/data/t070331/Software/novocraft/ --BAMDir /home/data/t070331/WES/LINXIN/LOHHLA_Result/ --HLAexonLoc /home/data/t070331/Software/LOHHLA/data/hla.dat 
        fi
        wait
        if [ -e /home/data/t070331/WES/LINXIN/LOHHLA_Result/lohhla_output/${tumor_file}/${tumor_file}.*.DNA.HLAlossPrediction_CI.*.tsv ]
        then
            rm /home/data/t070331/WES/LINXIN/LOHHLA_Result/${normal_file}.bam
            rm /home/data/t070331/WES/LINXIN/LOHHLA_Result/${tumor_file}.bam
        fi
    fi
done
