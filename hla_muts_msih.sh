#!/bin/bash
cd /home/data/t070331/WES/LINXIN/hla_muts
cat /home/data/t070331/WES/LINXIN/Results/MSIH/Coverage/tasks.tsv | while read task_file
do
    arr=(${task_file})
    normal_file=${arr[1]}
    tumor_file=${arr[0]}
    # if [ ! -e /home/data/t070331/WES/LINXIN/LOHHLA_Result/lohhla_output/${tumor_file}/${tumor_file}.*.DNA.HLAlossPrediction_CI.*.tsv ]
    if [ ! -e /home/data/t070331/WES/LINXIN/hla_muts/${tumor_file}/${tumor_file}.mutect.unfiltered.annotated ]
    then
        if [ ! -f /home/data/t070331/WES/LINXIN/hla_muts/$tumor_file.bam ]
        then
            genounzip /home/data/t070331/WES/LINXIN/Results/Pending_BAMS/${tumor_file}.bam.genozip --output /home/data/t070331/WES/LINXIN/hla_muts/${tumor_file}.bam -@ 35 --index
        fi
        # if [ -f /home/data/t070331/WES/LINXIN/hla_muts/$tumor_file.bam ] && [ ! -f /home/data/t070331/WES/LINXIN/hla_muts/$tumor_file.bam.bai ]
        # then
        #     samtools index -@ 35 /home/data/t070331/WES/LINXIN/hla_muts/${tumor_file}.bam
        # fi
        if [ ! -f /home/data/t070331/WES/LINXIN/hla_muts/$normal_file.bam ]
        then
            genounzip /home/data/t070331/WES/LINXIN/Results/Pending_BAMS/${normal_file}.bam.genozip --output /home/data/t070331/WES/LINXIN/hla_muts/${normal_file}.bam -@ 35 --index
            # samtools index -@ 35 ${normal_file}.bam
        fi
        # if [ -f /home/data/t070331/WES/LINXIN/hla_muts/$normal_file.bam ] && [ ! -f /home/data/t070331/WES/LINXIN/hla_muts/$normal_file.bam.bai ]
        # then
        #     samtools index -@ 35 /home/data/t070331/WES/LINXIN/hla_muts/${normal_file}.bam
        # fi
        if [ -f /home/data/t070331/WES/LINXIN/hla_muts/${normal_file}.bam ] && [ -f /home/data/t070331/WES/LINXIN/hla_muts/${tumor_file}.bam ] && [ -f /home/data/t070331/WES/LINXIN/hla_muts/${normal_file}.bam.bai ] && [ -f /home/data/t070331/WES/LINXIN/hla_muts/${tumor_file}.bam.bai ]
        then
        mkdir /home/data/t070331/WES/LINXIN/hla_muts/${tumor_file}
        shell_call_hla_mutations_from_type /home/data/t070331/WES/LINXIN/hla_muts/${normal_file}.bam /home/data/t070331/WES/LINXIN/hla_muts/${tumor_file}.bam /home/data/t070331/WES/LINXIN/Results/MSIH/HLA_Prediction/${normal_file}/winners.hla.txt hg19 STDFQ /home/data/t070331/WES/LINXIN/hla_muts/${tumor_file} \
        1>>/home/data/t070331/WES/LINXIN/hla_muts/${tumor_file}/${tumor_file}_hlamuts.log 2>&1
        shell_annotate_hla_mutations ${tumor_file} /home/data/t070331/WES/LINXIN/hla_muts/${tumor_file} \
        1>>/home/data/t070331/WES/LINXIN/hla_muts/${tumor_file}/${tumor_file}_hlamuts_annotated.log 2>&1
        fi
        if [ -e /home/data/t070331/WES/LINXIN/hla_muts/${tumor_file}/${tumor_file}.mutect.unfiltered.annotated ]
        then
            rm /home/data/t070331/WES/LINXIN/hla_muts/${normal_file}.bam
            rm /home/data/t070331/WES/LINXIN/hla_muts/${tumor_file}.bam
        fi
    fi
done
