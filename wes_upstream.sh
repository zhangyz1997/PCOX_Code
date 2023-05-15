#!/bin/bash
# conda activate wes
# sed 's/^/chr/' /home/data/refdir/server/annotation/CCDS/human/hg19_exon.bed > ~/WES/LINXIN/hg19_exon.bed
wes_dir=/home/data/t070331/WES/LINXIN/Results/MSIH/
INDEX=~/Ref/bwa2_index/hg19.fa
snp=/home/data/t070331/WES/LINXIN/dbsnp_138.hg19.vcf.gz
indel=/home/data/t070331/WES/LINXIN/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
ref=~/WES/LINXIN/ucsc.hg19.fasta
small_vcf=~/WES/LINXIN/small_exac_common_3_hg19.vcf
GATK=/home/data/t070331/Software/gatk-4.3.0.0/gatk
tmp_dir=/home/data/t070331/WES/LINXIN/gatk_tmp/
# export PATH="/usr/lib/jvm/java-1.8.0-openjdk-amd64/bin:$PATH"
cd $wes_dir
cat ./normals.txt | while read file
# conda activate wes
do
    cd $wes_dir$file
    normal_file=`cat normal.txt`
    cd "${wes_dir}${file}/${normal_file}"
    if [ ! -f ${wes_dir}${file}/${normal_file}_bqsr.bam ]
    then
        if [ ! -f ${wes_dir}${file}/${normal_file}.bam ]
        then
            wes_file_1=$(ls | grep -E "R1.fastq.genozip|1.fq.genozip")
            wes_file_2=$(ls | grep -E "R2.fastq.genozip|2.fq.genozip")
            genounzip $wes_file_1 --output 1.fastq.gz -e ~/WES/LINXIN/hg38.ref.genozip
            genounzip $wes_file_2 --output 2.fastq.gz -e ~/WES/LINXIN/hg38.ref.genozip
            # 质控,去接头
            fastp -i 1.fastq.gz -I 2.fastq.gz -o 1.fastped.fastq.gz -O 2.fastped.fastq.gz --html ${normal_file}.html --json ${normal_file}.json
            # BWA-MEM2比对
            rm 1.fastq.gz 2.fastq.gz
            bwa-mem2 mem -M -t 20 -R "@RG\tID:${normal_file}\tSM:${normal_file}\tLB:WXS\tPL:Illumina" ${INDEX} "1.fastped.fastq.gz" "2.fastped.fastq.gz" | samtools sort -@ 14 -m 1G -o  "${wes_dir}${file}/${normal_file}.bam" -
            samtools index ${wes_dir}${file}/${normal_file}.bam ${wes_dir}${file}/${normal_file}.bam.bai -@ 30
            rm 1.fastped.fastq.gz 2.fastped.fastq.gz
            if [ ! -f /home/data/t070331/WES/LINXIN/Results/MSIH/HLA_Prediction/${normal_file}/winners.hla.txt ]
            then
                conda run -n polysolver shell_call_hla_type ${wes_dir}${file}/${normal_file}.bam Asian 1 hg19 STDFQ 0 /home/data/t070331/WES/LINXIN/Results/MSIH/HLA_Prediction/${normal_file}
            fi
            if [ ! -f /home/data/t070331/WES/LINXIN/Results/MSIH/Coverage/${normal_file}_coverage_loess.txt.gz ]
            then
            singularity exec /home/data/t070331/WES/LINXIN/PureCN/purecn_latest.sif Rscript /usr/local/lib/R/site-library/PureCN/extdata/Coverage.R \
            --out-dir /home/data/t070331/WES/LINXIN/Results/MSIH/Coverage \
            --bam ${wes_dir}${file}/${normal_file}.bam \
            --intervals /home/data/t070331/WES/LINXIN/PureCN/reference/baits_hg19_intervals.txt  --cores 8
            fi
            if [ -f ${wes_dir}${file}/${normal_file}.bam ] && [ -f /home/data/t070331/WES/LINXIN/Results/MSIH/HLA_Prediction/${normal_file}/winners.hla.txt ] && [ -f /home/data/t070331/WES/LINXIN/Results/MSIH/Coverage/${normal_file}_coverage_loess.txt.gz ]
            then
                rm 1.fastq.gz 2.fastq.gz 1.fastped.fastq.gz 2.fastped.fastq.gz $wes_file_1 $wes_file_2
            fi
        fi
        mkdir ~/WES/LINXIN/Results/$normal_file
        mv ${wes_dir}${file}/${normal_file}/${normal_file}.html ~/WES/LINXIN/Results/${normal_file}/${normal_file}.html
        mv ${wes_dir}${file}/${normal_file}/${normal_file}.json ~/WES/LINXIN/Results/${normal_file}/${normal_file}.json
        if [ ! -f ${wes_dir}${file}/${normal_file}_marked.status ]
        then
            echo "start MarkDuplicates for ${normal_file}" `date`
            $GATK --java-options "-Xmx20G -Djava.io.tmpdir=${tmp_dir}" MarkDuplicatesSpark \
            -I ${wes_dir}${file}/${normal_file}.bam \
            --remove-all-duplicates \
            -O ${wes_dir}${file}/${normal_file}_marked.bam \
            -M ${wes_dir}${file}/${normal_file}.metrics --spark-master local[8] \
            1>${wes_dir}${file}/${normal_file}_log.mark 2>&1 
            
            if [ $? -eq 0 ]
            then
                touch ${wes_dir}${file}/${normal_file}_marked.status
            fi
            echo "end MarkDuplicates for ${normal_file}" `date`
            samtools index -@ 30 -m 4G -b ${wes_dir}${file}/${normal_file}_marked.bam ${wes_dir}${file}/${normal_file}_marked.bai
        fi
        # BQSR
        $GATK --java-options "-Xmx20G -Djava.io.tmpdir=${tmp_dir} -XX:ParallelGCThreads=2"  BaseRecalibratorSpark \
        -R $ref  \
        -I "${wes_dir}${file}/${normal_file}_marked.bam"  \
        --known-sites ${snp} \
        --known-sites ${indel} \
        -O "${wes_dir}${file}/${normal_file}_recal.table" --spark-master local[8]
        $GATK --java-options "-Xmx20G -Djava.io.tmpdir=${tmp_dir} -XX:ParallelGCThreads=2"  ApplyBQSRSpark \
        -R $ref  \
        -I "${wes_dir}${file}/${normal_file}_marked.bam"  \
        -bqsr "${wes_dir}${file}/${normal_file}_recal.table" \
        -O "${wes_dir}${file}/${normal_file}_bqsr.bam" --spark-master local[8]
    fi
    if [ -f ${wes_dir}${file}/${normal_file}_bqsr.bam ]
    then
        rm ${wes_dir}${file}/${normal_file}_marked.bam
    fi
    if [ ! -f ${wes_dir}${file}/${normal_file}_getpileupsummaries.table ]
    then
        $GATK --java-options "-Xmx20G -Djava.io.tmpdir=${tmp_dir}" GetPileupSummaries \
        -I ${wes_dir}${file}/${normal_file}_bqsr.bam \
        -V ${small_vcf}            \
        -L ~/WES/LINXIN/S31285117_Regions.bed                  \
        -O ${wes_dir}${file}/${normal_file}_getpileupsummaries.table
    fi
    echo "Create Result Folder Complete"
    cd $wes_dir$file
    echo "Reading Tumor File"
    cat ./tumor.txt
    cat ./tumor.txt | while read tumor_file
    do
        echo "Starting analysis of tumor BAM: ${tumor_file}"
        cd "${wes_dir}${file}/${tumor_file}"
        if [ ! -f ${wes_dir}${file}/${tumor_file}_bqsr.bam ]
        then
            if [ ! -f ${wes_dir}${file}/${tumor_file}.bam ]
            then
                wes_file_1=$(ls | grep -E "R1.fastq.genozip|1.fq.genozip")
                wes_file_2=$(ls | grep -E "R2.fastq.genozip|2.fq.genozip")
                genounzip $wes_file_1 --output 1.fastq.gz -e ~/WES/LINXIN/hg38.ref.genozip
                genounzip $wes_file_2 --output 2.fastq.gz -e ~/WES/LINXIN/hg38.ref.genozip
                # 质控,去接头
                fastp -i 1.fastq.gz -I 2.fastq.gz -o 1.fastped.fastq.gz -O 2.fastped.fastq.gz --html ${tumor_file}.html --json ${tumor_file}.json
                # BWA-MEM2比对
                rm 1.fastq.gz 2.fastq.gz
                bwa-mem2 mem -M -t 20 -R "@RG\tID:${tumor_file}\tSM:${tumor_file}\tLB:WXS\tPL:Illumina" ${INDEX} "1.fastped.fastq.gz" "2.fastped.fastq.gz" | samtools sort -@ 14 -m 1G -o  "${wes_dir}${file}/${tumor_file}.bam" -
                samtools index ${wes_dir}${file}/${tumor_file}.bam ${wes_dir}${file}/${tumor_file}.bam.bai -@ 30
                rm 1.fastped.fastq.gz 2.fastped.fastq.gz
            fi
            if [ ! -f /home/data/t070331/WES/LINXIN/Results/MSIH/Coverage/${tumor_file}_coverage_loess.txt.gz ]
            then
                singularity exec /home/data/t070331/WES/LINXIN/PureCN/purecn_latest.sif Rscript /usr/local/lib/R/site-library/PureCN/extdata/Coverage.R \
                --out-dir /home/data/t070331/WES/LINXIN/Results/MSIH/Coverage \
                --bam ${wes_dir}${file}/${tumor_file}.bam \
                --intervals /home/data/t070331/WES/LINXIN/PureCN/reference/baits_hg19_intervals.txt  --cores 8
            fi
            if [ -f ${wes_dir}${file}/${tumor_file}.bam ] && [ -f /home/data/t070331/WES/LINXIN/Results/MSIH/Coverage/${tumor_file}_coverage_loess.txt.gz ]
            then
                rm 1.fastq.gz 2.fastq.gz 1.fastped.fastq.gz 2.fastped.fastq.gz $wes_file_1 $wes_file_2
            fi
            if [ ! -f ${wes_dir}${file}/${tumor_file}_marked.status ]
            then
                echo "start MarkDuplicates for ${tumor_file}" `date`
                $GATK --java-options "-Xmx20G -Djava.io.tmpdir=${tmp_dir}" MarkDuplicatesSpark \
                -I ${wes_dir}${file}/${tumor_file}.bam \
                --remove-all-duplicates \
                -O ${wes_dir}${file}/${tumor_file}_marked.bam \
                -M ${wes_dir}${file}/${tumor_file}.metrics --spark-master local[8] \
                1>${wes_dir}${file}/${tumor_file}_log.mark 2>&1 
                
                if [ $? -eq 0 ]
                then
                    touch ${wes_dir}${file}/${tumor_file}_marked.status
                fi
                echo "end MarkDuplicates for ${tumor_file}" `date`
                samtools index -@ 30 -m 4G -b ${wes_dir}${file}/${tumor_file}_marked.bam ${wes_dir}${file}/${tumor_file}_marked.bai
            fi
            # BQSR
            $GATK --java-options "-Xmx20G -Djava.io.tmpdir=${tmp_dir} -XX:ParallelGCThreads=2"  BaseRecalibratorSpark \
            -R $ref  \
            -I "${wes_dir}${file}/${tumor_file}_marked.bam"  \
            --known-sites ${snp} \
            --known-sites ${indel} \
            -O "${wes_dir}${file}/${tumor_file}_recal.table" --spark-master local[8]
            $GATK --java-options "-Xmx20G -Djava.io.tmpdir=${tmp_dir} -XX:ParallelGCThreads=2"  ApplyBQSRSpark \
            -R $ref  \
            -I "${wes_dir}${file}/${tumor_file}_marked.bam"  \
            -bqsr "${wes_dir}${file}/${tumor_file}_recal.table" \
            -O "${wes_dir}${file}/${tumor_file}_bqsr.bam" --spark-master local[8]
        fi
        if [ ! -f ${wes_dir}${file}/${tumor_file}_getpileupsummaries.table ]
        then
            $GATK --java-options "-Xmx20G -Djava.io.tmpdir=${tmp_dir}" GetPileupSummaries \
            -I ${wes_dir}${file}/${tumor_file}_bqsr.bam \
            -V ${small_vcf}            \
            -L ~/WES/LINXIN/S31285117_Regions.bed                  \
            -O ${wes_dir}${file}/${tumor_file}_getpileupsummaries.table
        fi
        if [ -f ${wes_dir}${file}/${tumor_file}_bqsr.bam ]
        then
            rm ${wes_dir}${file}/${tumor_file}_marked.bam
        fi
        if [ ! -f ${wes_dir}${file}/${tumor_file}_mutect2.vcf ]
        then
            $GATK  --java-options "-Xmx20G -Djava.io.tmpdir=${tmp_dir}"  Mutect2 -R ${ref} \
            -I "${wes_dir}${file}/${tumor_file}_bqsr.bam" \
            -I "${wes_dir}${file}/${normal_file}_bqsr.bam" \
            -L ~/WES/LINXIN/S31285117_Regions.bed  \
            -normal ${normal_file}  \
            -tumor ${tumor_file}  \
            --panel-of-normals /home/data/t070331/WES/LINXIN/NeoAntigen/Pending_BAM/linxin.pon.vcf.gz    \
            --germline-resource ~/WES/LINXIN/af-only-gnomad.raw.sites.hg19.vcf.gz  \
            --genotype-germline-sites true \
            --genotype-pon-sites true  --native-pair-hmm-threads 20  \
            --f1r2-tar-gz ${wes_dir}${file}/${tumor_file}_f1r2.tar.gz \
            -O ${wes_dir}${file}/${tumor_file}_mutect2.vcf
        fi
        if [ ! -f ${wes_dir}${file}/${tumor_file}_contamination.table ]
        then
            $GATK --java-options "-Xmx20G -Djava.io.tmpdir=${tmp_dir}" CalculateContamination \
            -I ${wes_dir}${file}/${tumor_file}_getpileupsummaries.table         \
            -matched ${wes_dir}${file}/${normal_file}_getpileupsummaries.table \
            --tumor-segmentation ${wes_dir}${file}/${tumor_file}_segments.table \
            -O ${wes_dir}${file}/${tumor_file}_contamination.table              \
            1>>${wes_dir}${file}/${tumor_file}_mutect.log 2>&1
        fi
        if [ ! -f ${wes_dir}${file}/${tumor_file}_read-orientation-model.tar.gz ]
        then
            $GATK --java-options "-Xmx20G -Djava.io.tmpdir=${tmp_dir}" LearnReadOrientationModel \
            -I ${wes_dir}${file}/${tumor_file}_f1r2.tar.gz                  \
            -O ${wes_dir}${file}/${tumor_file}_read-orientation-model.tar.gz \
            1>>${wes_dir}${file}/${tumor_file}_mutect.log 2>&1
        fi
        if [ ! -f ${wes_dir}${file}/${tumor_file}_mutect2_somatic_filtered.vcf ]
        then
            $GATK --java-options "-Xmx20G -Djava.io.tmpdir=${tmp_dir}" FilterMutectCalls \
            -R ${ref} \
            -V ${wes_dir}${file}/${tumor_file}_mutect2.vcf \
            -O ${wes_dir}${file}/${tumor_file}_mutect2_somatic.vcf \
            --tumor-segmentation  ${wes_dir}${file}/${tumor_file}_segments.table      \
            --contamination-table ${wes_dir}${file}/${tumor_file}_contamination.table \
            --ob-priors ${wes_dir}${file}/${tumor_file}_read-orientation-model.tar.gz \
            1>>${wes_dir}${file}/${tumor_file}_mutect.log 2>&1

            ${GATK} --java-options "-Xmx20G -Djava.io.tmpdir=${tmp_dir}"  SelectVariants \
            -R ${ref}                     \
            -V ${wes_dir}${file}/${tumor_file}_mutect2_somatic.vcf \
            --exclude-filtered               \
            -O ${wes_dir}${file}/${tumor_file}_mutect2_somatic_filtered.vcf    \
            1>>${wes_dir}${file}/${tumor_file}_mutect.log 2>&1
        fi
        cat ${wes_dir}${file}/${tumor_file}_mutect2_somatic.vcf | \
        awk 'BEGIN{FS="\t";OFS="\t"} {if($0~"^#") {print $0} else {if($7!="PASS"){$7="REJECT"} print $0} }'  > ${wes_dir}${file}/${tumor_file}_mutect2_somatic_cnvkit.vcf
        sed -i '2 i ##FILTER=<ID=REJECT,Description="Not somatic due to normal call frequency or phred likelihoods: tumor: 35, normal 35.">' ${wes_dir}${file}/${tumor_file}_mutect2_somatic_cnvkit.vcf
        sed -i "2 i ##PEDIGREE=<Derived=${tumor_file},Original=${normal_file}>" ${wes_dir}${file}/${tumor_file}_mutect2_somatic_cnvkit.vcf
        if [ ! -f ${wes_dir}${file}/${normal_file}.cnn ]
        then
            cnvkit.py batch "${wes_dir}${file}/${tumor_file}.bam"  \
            --normal  "${wes_dir}${file}/${normal_file}.bam" \
            --targets  "/home/data/t070331/WES/LINXIN/S31285117_Regions.bed" \
            --fasta "/home/data/t070331/WES/LINXIN/ucsc.hg19.fasta"  \
            --annotate "/home/data/t070331/WES/LINXIN/refFlat.txt"  \
            --scatter --diagram --access "/home/data/t070331/WES/LINXIN/access-5k-mappable.hg19.bed" \
            --output-reference "${wes_dir}${file}/${normal_file}.cnn" --drop-low-coverage -p 30 --output-dir ./cnvkit
        else
            cnvkit.py batch "${wes_dir}${file}/${tumor_file}.bam" -r "${wes_dir}${file}/${normal_file}.cnn" --drop-low-coverage --scatter --diagram --method amplicon -p 8 --output-dir ./cnvkit
        fi
        cd cnvkit
        cnvkit.py export seg ${tumor_file}.cns -o "${wes_dir}${file}/${tumor_file}_${normal_file}.seg"
        cnvkit.py export seg ${tumor_file}.call.cns -o "${wes_dir}${file}/${tumor_file}_${normal_file}.call.seg"
        cnvkit.py export bed ${tumor_file}.call.cns -o "${wes_dir}${file}/${tumor_file}_${normal_file}.bed"
        cnvkit.py call ${tumor_file}.cns -v ${wes_dir}${file}/${tumor_file}_mutect2_somatic_cnvkit.vcf -o ${tumor_file}_${normal_file}_cn2_call.cns --sample-id  ${tumor_file} --normal-id ${normal_file}
        cd "${wes_dir}${file}"
        # sed 's/_bqsr//' gistic.segments
        # awk '{print FILENAME"\t"$0}' *.cns  | grep -v chromosome |sed 's/.cns//g' |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8"\t"$6}' > "${wes_dir}${file}/${tumor_file}.seg"
        genozip --reference ~/WES/LINXIN/ucsc.hg19.ref.genozip --best -@ 12 --replace "${wes_dir}${file}/${tumor_file}.bam"
        genozip --reference ~/WES/LINXIN/ucsc.hg19.ref.genozip --best -@ 12 --force "${wes_dir}${file}/${tumor_file}_bqsr.bam"
        mv ${wes_dir}${file}/${tumor_file}_bqsr.bam.genozip ~/WES/LINXIN/Results/${normal_file}/${tumor_file}_bqsr.bam.genozip
        mv ${wes_dir}${file}/${tumor_file}_${normal_file}.seg ~/WES/LINXIN/Results/${normal_file}/${tumor_file}_${normal_file}.seg
        mv ${wes_dir}${file}/${tumor_file}_${normal_file}.call.seg ~/WES/LINXIN/Results/${normal_file}/${tumor_file}_${normal_file}.call.seg
        mv ${wes_dir}${file}/${tumor_file}_${normal_file}.bed ~/WES/LINXIN/Results/${normal_file}/${tumor_file}_${normal_file}.bed
        mv ${wes_dir}${file}/${tumor_file}.bam.genozip ~/WES/LINXIN/Results/${normal_file}/${tumor_file}.bam.genozip
        mv ${wes_dir}${file}/${tumor_file}/cnvkit ~/WES/LINXIN/Results/${normal_file}/${tumor_file}_${normal_file}_cnvkit
        mv ${wes_dir}${file}/${tumor_file}/${tumor_file}.html ~/WES/LINXIN/Results/${normal_file}/${tumor_file}.html
        mv ${wes_dir}${file}/${tumor_file}/${tumor_file}.json ~/WES/LINXIN/Results/${normal_file}/${tumor_file}.json
        mv ${wes_dir}${file}/${tumor_file}_mutect2.vcf ~/WES/LINXIN/Results/${normal_file}/${tumor_file}_mutect2.vcf
        mv ${wes_dir}${file}/${tumor_file}_mutect2_somatic_filtered.vcf ~/WES/LINXIN/Results/${normal_file}/${tumor_file}_mutect2_somatic_filtered.vcf
        mv ${wes_dir}${file}/${tumor_file}_mutect2_somatic.vcf ~/WES/LINXIN/Results/${normal_file}/${tumor_file}_mutect2_somatic.vcf
        mv `ls|grep ".vcf"` ~/WES/LINXIN/Results/${normal_file}/
        cd ${wes_dir}${file}/${tumor_file}
        if [ -f ~/WES/LINXIN/Results/${normal_file}/${tumor_file}_bqsr.bam.genozip ] && [ -f ~/WES/LINXIN/Results/${normal_file}/${tumor_file}.bam.genozip ]
        then
            cd "${wes_dir}${file}"
            rm ${tumor_file}.bam
            rm ${wes_dir}${file}/${tumor_file}_bqsr.bam ${wes_dir}${file}/${tumor_file}_marked.bam
        fi
    done
    genozip --reference ~/WES/LINXIN/ucsc.hg19.ref.genozip --best -@ 12 --replace "${wes_dir}${file}/${normal_file}.bam"
    if [ ! -f ~/WES/LINXIN/Results/${normal_file}/${normal_file}_bqsr.bam.genozip ]
    then
        cd ${wes_dir}${file}
        genozip --reference ~/WES/LINXIN/ucsc.hg19.ref.genozip --best -@ 12 --force "${wes_dir}${file}/${normal_file}_bqsr.bam"
        mv ${wes_dir}${file}/${normal_file}_bqsr.bam.genozip ~/WES/LINXIN/Results/${normal_file}/${normal_file}_bqsr.bam.genozip
    fi
    mv ${wes_dir}${file}/${normal_file}.bam.genozip ~/WES/LINXIN/Results/${normal_file}/${normal_file}.bam.genozip
    mv ${wes_dir}${file}/${normal_file}.cnn ~/WES/LINXIN/Results/${normal_file}/${normal_file}.cnn
    cp ${wes_dir}${file}/tumor.txt ~/WES/LINXIN/Results/${normal_file}/tumor.txt
    cp ${wes_dir}${file}/normal.txt ~/WES/LINXIN/Results/${normal_file}/normal.txt
    if [ -f ~/WES/LINXIN/Results/${normal_file}/${normal_file}_bqsr.bam.genozip ] && [ -f /home/data/t070331/WES/LINXIN/Results/MSIH/HLA_Prediction/${normal_file}/winners.hla.txt ]
    then
        rm ${wes_dir}${file}/${normal_file}.bam
        rm ${wes_dir}${file}/${normal_file}_bqsr.bam ${wes_dir}${file}/${normal_file}_marked.bam
    fi
cd $wes_dir
done
cd $wes_dir