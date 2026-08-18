#!/bin/bash

# 激活 Conda 环境（如果已经在该环境，可注释掉）
source ~/anaconda3/bin/activate salmon

# **指定 GSE 目录**
GSE_DIR=$(pwd)  # 当前目录作为 GSE 目录
echo "Working in directory: $GSE_DIR"

# **创建结果文件夹**
mkdir -p ./fastq ./quants ./readcount ./TPM

# **读取 SRR 号列表**
SRR_LIST=( $(cat ./SRR_Acc_List.txt) )

# **循环处理每个样本**
for sample in "${SRR_LIST[@]}"; do
    echo "Processing sample: $sample"

    # 使用 prefetch 下载 SRA 文件
    prefetch ${sample}
    if [ $? -ne 0 ]; then
        echo "Error: Failed to download ${sample} with prefetch." >> error.log
        continue
    fi

    # 使用 fasterq-dump 转换 FASTQ 格式
    fasterq-dump ${sample} --split-files --outdir ./fastq --verbose
    if [ $? -ne 0 ]; then
        echo "Error: Failed to convert ${sample} to FASTQ format." >> error.log
        continue
    fi
    
    
       # 删除以样本名命名的文件夹（prefetch 下载后的文件夹）
    echo "Deleting folder: ${sample}"
    rm -rf ${sample}
    if [ $? -ne 0 ]; then
        echo "Error: Failed to delete folder ${sample}" >> error.log
        continue
    fi
    
# **运行 Salmon 进行定量分析**
    salmon quant -i /lustre/home/users/exr/pangenome/Ecoli/Ecoli_index/ -l A \
        -r ./fastq/${sample}*.fastq \
        -p 8 -o ./quants/${sample}_quant        
    if [ $? -ne 0 ]; then
        echo "Error: Salmon quantification failed for $sample"
        continue
    fi
    
    # **提取 NumReads 和 TPM**
    cat ./quants/${sample}_quant/quant.sf | cut -f 5 | awk -v samp=${sample} 'NR==1 {$1=samp}1' > ./readcount/${sample}_readcount.txt
    cat ./quants/${sample}_quant/quant.sf | cut -f 4 | awk -v samp=${sample} 'NR==1 {$1=samp}1' > ./TPM/${sample}_TPM.txt
    cat ./quants/${sample}_quant/quant.sf | cut -f 1 > info.txt

    # **合并结果**
    paste ./info.txt ./readcount/*.txt > readcount.txt
    paste ./info.txt ./TPM/*.txt > TPM.txt

    # **清理 SRA 文件**
    rm -rf ~/ncbi/public/sra/${sample}.sra
    rm -rf ~/ncbi/public/sra/${sample}.sra.cache

    echo "Processing complete for $sample."
done

echo "All samples processed."


# **删除临时文件夹**
rm -rf ./fastq ./quants ./readcount ./TPM
echo "Temporary files deleted."