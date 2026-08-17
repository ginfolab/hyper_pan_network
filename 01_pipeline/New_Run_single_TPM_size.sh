#!/bin/bash

# 遍历当前目录下所有以 "GSE" 开头的文件夹
for dir in GSE*/; do
    echo "Checking directory: $dir"

    # 只处理包含 single_Ecoli.sh 的文件夹
    if [ ! -f "$dir/single_Ecoli.sh" ]; then
        echo "Skipping $dir (no single_Ecoli.sh found)"
        continue
    fi

    # 初始化执行标志
    run_script=false

    # 检查是否存在 TPM.txt
    if [ -f "$dir/TPM.txt" ]; then
        # 获取 TPM.txt 的大小（以字节为单位）
        file_size=$(stat -c %s "$dir/TPM.txt")

        # 如果文件大小小于 1MB（即 1048576 字节），则执行脚本
        if [ "$file_size" -lt 1048576 ]; then
            echo "TPM.txt in $dir is smaller than 1MB ($((file_size / 1024)) KB), executing script."
            run_script=true
        else
            echo "Skipping $dir (TPM.txt is larger than or equal to 1MB: $((file_size / 1024)) KB)"
        fi
    else
        echo "No TPM.txt found in $dir, executing script."
        run_script=true
    fi

    # 仅当满足条件时执行脚本
    if [ "$run_script" = true ]; then
        # 切换到 GSE 目录
        cd "$dir" || { echo "Failed to enter $dir"; continue; }

        # 执行 single_Ecoli.sh
        echo "Executing $dir/single_Ecoli.sh"
        bash "single_Ecoli.sh"
        echo "$dir/single_Ecoli.sh completed."

        # 返回上级目录
        cd .. || { echo "Failed to return to parent directory"; exit 1; }
        
        echo "Finished processing: $dir"
    fi
done

echo "All scripts executed."